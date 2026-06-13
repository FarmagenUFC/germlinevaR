# Package-private in-memory cache for the InterPro auto-fetch.
# Lives at file scope so successive `gvr_lollipop()` calls share it within
# one R session. Keyed by sprintf("%s|%s", gene, organism).
# Placed BEFORE the roxygen block so it does not capture the @export tag.
.gvr_domain_mem_cache <- new.env(parent = emptyenv())

# Internal: retrying HTTP GET helper. Some UniProt/InterPro requests come back
# with transient HTTP 502/503/504 under load, or stall during TLS handshake.
# Retries up to `tries` times on those statuses (and on `httr::GET` errors such
# as timeout / connection reset), with exponential backoff. Returns the final
# response object on success/error-after-retries; the caller still inspects
# `httr::status_code()` so non-retryable 4xx responses surface as before.
.gvr_http_get_retry <- function(url, timeout_s = 30, tries = 3) {
  delays <- c(0.5, 1.5, 3.0)  # seconds between attempts
  last_err <- NULL
  for (i in seq_len(tries)) {
    resp <- tryCatch(
      httr::GET(url, httr::timeout(timeout_s)),
      error = function(e) { last_err <<- e; NULL }
    )
    if (!is.null(resp)) {
      sc <- httr::status_code(resp)
      if (!(sc %in% c(429L, 500L, 502L, 503L, 504L))) return(resp)
      # transient HTTP status: fall through to retry
    }
    if (i < tries) Sys.sleep(delays[i])
  }
  if (!is.null(resp)) return(resp)
  # All attempts errored out; re-throw the last error for the caller's tryCatch
  stop(last_err)
}

# Internal: abbreviate verbose InterPro domain names for plot readability.
# Applies targeted substitutions (NOT generic truncation) so the meaning is
# preserved. Returns input unchanged when `abbrev = FALSE`.
.gvr_abbrev_domain <- function(name, abbrev = TRUE) {
  if (!isTRUE(abbrev) || length(name) == 0L) return(name)
  # Apply substitutions in order. Each pair is c(pattern, replacement).
  rules <- list(
    c("von Willebrand factor",                    "VWF"),
    c("Cellular tumor antigen p53",               "p53"),
    c("VWF/SSPO/Zonadhesin-like",                  "VWF/SSPO/Zon"),
    c("cysteine-rich domain",                      "CR domain"),
    c("cysteine rich domain",                      "CR domain"),
    c("Trypsin Inhibitor-like",                    "Trypsin-Inh"),
    c("Otogelin-like/Mucin, TIL domain",           "Otogelin-like TIL"),
    c(", C-terminal",                              " C-term"),
    c(", N-terminal",                              " N-term"),
    c("transactivation domain",                    "TAD"),
    c("DNA-binding domain",                        "DBD"),
    c("tetramerisation domain",                    "tetram."),
    c("tetramerization domain",                    "tetram."),
    c(", type D domain",                           " type-D"),
    c("Cystine knot",                              "Cys-knot")
  )
  out <- name
  for (r in rules) out <- gsub(r[[1]], r[[2]], out, fixed = TRUE)
  out
}

#' Per-gene amino-acid lollipop plot for a germline MAF
#'
#' @description
#' Builds a ggplot2 lollipop plot of every protein-altering variant of a single
#' gene in a `read.gvr()` / `gvr_filter()` / `gvr_novel()` MAF. Each amino-acid
#' position is a stem; each sample-variant carried at that position is one dot
#' stacked on the stem; dots are coloured by `Variant_Classification` using the
#' package's colourblind-safe palette. A horizontal protein-body bar is drawn
#' along the x-axis (`maftools::lollipopPlot` style); protein domains are
#' overlaid as coloured rectangles on top of that bar. By default, domains are
#' fetched automatically from the EBI InterPro REST API (`domains = "auto"`).
#' Pass `domains = NULL` for a plain bar with no domains, or supply a
#' `data.frame` for custom domain coordinates. By default writes both an
#' `.svg` and a `.png` to disk (matching [gvr_plot()] behaviour).
#'
#' @details
#' Variant classes are filtered down to the canonical protein-altering set (the
#' same 9 classes used by [gvr_filter()] when `vc_nonSyn = TRUE`):
#' `Frame_Shift_Del`, `Frame_Shift_Ins`, `Splice_Site`, `Translation_Start_Site`,
#' `Nonsense_Mutation`, `Nonstop_Mutation`, `In_Frame_Del`, `In_Frame_Ins`,
#' `Missense_Mutation`. The default can be overridden with the `vc_keep`
#' argument.
#'
#' Amino-acid positions are extracted from `HGVSp_Short` using the regex
#' `^p\\.[A-Z*]([0-9]+)`. Rows whose `HGVSp_Short` does not match (complex
#' indels with `p.M1?`, splice variants without `p.`, etc.) are dropped with a
#' warning showing how many rows were dropped. Multi-amino-acid descriptors
#' (e.g. `p.D1850_E1852del`) match the first amino-acid number, which is the
#' standard 5'-most position convention.
#'
#' Protein length (the x-axis upper bound) is derived from the
#' `Protein_position` column: each value is a string `"<current>/<total>"`. The
#' most-common `total` across the gene's rows is used. If no row has a
#' parseable `total`, the x-axis falls back to `max(position) * 1.1`. The
#' `protein_length` argument overrides this.
#'
#' Stacking: each surviving sample-variant becomes one dot. Dots sharing the
#' same amino-acid position are stacked vertically, so the height of a stack
#' equals the count of sample-variants at that position. Stems run from the
#' top of the protein-body bar up to the top of the stack.
#'
#' Labels: the top-`label_top` positions by count are annotated with their
#' `HGVSp_Short` (e.g. `p.R175H`). When multiple distinct `HGVSp_Short` share
#' a position, the most-common one is shown followed by `(+N more)`. Set
#' `label_top = 0` to disable labels, `label_top = Inf` to label every
#' position. If \pkg{ggrepel} is installed, labels are placed via
#' `ggrepel::geom_text_repel()`; otherwise a plain `geom_text()` is used.
#'
#' Protein-body bar: a thin horizontal rectangle is drawn straddling `y = 0`
#' from x = 0 to x = `protein_length`. Bar height scales with the visible
#' stack height (about 4 percent of the y range, clamped to `[0.15, 0.5]`).
#' The y-axis is expanded slightly below zero to accommodate the bar.
#'
#' Domains: `domains` accepts three forms:
#' \itemize{
#'   \item `"auto"` (default) - the function calls the EBI InterPro REST API
#'         for the `gene` (gene -> UniProt accession via `rest.uniprot.org`,
#'         then UniProt -> InterPro domains via `www.ebi.ac.uk/interpro`),
#'         keeps only InterPro-integrated entries (those with
#'         `source_database == "interpro"`, e.g. `IPR011615`), caches the
#'         result, and converts to the same `data.frame` shape as below.
#'         Requires the optional packages \pkg{httr} and \pkg{jsonlite}; if
#'         either is missing, or any HTTP / parsing step fails, a warning is
#'         emitted and the plot still renders with the plain bar.
#'   \item `NULL` - no domain rectangles, just the plain bar.
#'   \item A `data.frame` with required columns `start`, `end` and optional
#'         columns `name` (label, shown when the rectangle is wide enough) and
#'         `color` (hex code or R colour name; auto-assigned from a
#'         colourblind-safe palette when missing). Rows entirely outside
#'         `[0, protein_length]` are dropped; partial overlaps are clipped.
#' }
#'
#' Cache (used when `domains = "auto"`, the default): the resolved cache
#' directory is picked from the first writable entry in this precedence chain:
#' \enumerate{
#'   \item `cache_dir` argument (explicit override; pass `FALSE` to disable
#'         on-disk caching entirely).
#'   \item Environment variable `GVR_CACHE_DIR` (useful for HPC `/scratch/$USER`).
#'   \item R option `getOption("germlinevaR.cache_dir")` (for `.Rprofile`).
#'   \item `tools::R_user_dir("germlinevaR", "cache")` (XDG-compliant default).
#'   \item `file.path(tempdir(), "germlinevaR_cache")` (last-resort, per-session).
#' }
#' Cached files: `<cache_dir>/domains_interpro_<GENE>_<ORG>.rds`. To force a
#' refresh, delete the file (or use the [gvr_domain_cache_clear()] helper).
#'
#' Empty-gene behaviour: if no row matches the gene, or none survives the
#' `vc_keep` filter, or none has a parseable amino-acid position, the function
#' issues a warning and returns a ggplot2 object with a single centered
#' "No protein-altering variants for <GENE>" annotation. No files are written
#' for empty-gene plots.
#'
#' File output: by default the function writes both SVG and PNG to `out_dir`
#' (default `"."`, the current working directory), using `out_prefix` (default
#' `gene`) as the filename prefix:
#' \itemize{
#'   \item `<out_dir>/<out_prefix>_lollipop.svg` - vector (lossless)
#'   \item `<out_dir>/<out_prefix>_lollipop.png` - raster at `dpi`
#' }
#' Both files are first rendered to `tempdir()` then `cp`'d to `out_dir`, so
#' S3-backed FUSE mounts (e.g. `/mnt/results/`) work without 0-byte issues.
#' Set `out_dir = NULL` to suppress file output and return only the ggplot2
#' object. The function always returns the ggplot2 object regardless of whether
#' files were written.
#'
#' @param maf A `data.table` / `data.frame` MAF from [read.gvr()], or any
#'   compatible table with at least `Hugo_Symbol`, `HGVSp_Short`,
#'   `Variant_Classification`, `Tumor_Sample_Barcode`, `Protein_position`.
#' @param gene Character(1). The gene symbol (matched against `Hugo_Symbol`).
#' @param vc_keep Character vector of `Variant_Classification` values to keep.
#'   `NULL` (default) uses the 9-class protein-altering set from
#'   [gvr_filter()].
#' @param protein_length Integer(1) or `NULL`. If `NULL` (default), derived
#'   from the `Protein_position` column (most-common total). Provide explicitly
#'   to fix the x-axis length.
#' @param domains `"auto"`, `NULL`, or a `data.frame`. Default `"auto"`:
#'   domains are fetched from the EBI InterPro REST API and cached. `NULL`
#'   draws a plain bar with no domain rectangles. A `data.frame` with columns
#'   `start`, `end` (and optional `name`, `color`) draws custom domain
#'   rectangles. See \strong{Details} for the three accepted forms.
#' @param organism Integer or string. NCBI taxonomy id used when
#'   `domains = "auto"` (passed through to the UniProt search as
#'   `organism_id`). Default `9606L` (human). Ignored when `domains` is not
#'   `"auto"`.
#' @param cache_dir `NULL`, a directory path, or `FALSE`. Controls the
#'   on-disk cache used for `domains = "auto"`. `NULL` (default) triggers the
#'   precedence chain documented in \strong{Details}. `FALSE` disables on-disk
#'   caching. Ignored when `domains` is not `"auto"`.
#' @param label_top Integer(1). Number of top-counted positions to label.
#'   `0` disables labels, `Inf` labels every position. Default `5L`.
#' @param domain_label_min_frac Numeric(1) in `[0, 1]`. Minimum domain width
#'   as a fraction of `protein_length` for a domain label to be rendered in
#'   `"name"` and `"id"` modes. Default `0.05` (5%). Lower this for long
#'   proteins (e.g. `0.01` so that even small domains get labelled). Ignored
#'   when `domain_label_mode = "number"` (all domains labelled) or `"none"`.
#' @param domain_label_mode One of `"name"` (default), `"id"`, `"number"`,
#'   or `"none"`. Controls how domain rectangles are labelled.
#'   * `"name"`: human-readable InterPro name, repelled below the bar with
#'     leader lines (uses [ggrepel::geom_text_repel()] when available).
#'   * `"id"`: InterPro accession (e.g. `IPR011615`); falls back to name
#'     when accession is missing (user-supplied data.frame without it).
#'   * `"number"`: numbers 1..N centered in each rectangle, with a
#'     companion `"Domains"` legend mapping number to full name. Requires
#'     the optional `ggnewscale` package; without it, only numbers render.
#'   * `"auto"`: opt-in heuristic. Use `"number"` when
#'     `protein_length > 2000` aa AND `>= 5` domains are drawn
#'     (long, densely-annotated proteins where in-plot names overlap),
#'     otherwise use `"name"`. With `verbose = TRUE`, the resolved mode
#'     and reason are printed.
#'   * `"none"`: no labels (rely on the variant legend only).
#' @param domain_name_abbrev Logical(1). When `TRUE` (default), apply a small
#'   set of substitutions to compress verbose InterPro domain names (e.g.
#'   "von Willebrand factor, type D domain" becomes "VWF type-D") so that
#'   labels fit inside narrower rectangles in `domain_label_mode = "name"`.
#'   Set `FALSE` to keep the raw InterPro names verbatim.
#' @param domain_label_position One of `"inside"` (default) or `"below"`.
#'   Controls where InterPro domain labels render relative to the domain
#'   rectangle. `"inside"` places labels centred inside each rectangle with
#'   an automatically chosen black/white text colour for contrast against the
#'   domain fill (WCAG luminance rule). `"below"` reproduces the legacy
#'   layout used in earlier versions (labels under the bar, leader lines from
#'   each rectangle).
#'
#'   When `domain_label_position = "inside"`, labels that don't fit inside
#'   their domain rectangle (estimated by character count vs rectangle width)
#'   fall back to the `"below"` style with a leader line on a per-domain
#'   basis so no information is lost. A verbose message reports the count of
#'   overflowing labels.
#' @param hotspot_window Integer(1). Sliding-window width (in amino
#'   acids) used to detect mutation hotspots. A hotspot is a region
#'   containing at least `hotspot_min_n` distinct variant positions
#'   within `hotspot_window` aa. Drawn as a soft translucent vertical
#'   band behind the bar/domains/stems. Default `20L` (publication-tight
#'   rule). Increase for noisier/exploratory cohorts.
#' @param hotspot_min_n Numeric(1). Minimum distinct variant positions
#'   inside a `hotspot_window`-wide region for it to be drawn as a
#'   hotspot band. Default `4`. Pass `Inf` to disable hotspot detection
#'   entirely. Counting uses unique aa positions (not sample counts), so
#'   a single recurrently-hit position does not by itself create a band.
#' @param stem_alpha Numeric(1) in `[0, 1]`. Opacity of lollipop stems.
#'   Default `0.6`. Lower this further for very dense plots; raise to `1`
#'   for the original solid stems.
#' @param point_size Numeric(1). Dot size. Default `3`.
#' @param stem_color Character(1). Stem (vertical line) colour. Default
#'   `"grey50"`.
#' @param bar_color Character(1). Protein-body bar fill colour. Default
#'   `"grey85"`.
#' @param bar_border Character(1) or `NA`. Protein-body bar border colour.
#'   Default `"grey40"`.
#' @param base_size Numeric(1). ggplot2 base font size. Default `12`.
#' @param axis_text_size Numeric(1). Font size of axis tick labels (aa
#'   positions on x, count on y). Default `11`, matching the existing visual
#'   output. Set to `12` or higher for publication-style figures.
#' @param axis_title_size Numeric(1). Font size of axis titles ("Amino acid
#'   position" / "Variant count"). Default `12`. Set to `14` for sibling-
#'   package publication defaults.
#' @param axis_text_color Character(1). Colour of axis tick labels. Default
#'   `"grey20"` (current visual). Set to `"black"` for sibling-package
#'   publication defaults.
#' @param axis_title_color Character(1). Colour of axis titles. Default
#'   `"black"`.
#' @param axis_line_color Character(1). Colour of axis lines and ticks.
#'   Default `"grey40"` (current visual). Set to `"black"` for sibling-
#'   package publication defaults.
#' @param axis_line_width Numeric(1). Line width (in `ggplot2` linewidth
#'   units) of axis lines and ticks. Default `0.4` (current visual). Set to
#'   `0.8` for sibling-package publication defaults.
#' @param out_dir Character(1) or `NULL`. Output directory for the SVG/PNG
#'   files. Created if it does not exist. Default `"."` (current working
#'   directory), matching [gvr_plot()]. Set to `NULL` to suppress file output.
#' @param out_prefix Character(1) or `NULL`. Filename prefix; output files are
#'   `<out_prefix>_lollipop.svg` / `<out_prefix>_lollipop.png`. Default `gene`
#'   (the gene symbol), so `gvr_lollipop(f, "TP53")` writes
#'   `TP53_lollipop.svg` and `TP53_lollipop.png`. Set to `NULL` to suppress
#'   file output.
#' @param width Numeric(1). Plot width in inches. Default `10`.
#' @param height Numeric(1). Plot height in inches. Default `4`.
#' @param dpi Numeric(1). PNG resolution. Default `300`.
#' @param variant_palette Character. Variant colour assignment. Accepts three
#'   forms. (1) A palette name (e.g. `"gvr"`, `"nature"`, `"set2"`); the
#'   default `"gvr"` is the semantic palette where, e.g., `Missense_Mutation`
#'   is always green and `Nonsense_Mutation` always orange regardless of which
#'   classes are present. Any other palette name is treated as positional
#'   ordinal (colours assigned in the order classes appear in the data).
#'   (2) A named character vector of explicit hex codes for specific
#'   `Variant_Classification` values, e.g.
#'   `c(Missense_Mutation = "#FF0000", Nonsense_Mutation = "#000000")`;
#'   classes not listed inherit the `"gvr"` defaults. (3) A mixed vector
#'   with one unnamed element used as the fallback palette name, e.g.
#'   `c(Missense_Mutation = "#FF00FF", "nature")` (Missense magenta, all
#'   other classes from the `nature` palette). Use [gvr_list_palettes()] to
#'   discover names and [gvr_color_palette()] to inspect colours.
#' @param domain_palette Character. Palette used for InterPro domain
#'   rectangle fills. Accepts the same three forms as `variant_palette`:
#'   palette name (default `"okabe_ito"`), named vector keyed by accession
#'   or domain name, or mixed vector with fallback palette name. Unlike
#'   `variant_palette` there is no `"gvr"` semantic special case; the
#'   default cycles a 9-colour colour-blind-safe palette across all domains.
#' @param verbose Logical(1). If `TRUE` (default), emit progress messages
#'   (counts, dropped rows, file paths, cache hits, network calls).
#'
#' @return A ggplot2 object. By default SVG and PNG files are also written to
#'   disk as a side effect. Set `out_dir = NULL` to suppress file output.
#'
#' @seealso [read.gvr()], [gvr_filter()], [gvr_novel()], [gvr_summary()],
#'   [gvr_plot()], [gvr_domain_cache_clear()]
#'
#' @examples
#' \dontrun{
#'   maf <- read.gvr("vcf_dir/", pattern = "\\.vep\\.vcf\\.gz$")
#'   f   <- gvr_filter(maf)
#'
#'   ## default: auto-fetches InterPro domains, auto-saves to "."
#'   p <- gvr_lollipop(f, "TP53")
#'
#'   ## second call uses cache, no network
#'   gvr_lollipop(f, "TP53")
#'
#'   ## plain bar, no domain rectangles
#'   gvr_lollipop(f, "TP53", domains = NULL)
#'
#'   ## with user-supplied domains (TP53 example - canonical UniProt P04637)
#'   tp53_dom <- data.frame(
#'     start = c(95,   323),
#'     end   = c(288,  356),
#'     name  = c("DNA-binding", "Tetramerization"),
#'     color = c("#FF9400", "#75A025")
#'   )
#'   gvr_lollipop(f, "TP53", domains = tp53_dom)
#'
#'   ## non-human cohort: pass the NCBI taxonomy id
#'   gvr_lollipop(f, "Trp53", organism = 10090) # mouse
#'
#'   ## HPC: send the cache to scratch
#'   Sys.setenv(GVR_CACHE_DIR = "/scratch/$USER/germlinevaR_cache")
#'   gvr_lollipop(f, "MUC16")
#'
#'   ## disable on-disk caching entirely (CI / shared scratch)
#'   gvr_lollipop(f, "BRCA1", cache_dir = FALSE)
#'
#'   ## suppress file output, return ggplot2 only
#'   p <- gvr_lollipop(f, "BRCA1", out_dir = NULL)
#'   print(p)
#'
#'   ## Long protein, two label strategies side-by-side
#'   ## MUC19 has 9-12 InterPro domains clustered in aa 470-1604, so
#'   ## name-mode labels overlap. Number-mode keeps the in-plot annotation
#'   ## compact and routes domain names into a side legend.
#'   gvr_lollipop(f, "MUC19", domain_label_mode = "name")
#'   gvr_lollipop(f, "MUC19", domain_label_mode = "number")
#'
#'   ## Opt-in heuristic chooses between "name" and "number" automatically
#'   ## (uses "number" when protein > 2000 aa AND >= 5 InterPro domains).
#'   gvr_lollipop(f, "MUC19", domain_label_mode = "auto")
#'
#'   ## Highlight hotspots: clusters of >= 4 variants within 20 aa
#'   ## (the default; set hotspot_min_n = Inf to disable).
#'   gvr_lollipop(f, "TP53", hotspot_window = 20, hotspot_min_n = 4)
#'
#'   ## Inside-rectangle domain labels (new default in Phase J)
#'   gvr_lollipop(f, "TP53", domain_label_position = "inside")
#'
#'   ## Legacy below-rectangle labels with leader lines
#'   gvr_lollipop(f, "TP53", domain_label_position = "below")
#'
#'   ## Sibling-package publication defaults for axes
#'   gvr_lollipop(f, "TP53",
#'                axis_text_size = 12, axis_title_size = 14,
#'                axis_text_color = "black", axis_title_color = "black",
#'                axis_line_color = "black", axis_line_width = 0.8)
#'
#'   ## Colour-blind-safe variant palette (positional, not semantic)
#'   gvr_lollipop(f, "TP53", variant_palette = "okabe_ito")
#'
#'   ## Override Missense colour only; rest stay GVR semantic
#'   gvr_lollipop(f, "TP53",
#'                variant_palette = c(Missense_Mutation = "#FF00FF"))
#'
#'   ## Override Missense + use Nature palette for the remaining classes
#'   gvr_lollipop(f, "TP53",
#'                variant_palette = c(Missense_Mutation = "#FF00FF", "nature"))
#'
#'   ## Domain rectangles coloured along the viridis gradient
#'   gvr_lollipop(f, "MUC19", domain_palette = "viridis")
#'
#'   ## List available palettes / inspect colours
#'   gvr_list_palettes()
#'   gvr_color_palette("nature", n = 5)
#' }
#'
#' @importFrom data.table as.data.table copy is.data.table setDT :=
#' @importFrom ggplot2 ggplot aes geom_segment geom_point geom_rect geom_text
#'   annotate scale_color_manual scale_fill_identity scale_x_continuous
#'   scale_y_continuous labs theme_classic theme element_text element_blank
#'   ggsave
#' @importFrom grDevices svg dev.off
#' @export
gvr_lollipop <- function(maf, gene,
                         vc_keep        = NULL,
                         protein_length = NULL,
                         domains        = "auto",
                         organism       = 9606L,
                         cache_dir      = NULL,
                         label_top      = 5L,
                         domain_label_min_frac = 0.05,
                         domain_label_mode     = c("name", "id", "number", "none", "auto"),
                         domain_name_abbrev    = TRUE,
                         domain_label_position = c("inside", "below"),
                         hotspot_window        = 20L,
                         hotspot_min_n         = 4,
                         stem_alpha            = 0.6,
                         point_size     = 3,
                         stem_color     = "grey50",
                         bar_color      = "grey85",
                         bar_border     = "grey40",
                         base_size      = 12,
                         axis_text_size   = 11,
                         axis_title_size  = 12,
                         axis_text_color  = "grey20",
                         axis_title_color = "black",
                         axis_line_color  = "grey40",
                         axis_line_width  = 0.4,
                         out_dir        = ".",
                         out_prefix     = gene,
                         width          = 10,
                         height         = 4,
                         dpi            = 300,
                         variant_palette = "gvr",
                         domain_palette  = "okabe_ito",
                         verbose        = TRUE) {

  # ---- Nested constants (byte-identical copies from gvr_filter.R / gvr_plot.R) ----
  .vc_nonSyn_default <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site",
                           "Translation_Start_Site", "Nonsense_Mutation",
                           "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins",
                           "Missense_Mutation")

  GVR_CLASS_COLORS <- c(
    "Translation_Start_Site" = "#000000", "Nonsense_Mutation" = "#D55E00",
    "Nonstop_Mutation" = "#882255", "Splice_Site" = "#CC79A7",
    "Frame_Shift_Del" = "#E69F00", "Frame_Shift_Ins" = "#F0E442",
    "In_Frame_Del" = "#56B4E9", "In_Frame_Ins" = "#0072B2",
    "Missense_Mutation" = "#009E73", "Splice_Region" = "#44AA99",
    "Protein_altering_variant" = "#117733", "Silent" = "#999933",
    "5'UTR" = "#AA4499", "3'UTR" = "#DDCC77", "5'Flank" = "#88CCEE",
    "3'Flank" = "#332288", "RNA" = "#BBBBBB", "Intron" = "#DDDDDD",
    "IGR" = "#777777", "Targeted_Region" = "#666666", "Other" = "#CCCCCC")

  # colourblind-safe palette for auto-assigned domain colours
  .GVR_DOMAIN_PALETTE <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                           "#0072B2", "#D55E00", "#CC79A7", "#999999",
                           "#88CCEE", "#44AA99", "#117733", "#DDCC77",
                           "#AA4499")


  # ---- Palette resolvers (variant + domain) ---------------------------
  # Accept three input forms:
  #   (a) single palette name ("gvr" / one of gvr_list_palettes())
  #   (b) named character vector of fixed overrides; rest = GVR default
  #   (c) mixed: named overrides + ONE unnamed palette-name fallback

  .resolve_variant_palette <- function(vp, classes_present) {
    classes_present <- as.character(classes_present)
    if (length(classes_present) == 0L) return(character(0))

    # Case (a): single character string with no names
    if (is.character(vp) && length(vp) == 1L && is.null(names(vp))) {
      if (identical(vp, "gvr")) {
        out <- GVR_CLASS_COLORS[classes_present]
        out[is.na(out)] <- GVR_CLASS_COLORS[["Other"]]
        names(out) <- classes_present
        return(out)
      }
      cols <- gvr_color_palette(vp, length(classes_present))
      return(setNames(cols, classes_present))
    }

    # Case (b) / (c): named vector (possibly with one unnamed fallback)
    if (!is.character(vp))
      stop("`variant_palette` must be a character string or character vector")
    nm <- names(vp); if (is.null(nm)) nm <- rep("", length(vp))
    overrides <- vp[nzchar(nm)]
    fallbacks <- vp[!nzchar(nm)]
    if (length(fallbacks) > 1L)
      stop("`variant_palette` accepts at most one unnamed element ",
           "(used as the fallback palette name); got ", length(fallbacks))
    fallback_name <- if (length(fallbacks) == 1L) unname(fallbacks) else "gvr"

    # Build base from fallback palette
    out <- if (identical(fallback_name, "gvr")) {
      tmp <- GVR_CLASS_COLORS[classes_present]
      tmp[is.na(tmp)] <- GVR_CLASS_COLORS[["Other"]]
      names(tmp) <- classes_present; tmp
    } else {
      setNames(gvr_color_palette(fallback_name, length(classes_present)),
               classes_present)
    }
    # Apply explicit overrides
    keep <- intersect(names(overrides), classes_present)
    if (length(keep)) out[keep] <- overrides[keep]
    out
  }

  .resolve_domain_palette <- function(dp, dom_df) {
    n_d <- if (is.null(dom_df)) 0L else nrow(dom_df)
    if (n_d == 0L) return(character(0))

    # Case (a): single palette name
    if (is.character(dp) && length(dp) == 1L && is.null(names(dp))) {
      return(gvr_color_palette(dp, n_d))
    }

    if (!is.character(dp))
      stop("`domain_palette` must be a character string or character vector")
    nm <- names(dp); if (is.null(nm)) nm <- rep("", length(dp))
    overrides <- dp[nzchar(nm)]
    fallbacks <- dp[!nzchar(nm)]
    if (length(fallbacks) > 1L)
      stop("`domain_palette` accepts at most one unnamed element ",
           "(used as the fallback palette name); got ", length(fallbacks))
    fallback_name <- if (length(fallbacks) == 1L) unname(fallbacks)
                     else "okabe_ito"

    # Base palette over all domains
    out <- gvr_color_palette(fallback_name, n_d)

    # Apply overrides keyed by accession (preferred) or name
    if (length(overrides)) {
      acc <- as.character(dom_df$accession)
      nms <- as.character(dom_df$name)
      for (k in names(overrides)) {
        idx <- which(acc == k)
        if (length(idx) == 0L) idx <- which(nms == k)
        if (length(idx) > 0L) out[idx] <- overrides[[k]]
      }
    }
    out
  }

  # ---- Contrast helper for inside-domain labels ------------------------
  # Returns "black" or "white" per input hex, by WCAG relative luminance.
  .pick_contrast_color <- function(fill_hex) {
    rgb <- grDevices::col2rgb(fill_hex) / 255
    # Linearize per WCAG (approximation: skip the sRGB->linear conversion
    # since our threshold of 0.5 is robust to the small bias)
    lum <- 0.2126 * rgb["red", ] + 0.7152 * rgb["green", ] +
           0.0722 * rgb["blue", ]
    ifelse(lum > 0.5, "black", "white")
  }

  # ---- Nested helpers ----
  .is_missing <- function(x) is.na(x) | x == ""

  .empty_plot <- function(gene_label) {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5,
                        label = sprintf("No protein-altering variants for %s", gene_label),
                        size = base_size * 0.4) +
      ggplot2::scale_x_continuous(limits = c(0, 1)) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::theme_classic(base_size = base_size) +
      ggplot2::theme(axis.line   = ggplot2::element_blank(),
                     axis.ticks  = ggplot2::element_blank(),
                     axis.text   = ggplot2::element_blank(),
                     axis.title  = ggplot2::element_blank()) +
      ggplot2::labs(title = sprintf("%s (no variants)", gene_label))
  }

  # ---- FUSE-safe writer (mirrors gvr_plot.R's .fuse_save_png) ----
  .fuse_save <- function(final_path, draw_fun) {
    tmp <- file.path(tempdir(), basename(final_path))
    ok <- tryCatch({ draw_fun(tmp); file.exists(tmp) && file.info(tmp)$size > 0 },
                   error = function(e) {
                     warning(sprintf("gvr_lollipop: render failed (%s): %s",
                                     basename(final_path), conditionMessage(e))); FALSE })
    if (!ok) return(NA_character_)
    system2("cp", c("-f", shQuote(tmp), shQuote(final_path)))
    if (!file.exists(final_path) || file.info(final_path)$size == 0) {
      warning(sprintf("gvr_lollipop: copy to '%s' failed; left at '%s'.", final_path, tmp))
      return(tmp)
    }
    final_path
  }

  # ---- Cache directory resolver (precedence chain) ----
  # Returns either a writable directory path, or NA_character_ if no disk
  # caching should be performed (cache_dir == FALSE).
  .gvr_resolve_cache_dir <- function(cache_dir_arg) {
    # Explicit FALSE disables on-disk caching
    if (isFALSE(cache_dir_arg)) return(NA_character_)

    candidates <- character(0)

    # 1. Explicit argument
    if (!is.null(cache_dir_arg) && is.character(cache_dir_arg) &&
        length(cache_dir_arg) == 1L && nzchar(cache_dir_arg)) {
      candidates <- c(candidates, cache_dir_arg)
    }

    # 2. Environment variable
    env_dir <- Sys.getenv("GVR_CACHE_DIR", unset = "")
    if (nzchar(env_dir)) candidates <- c(candidates, env_dir)

    # 3. R option
    opt_dir <- getOption("germlinevaR.cache_dir", default = NULL)
    if (!is.null(opt_dir) && is.character(opt_dir) && length(opt_dir) == 1L &&
        nzchar(opt_dir)) {
      candidates <- c(candidates, opt_dir)
    }

    # 4. XDG-compliant default via tools::R_user_dir (base R >= 4.0)
    candidates <- c(candidates, tools::R_user_dir("germlinevaR", which = "cache"))

    # 5. Last resort: tempdir
    candidates <- c(candidates, file.path(tempdir(), "germlinevaR_cache"))

    # Try each in order; return first that is creatable/writable
    for (cand in candidates) {
      ok <- tryCatch({
        if (!dir.exists(cand))
          dir.create(cand, recursive = TRUE, showWarnings = FALSE)
        # writability test: file.access mode 2 = write
        dir.exists(cand) && file.access(cand, mode = 2L) == 0L
      }, error = function(e) FALSE, warning = function(w) FALSE)
      if (isTRUE(ok)) return(cand)
    }
    NA_character_   # nothing worked
  }

  # tiny operator helper for null-coalesce (used inside auto-fetch)
  `%||%` <- function(a, b) if (is.null(a)) b else a

  # ---- InterPro auto-fetch (gene -> UniProt -> InterPro domains) ----
  # Returns a data.frame of (start, end, name, color, accession, source).
  # Any failure path returns a 0-row data.frame with a warning().
  .gvr_interpro_get <- function(gene_sym, org, cache_dir_arg, verbose) {
    empty_df <- data.frame(
      start = integer(0), end = integer(0),
      name = character(0), color = character(0),
      accession = character(0), source = character(0),
      stringsAsFactors = FALSE
    )

    # Soft-fail on missing optional deps
    if (!requireNamespace("httr", quietly = TRUE) ||
        !requireNamespace("jsonlite", quietly = TRUE)) {
      warning("gvr_lollipop: domains='auto' requires packages 'httr' and ",
              "'jsonlite'; install with install.packages(c('httr','jsonlite')). ",
              "Falling back to plain bar.", call. = FALSE)
      return(empty_df)
    }

    # In-memory cache lookup (session-scoped)
    mem_key <- sprintf("%s|%s", gene_sym, as.character(org))
    if (exists(mem_key, envir = .gvr_domain_mem_cache, inherits = FALSE)) {
      if (isTRUE(verbose))
        message(sprintf("gvr_lollipop: cached domains for %s (org %s) (session memory)",
                        gene_sym, as.character(org)))
      return(get(mem_key, envir = .gvr_domain_mem_cache, inherits = FALSE))
    }

    # Resolve cache_dir (may be NA if disk caching is disabled or unwritable)
    cdir <- .gvr_resolve_cache_dir(cache_dir_arg)
    cache_file <- if (!is.na(cdir))
      file.path(cdir, sprintf("domains_interpro_%s_%s.rds", gene_sym, as.character(org)))
    else NA_character_

    # On-disk cache lookup
    if (!is.na(cache_file) && file.exists(cache_file)) {
      cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
      if (is.data.frame(cached)) {
        # Phase J: drop legacy color column so domain_palette resolves fresh.
        if ("color" %in% names(cached)) cached$color <- NULL
        if (isTRUE(verbose))
          message(sprintf("gvr_lollipop: cached domains for %s (org %s) from %s",
                          gene_sym, as.character(org), cache_file))
        assign(mem_key, cached, envir = .gvr_domain_mem_cache)
        return(cached)
      }
      # Corrupt file: fall through and re-fetch
      if (isTRUE(verbose))
        message(sprintf("gvr_lollipop: cache file unreadable, re-fetching: %s",
                        cache_file))
    }

    # ---- Stage A: gene -> UniProt accession (single HTTP GET) ----
    if (isTRUE(verbose))
      message(sprintf("gvr_lollipop: fetching UniProt accession for %s (org %s) ...",
                      gene_sym, as.character(org)))

    uniprot_url <- sprintf(
      "https://rest.uniprot.org/uniprotkb/search?query=gene_exact:%s+AND+organism_id:%s+AND+reviewed:true&format=json&fields=accession,gene_names,length&size=5",
      utils::URLencode(gene_sym), utils::URLencode(as.character(org))
    )

    acc <- tryCatch({
      resp <- .gvr_http_get_retry(uniprot_url, timeout_s = 30, tries = 3)
      if (httr::status_code(resp) >= 400L) {
        warning(sprintf("gvr_lollipop: UniProt search returned HTTP %d for '%s'. Falling back to plain bar.",
                        httr::status_code(resp), gene_sym), call. = FALSE)
        return(empty_df)
      }
      txt <- httr::content(resp, as = "text", encoding = "UTF-8")
      js <- jsonlite::fromJSON(txt, simplifyVector = FALSE)
      if (length(js$results) == 0L) {
        warning(sprintf("gvr_lollipop: no reviewed UniProt entry for gene '%s' in organism %s. Falling back to plain bar.",
                        gene_sym, as.character(org)), call. = FALSE)
        return(empty_df)
      }
      if (length(js$results) > 1L && isTRUE(verbose)) {
        all_accs <- vapply(js$results, function(r) r$primaryAccession %||% NA_character_,
                            character(1))
        message(sprintf("gvr_lollipop: UniProt returned %d entries for '%s'; using first (%s). Others: %s",
                        length(js$results), gene_sym, all_accs[1],
                        paste(all_accs[-1], collapse = ", ")))
      }
      js$results[[1]]$primaryAccession
    }, error = function(e) {
      warning(sprintf("gvr_lollipop: UniProt fetch for '%s' failed: %s. Falling back to plain bar.",
                      gene_sym, conditionMessage(e)), call. = FALSE)
      NULL
    })
    if (is.null(acc) || !is.character(acc) || length(acc) == 0L) return(empty_df)
    # acc may be empty data.frame if Stage A returned via the warning() path
    if (is.data.frame(acc)) return(acc)

    # ---- Stage B: UniProt accession -> InterPro domains (single HTTP GET) ----
    if (isTRUE(verbose))
      message(sprintf("gvr_lollipop: fetching InterPro domains for %s (UniProt %s) ...",
                      gene_sym, acc))

    interpro_url <- sprintf(
      "https://www.ebi.ac.uk/interpro/api/entry/all/protein/UniProt/%s/?type=domain&page_size=100",
      utils::URLencode(acc)
    )

    ip <- tryCatch({
      resp <- .gvr_http_get_retry(interpro_url, timeout_s = 30, tries = 3)
      sc <- httr::status_code(resp)
      if (sc == 204L || sc == 404L) {
        warning(sprintf("gvr_lollipop: no InterPro domains for '%s' (UniProt %s)",
                        gene_sym, acc), call. = FALSE)
        return(empty_df)
      }
      if (sc >= 400L) {
        warning(sprintf("gvr_lollipop: InterPro returned HTTP %d for '%s'. Falling back to plain bar.",
                        sc, gene_sym), call. = FALSE)
        return(empty_df)
      }
      txt <- httr::content(resp, as = "text", encoding = "UTF-8")
      jsonlite::fromJSON(txt, simplifyVector = FALSE)
    }, error = function(e) {
      warning(sprintf("gvr_lollipop: InterPro fetch for '%s' failed: %s. Falling back to plain bar.",
                      gene_sym, conditionMessage(e)), call. = FALSE)
      NULL
    })
    if (is.null(ip)) return(empty_df)
    if (is.data.frame(ip)) return(ip)
    if (is.null(ip$results) || length(ip$results) == 0L) {
      warning(sprintf("gvr_lollipop: no InterPro domains for '%s' (UniProt %s)",
                      gene_sym, acc), call. = FALSE)
      return(empty_df)
    }

    # ---- Stage C: deduplicate + assemble ----
    rows <- list()
    for (entry in ip$results) {
      md <- entry$metadata
      if (is.null(md$source_database) || md$source_database != "interpro") next
      protein_locs <- entry$proteins[[1]]$entry_protein_locations
      if (is.null(protein_locs) || length(protein_locs) == 0L) next
      for (loc in protein_locs) {
        if (is.null(loc$fragments) || length(loc$fragments) == 0L) next
        for (frag in loc$fragments) {
          s <- suppressWarnings(as.integer(frag$start))
          e <- suppressWarnings(as.integer(frag$end))
          if (is.na(s) || is.na(e)) next
          rows[[length(rows) + 1L]] <- list(
            start     = s,
            end       = e,
            name      = md$name %||% "",
            accession = md$accession %||% NA_character_,
            source    = md$source_database
          )
        }
      }
    }
    if (length(rows) == 0L) {
      warning(sprintf("gvr_lollipop: no integrated InterPro domains for '%s' (UniProt %s)",
                      gene_sym, acc), call. = FALSE)
      return(empty_df)
    }

    df <- data.frame(
      start     = vapply(rows, function(r) r$start, integer(1)),
      end       = vapply(rows, function(r) r$end,   integer(1)),
      name      = vapply(rows, function(r) r$name,  character(1)),
      accession = vapply(rows, function(r) r$accession %||% NA_character_, character(1)),
      source    = vapply(rows, function(r) r$source,    character(1)),
      stringsAsFactors = FALSE
    )

    # Write cache (graceful on FS failure)
    if (!is.na(cache_file)) {
      ok <- tryCatch({ saveRDS(df, cache_file); TRUE },
                     error = function(e) {
                       warning(sprintf("gvr_lollipop: could not write cache file '%s': %s",
                                       cache_file, conditionMessage(e)), call. = FALSE)
                       FALSE })
      if (isTRUE(ok) && isTRUE(verbose))
        message(sprintf("gvr_lollipop: wrote cache %s", cache_file))
    }

    # Populate session cache
    assign(mem_key, df, envir = .gvr_domain_mem_cache)
    if (isTRUE(verbose))
      message(sprintf("gvr_lollipop: fetched %d InterPro domain(s) for %s (UniProt %s)",
                      nrow(df), gene_sym, acc))

    df
  }

  # ---- Input validation ----
  if (!is.data.frame(maf)) {
    stop("gvr_lollipop: 'maf' must be a data.frame / data.table.")
  }
  if (!is.character(gene) || length(gene) != 1L || !nzchar(gene)) {
    stop("gvr_lollipop: 'gene' must be a single non-empty string.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("gvr_lollipop requires the 'ggplot2' package.")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("gvr_lollipop requires the 'data.table' package.")
  }

  req_cols <- c("Hugo_Symbol", "HGVSp_Short", "Variant_Classification",
                "Tumor_Sample_Barcode", "Protein_position")
  miss_req <- req_cols[!req_cols %in% names(maf)]
  if (length(miss_req) > 0L) {
    stop(sprintf("gvr_lollipop: required column(s) not found: %s",
                 paste(miss_req, collapse = ", ")))
  }

  # ---- Auto-fetch dispatch (must happen BEFORE the existing data.frame validator) ----
  # If domains is the literal string "auto", call the InterPro fetcher and replace
  # domains with its result (a data.frame). Empty result -> fall through to the
  # NULL-equivalent path (plain bar, no domain rectangles).
  if (!is.null(domains) && is.character(domains) && length(domains) == 1L &&
      tolower(domains) == "auto") {
    domains <- .gvr_interpro_get(gene, organism, cache_dir, verbose)
    if (is.data.frame(domains) && nrow(domains) == 0L) {
      # Treat as no-domains: existing downstream code handles is.null(domains)
      domains <- NULL
    }
  }

  # ---- Resolve vc_keep ----
  vc_use <- if (is.null(vc_keep)) .vc_nonSyn_default else as.character(vc_keep)
  if (length(vc_use) == 0L) {
    stop("gvr_lollipop: 'vc_keep' must be NULL or a non-empty character vector.")
  }

  # ---- Subset to gene + protein-altering ----
  dt <- if (data.table::is.data.table(maf)) data.table::copy(maf) else data.table::as.data.table(maf)
  .gene_sym <- gene
  dt <- dt[Hugo_Symbol == .gene_sym]

  if (nrow(dt) == 0L) {
    warning(sprintf("gvr_lollipop: no rows for gene '%s' in input MAF", gene))
    return(.empty_plot(gene))
  }
  if (isTRUE(verbose)) message(sprintf("gvr_lollipop: '%s' - %d rows in gene", gene, nrow(dt)))

  dt <- dt[Variant_Classification %in% vc_use]
  if (nrow(dt) == 0L) {
    warning(sprintf("gvr_lollipop: no protein-altering variants for '%s' in input MAF", gene))
    return(.empty_plot(gene))
  }
  if (isTRUE(verbose)) message(sprintf("gvr_lollipop: '%s' - %d rows after VC filter", gene, nrow(dt)))

  # ---- Parse amino-acid position from HGVSp_Short (^p\.[A-Z*]([0-9]+)) ----
  aa_pos <- suppressWarnings(as.integer(sub("^p\\.[A-Z*]([0-9]+).*", "\\1", dt$HGVSp_Short)))
  # rows with no parseable HGVSp_Short position
  bad_idx <- is.na(aa_pos) | !grepl("^p\\.[A-Z*][0-9]+", dt$HGVSp_Short)
  if (any(bad_idx)) {
    if (isTRUE(verbose)) message(sprintf("gvr_lollipop: '%s' - dropping %d rows with unparseable HGVSp_Short (e.g. %s)",
                                          gene, sum(bad_idx),
                                          paste(utils::head(unique(dt$HGVSp_Short[bad_idx]), 3), collapse = ", ")))
    dt     <- dt[!bad_idx]
    aa_pos <- aa_pos[!bad_idx]
  }
  if (nrow(dt) == 0L) {
    warning(sprintf("gvr_lollipop: no parseable HGVSp_Short positions for '%s'", gene))
    return(.empty_plot(gene))
  }

  dt$.__aa_pos__ <- aa_pos

  # ---- Phase L: capture whether the user supplied a length (precedence input) ----
  .user_supplied_length <- !is.null(protein_length)

  # ---- Derive protein length from Protein_position (regex hardened: requires '/') ----
  if (is.null(protein_length)) {
    pp <- as.character(dt$Protein_position)
    with_slash <- grepl("/", pp, fixed = TRUE)
    totals <- suppressWarnings(as.integer(sub("^[^/]*/", "", pp[with_slash])))
    totals <- totals[!is.na(totals) & totals > 0]
    if (length(totals) > 0L) {
      tab <- sort(table(totals), decreasing = TRUE)
      protein_length <- as.integer(names(tab)[1])
    } else {
      protein_length <- as.integer(
        ceiling(max(dt$.__aa_pos__, na.rm = TRUE) * 1.1))
    }
  }
  protein_length <- max(protein_length, max(dt$.__aa_pos__, na.rm = TRUE))

  # ---- Build per-dot table: one dot per surviving sample-variant ----
  # stack y = within-position rank
  data.table::setDT(dt)
  dt[, .__hgvsp__ := HGVSp_Short]
  dt[, .__vc__    := Variant_Classification]
  dt[, .__sample__ := as.character(Tumor_Sample_Barcode)]
  # within each aa position, assign sequential y from 1
  dt <- dt[order(.__aa_pos__, .__vc__, .__sample__)]
  dt[, .__y__ := seq_len(.N), by = .__aa_pos__]

  # stem height = stack height per position
  pos_height <- dt[, list(.__top__ = max(.__y__)), by = .__aa_pos__]

  # ---- Top label table ----
  label_n <- if (is.infinite(label_top)) Inf else as.integer(label_top)
  if (label_n > 0L) {
    pos_counts <- dt[, list(.n = .N), by = .__aa_pos__]
    pos_counts <- pos_counts[order(-.n)]
    keep_n <- if (is.infinite(label_n)) nrow(pos_counts) else min(label_n, nrow(pos_counts))
    top_positions <- pos_counts$.__aa_pos__[seq_len(keep_n)]

    # For each top position, pick the most-common HGVSp_Short + count any other distinct ones
    label_dt <- dt[.__aa_pos__ %in% top_positions, list(.__hgvsp__, .__aa_pos__)]
    label_rows <- vector("list", length(top_positions))
    for (i in seq_along(top_positions)) {
      p <- top_positions[i]
      sub_h <- table(label_dt$.__hgvsp__[label_dt$.__aa_pos__ == p])
      sub_h <- sort(sub_h, decreasing = TRUE)
      top_h <- names(sub_h)[1]
      n_more <- length(sub_h) - 1L
      lbl <- if (n_more > 0L) sprintf("%s (+%d more)", top_h, n_more) else top_h
      label_rows[[i]] <- data.frame(.__aa_pos__ = p,
                                    .__top__   = pos_height$.__top__[pos_height$.__aa_pos__ == p],
                                    .__label__ = lbl,
                                    stringsAsFactors = FALSE)
    }
    label_df <- do.call(rbind, label_rows)
  } else {
    label_df <- NULL
  }

  # ---- Class palette: resolved via .resolve_variant_palette ----
  classes_present <- sort(unique(dt$.__vc__))
  col_map <- .resolve_variant_palette(variant_palette, classes_present)

  # ---- x-axis breaks: aim for ~10-12 ticks regardless of protein length ----
  .nice_step <- function(L) {
    raw <- L / 10
    pow <- 10 ^ floor(log10(raw))
    base <- raw / pow
    nice <- if (base < 1.5) 1 else if (base < 3) 2 else if (base < 7) 5 else 10
    nice * pow
  }
  brk_step <- .nice_step(protein_length)
  x_breaks <- seq(0, ceiling(protein_length / brk_step) * brk_step, by = brk_step)

  # ---- Protein-body bar geometry (overlapping y=0) ----
  y_max_stack <- max(pos_height$.__top__)
  # bar half-thickness = ~4% of stack height, clamped to [0.15, 0.5]
  bar_half <- max(0.22, min(0.5, y_max_stack * 0.05))
  # Top y-bound (does not depend on label fit; we compute y_lower later).
  # Resolve user-requested mode + position here; .dlm_resolved comes later.
  .dlm_requested <- match.arg(domain_label_mode)
  domain_label_position <- match.arg(domain_label_position)
  # Give a tad more head-room so the tallest stack doesn't touch the top edge.
  y_upper  <- y_max_stack + max(1, y_max_stack * 0.10)

  bar_df <- data.frame(xmin = 0, xmax = protein_length,
                       ymin = -bar_half, ymax = bar_half,
                       stringsAsFactors = FALSE)

  # ---- Optional domain rectangles ----
  domains_df <- NULL
  if (!is.null(domains)) {
    if (!is.data.frame(domains)) {
      stop("gvr_lollipop: 'domains' must be a data.frame / data.table.")
    }
    dom <- as.data.frame(domains, stringsAsFactors = FALSE)
    miss_dom <- setdiff(c("start", "end"), names(dom))
    if (length(miss_dom) > 0L) {
      stop(sprintf("gvr_lollipop: 'domains' missing required column(s): %s",
                   paste(miss_dom, collapse = ", ")))
    }
    dom$start <- suppressWarnings(as.numeric(dom$start))
    dom$end   <- suppressWarnings(as.numeric(dom$end))
    # drop invalid intervals
    bad_iv <- is.na(dom$start) | is.na(dom$end) | dom$end < dom$start
    if (any(bad_iv)) {
      if (isTRUE(verbose))
        message(sprintf("gvr_lollipop: domains - dropping %d row(s) with NA / inverted start..end",
                        sum(bad_iv)))
      dom <- dom[!bad_iv, , drop = FALSE]
    }
    # drop intervals entirely outside protein
    out_iv <- dom$end < 0 | dom$start > protein_length
    if (any(out_iv)) {
      if (isTRUE(verbose))
        message(sprintf("gvr_lollipop: domains - dropping %d row(s) outside [0, %d]",
                        sum(out_iv), protein_length))
      dom <- dom[!out_iv, , drop = FALSE]
    }
    if (nrow(dom) > 0L) {
      # clip to protein bounds
      dom$start <- pmax(dom$start, 0)
      dom$end   <- pmin(dom$end,   protein_length)

      # name column (optional) -- character or empty
      if (!"name" %in% names(dom)) {
        dom$name <- rep("", nrow(dom))
      } else {
        dom$name <- as.character(dom$name)
        dom$name[is.na(dom$name)] <- ""
      }

      # color column (optional) -- if missing/empty, fill from domain_palette
      .auto_dom_colors <- .resolve_domain_palette(domain_palette, dom)
      if (!"color" %in% names(dom)) {
        dom$color <- .auto_dom_colors
      } else {
        dom$color <- as.character(dom$color)
        empty_col <- is.na(dom$color) | dom$color == ""
        if (any(empty_col)) {
          dom$color[empty_col] <- .auto_dom_colors[empty_col]
        }
      }

      # build geom_rect-friendly frame; domains slightly taller than bar
      dom_half <- bar_half * 2.4  # clearly taller than bar (maftools-style highlight)
      domains_df <- data.frame(xmin = dom$start, xmax = dom$end,
                               ymin = -dom_half, ymax = dom_half,
                               fill = dom$color,
                               name = dom$name,
                               stringsAsFactors = FALSE)
      # Label visible if rectangle width is at least the smaller of:
      #   (a) domain_label_min_frac * protein_length  (default 5%), and
      #   (b) a 50 aa absolute floor (so long proteins still label small domains).
      # The intent: percentage rule keeps short proteins from labelling tiny dots,
      # while absolute floor keeps long proteins from suppressing every domain.
      .min_frac <- if (is.null(domain_label_min_frac) || !is.finite(domain_label_min_frac))
                     0.05 else max(0, min(1, as.numeric(domain_label_min_frac)))
      .min_aa   <- min(protein_length * .min_frac, 50)
      domains_df$show_label <- (domains_df$xmax - domains_df$xmin) >=
                               .min_aa & nzchar(domains_df$name)
      domains_df$xmid <- (domains_df$xmin + domains_df$xmax) / 2
      # carry accession through (NA when user passed a data.frame without it)
      domains_df$accession <- if ("accession" %in% names(dom)) as.character(dom$accession) else NA_character_

      if (isTRUE(verbose))
        message(sprintf("gvr_lollipop: '%s' - drawing %d domain(s)", gene, nrow(domains_df)))
        # Hint: if any domains exist but all are below the name-label threshold,
        # tell the user how to make labels appear.
        if (nrow(domains_df) > 0L && !any(domains_df$show_label))
          message(sprintf("gvr_lollipop: '%s' - all %d domain(s) are below the label-width threshold (%.0f aa). Try domain_label_mode = \"number\" or lower domain_label_min_frac.",
                          gene, nrow(domains_df), .min_aa))
    } else {
      if (isTRUE(verbose))
        message(sprintf("gvr_lollipop: '%s' - no valid domain rows after filtering",
                        gene))
    }
  }

  # ---- Resolve domain_label_mode (handles "auto" heuristic) ----
  # The first dispatch site (above) only knew the user-requested mode.
  # Now that domains_df is finalized we know n_dom and can resolve "auto".
  .n_dom <- if (is.null(domains_df)) 0L else nrow(domains_df)
  .dlm_resolved <- if (.dlm_requested == "auto") {
    if (protein_length > 2000L && .n_dom >= 5L) {
      if (isTRUE(verbose))
        message(sprintf(
          "gvr_lollipop: '%s' - auto-selecting domain_label_mode = \"number\" (protein length %d aa, %d domain(s))",
          gene, as.integer(protein_length), .n_dom))
      "number"
    } else {
      if (isTRUE(verbose))
        message(sprintf(
          "gvr_lollipop: '%s' - auto-selecting domain_label_mode = \"name\" (protein length %d aa, %d domain(s))",
          gene, as.integer(protein_length), .n_dom))
      "name"
    }
  } else {
    .dlm_requested
  }

  # ---- Precompute label-fit decisions and choose tight vs generous y_lower ----
  # Phase K: y_lower must reflect actual overflow, not just requested mode.
  # We mimic the renderer's per-domain fit heuristic here, then reuse the
  # same .fits_preview decisions when the renderer draws labels.
  .text_size_preview <- if (.dlm_resolved == "number") base_size * 0.28 else base_size * 0.23

  if (.dlm_resolved == "none" || is.null(domains_df) || nrow(domains_df) == 0L) {
    .n_overflow_preview <- 0L
  } else if (.dlm_resolved == "number") {
    # Numbers (1, 2, ...) always fit; mark every row as fit
    domains_df$.fits_preview <- rep(TRUE, nrow(domains_df))
    .n_overflow_preview <- 0L
  } else {
    # "name" or "id" -- estimate fit for labeled rows only
    domains_df$.fits_preview <- rep(NA, nrow(domains_df))
    .lab_idx <- which(domains_df$show_label)
    if (length(.lab_idx) == 0L) {
      .n_overflow_preview <- 0L
    } else {
      if (.dlm_resolved == "id") {
        .lbl_preview <- ifelse(
          is.na(domains_df$accession[.lab_idx]) | !nzchar(domains_df$accession[.lab_idx]),
          domains_df$name[.lab_idx],
          domains_df$accession[.lab_idx])
      } else {  # "name"
        .lbl_preview <- .gvr_abbrev_domain(domains_df$name[.lab_idx],
                                          abbrev = isTRUE(domain_name_abbrev))
      }
      if (domain_label_position == "inside") {
        .char_aa  <- (protein_length / 80) * (.text_size_preview / 11)
        .est_w    <- nchar(.lbl_preview) * .char_aa
        .rect_w   <- domains_df$xmax[.lab_idx] - domains_df$xmin[.lab_idx]
        .fits_vec <- .est_w <= .rect_w * 0.9
        domains_df$.fits_preview[.lab_idx] <- .fits_vec
        .n_overflow_preview <- sum(!.fits_vec)
      } else {  # "below"
        # In below mode every labeled row renders below by construction
        domains_df$.fits_preview[.lab_idx] <- FALSE
        .n_overflow_preview <- length(.lab_idx)
      }
    }
  }

  # Decide y_lower based on whether any label needs the below-band
  .dom_floor <- bar_half * 2.4
  .need_below_space <-
    (domain_label_position == "below"  && .dlm_resolved %in% c("name", "id")) ||
    (domain_label_position == "inside" && .n_overflow_preview > 0L)
  y_lower <- if (.need_below_space)
               -.dom_floor - max(0.50, y_max_stack * 0.35)   # generous (legacy)
             else
               -.dom_floor - max(0.05, y_max_stack * 0.02)   # tight (Phase K)


  # ---- Hotspot detection (clusters of >= hotspot_min_n positions within hotspot_window aa) ----
  # Validate inputs: any non-finite / <=0 value disables hotspots safely.
  .hw <- suppressWarnings(as.numeric(hotspot_window))
  .hm <- suppressWarnings(as.numeric(hotspot_min_n))
  .hotspots_enabled <- isTRUE(is.finite(.hw) && .hw > 0) &&
                       isTRUE(!is.na(.hm) && .hm > 0)  # Inf passes >0 test
  if (isTRUE(verbose)) {
    if (!is.null(hotspot_window) && (!is.finite(.hw) || .hw <= 0))
      message(sprintf("gvr_lollipop: invalid hotspot_window (%s); hotspot detection disabled.",
                      as.character(hotspot_window)))
    if (!is.null(hotspot_min_n) && (is.na(.hm) || .hm <= 0))
      message(sprintf("gvr_lollipop: invalid hotspot_min_n (%s); hotspot detection disabled.",
                      as.character(hotspot_min_n)))
  }
  # If hotspot_min_n == Inf, detection is officially disabled (no message).
  hotspot_df <- data.frame(xmin = numeric(0), xmax = numeric(0),
                           xmid = numeric(0), n_in_window = integer(0),
                           stringsAsFactors = FALSE)
  if (.hotspots_enabled && is.finite(.hm)) {
    # `pos_counts` holds unique aa positions (one row per .__aa_pos__).
    .unique_pos <- sort(unique(pos_counts$.__aa_pos__))
    if (length(.unique_pos) >= .hm) {
      .half_w <- .hw / 2
      # For each position, count distinct positions inside its centered window.
      .in_win <- vapply(.unique_pos, function(p) {
        sum(.unique_pos >= (p - .half_w) & .unique_pos <= (p + .half_w))
      }, integer(1))
      .seed_pos <- .unique_pos[.in_win >= .hm]
      if (length(.seed_pos) > 0L) {
        # Build candidate windows (one per seed), then greedy-merge overlaps.
        .cand <- data.frame(
          xmin        = pmax(0, .seed_pos - .half_w),
          xmax        = pmin(protein_length, .seed_pos + .half_w),
          n_in_window = .in_win[.in_win >= .hm],
          stringsAsFactors = FALSE
        )
        .cand <- .cand[order(.cand$xmin), , drop = FALSE]
        # Greedy merge: walk sorted candidates, extend the current band while
        # the next xmin <= current xmax (i.e. overlap or touching).
        .merged <- vector("list", nrow(.cand))
        .cur <- .cand[1, , drop = FALSE]
        .k <- 1L
        if (nrow(.cand) > 1L) {
          for (.i in 2:nrow(.cand)) {
            if (.cand$xmin[.i] <= .cur$xmax) {
              .cur$xmax        <- max(.cur$xmax, .cand$xmax[.i])
              .cur$n_in_window <- max(.cur$n_in_window, .cand$n_in_window[.i])
            } else {
              .merged[[.k]] <- .cur; .k <- .k + 1L
              .cur <- .cand[.i, , drop = FALSE]
            }
          }
        }
        .merged[[.k]] <- .cur
        hotspot_df <- do.call(rbind, .merged[seq_len(.k)])
        hotspot_df$xmid <- (hotspot_df$xmin + hotspot_df$xmax) / 2
        # Reorder columns to canonical layout
        hotspot_df <- hotspot_df[, c("xmin", "xmax", "xmid", "n_in_window"),
                                 drop = FALSE]
        rownames(hotspot_df) <- NULL
        if (isTRUE(verbose))
          message(sprintf("gvr_lollipop: '%s' - %d hotspot band(s) detected (window %g aa, min %g positions)",
                          gene, nrow(hotspot_df), .hw, .hm))
      }
    }
  }

  # ---- Build ggplot ----
  # convert internal dot-prefixed columns to plain names ggplot can see cleanly
  dot_df <- data.frame(aa_pos = dt$.__aa_pos__,
                       y      = dt$.__y__,
                       vc     = dt$.__vc__,
                       hgvsp  = dt$.__hgvsp__,
                       stringsAsFactors = FALSE)
  stem_df <- data.frame(aa_pos = pos_height$.__aa_pos__,
                        top    = pos_height$.__top__,
                        stringsAsFactors = FALSE)

  p <- ggplot2::ggplot() +
    # protein-body bar (drawn first so stems/dots/domains sit on top)
    ggplot2::geom_rect(data = bar_df,
                       ggplot2::aes(xmin = xmin, xmax = xmax,
                                    ymin = ymin, ymax = ymax),
                       fill = bar_color,
                       color = if (is.na(bar_border)) NA else bar_border,
                       linewidth = 0.3)

  # mutation-hotspot bands (drawn on top of the bar, behind domains/stems/dots).
  # Soft full-height vertical wash highlights clustered-mutation regions
  # without competing with the variant-class colour scale.
  if (nrow(hotspot_df) > 0L) {
    p <- p + ggplot2::geom_rect(
      data        = hotspot_df,
      ggplot2::aes(xmin = xmin, xmax = xmax,
                   ymin = -bar_half * 2.4,
                   ymax =  y_upper),
      fill        = "#FD9BED",  # phylo magenta
      color       = NA,
      alpha       = 0.15,
      inherit.aes = FALSE)
  }

  # optional domain rectangles (on top of bar, but below stems)
  if (!is.null(domains_df) && nrow(domains_df) > 0L) {
    # Use a single thin dark outline so adjacent same-coloured domains stay
    # visually separable. (Per-domain outlines would collide with the variant
    # color scale used by geom_point.)
    p <- p + ggplot2::geom_rect(data = domains_df,
                                ggplot2::aes(xmin = xmin, xmax = xmax,
                                             ymin = ymin, ymax = ymax,
                                             fill = fill),
                                color = "grey25",
                                linewidth = 0.4) +
             ggplot2::scale_fill_identity()
    # domain labels: behaviour controlled by domain_label_mode + domain_label_position
    .dlm <- .dlm_resolved
    if (.dlm != "none") {
      # Build lab_dom + the label string in .lbl per mode
      if (.dlm == "number") {
        lab_dom <- domains_df
        lab_dom$.lbl <- as.character(seq_len(nrow(lab_dom)))
      } else if (.dlm == "id") {
        lab_dom <- domains_df[domains_df$show_label, , drop = FALSE]
        lab_dom$.lbl <- ifelse(is.na(lab_dom$accession) | !nzchar(lab_dom$accession),
                               lab_dom$name, lab_dom$accession)
      } else {  # "name"
        lab_dom <- domains_df[domains_df$show_label, , drop = FALSE]
        lab_dom$.lbl <- .gvr_abbrev_domain(lab_dom$name, abbrev = isTRUE(domain_name_abbrev))
      }

      if (nrow(lab_dom) > 0L) {
        # ---- Decide inside vs below per domain ----
        # For position="below", everything goes below (legacy path).
        # For position="inside", estimate fit per domain; overflows -> below.
        if (.dlm == "number") {
          # numbers are short; default to inside unless explicitly below
          .text_size <- base_size * 0.28
        } else {
          .text_size <- base_size * 0.23
        }

        # Phase K: reuse .fits_preview attached during the y_lower precompute.
        # Identical inputs guarantee identical decisions; this also lets us
        # avoid the same heuristic getting computed twice.
        if (".fits_preview" %in% names(domains_df)) {
          if (.dlm == "number") {
            .fits <- as.logical(domains_df$.fits_preview)
          } else {
            .fits <- as.logical(domains_df$.fits_preview[domains_df$show_label])
          }
          # Defensive: if NAs ended up here (e.g. unlabeled rows leaked through),
          # treat NA as not-fit so they render below rather than failing silently.
          .fits[is.na(.fits)] <- FALSE
        } else {
          # Safety fallback: recompute identically to the precompute heuristic.
          if (domain_label_position == "inside" && nrow(lab_dom) > 0L) {
            .char_aa <- (protein_length / 80) * (.text_size / 11)
            .est_w   <- nchar(lab_dom$.lbl) * .char_aa
            .rect_w  <- lab_dom$xmax - lab_dom$xmin
            .fits    <- .est_w <= .rect_w * 0.9
          } else {
            .fits <- rep(FALSE, nrow(lab_dom))
          }
        }

        inside_dom <- lab_dom[ .fits, , drop = FALSE]
        below_dom  <- lab_dom[!.fits, , drop = FALSE]

        n_overflow <- nrow(below_dom)
        if (domain_label_position == "inside" && n_overflow > 0L && isTRUE(verbose)) {
          message(sprintf(
            "gvr_lollipop: '%s' - %d domain label(s) overflow rectangle; rendering below",
            gene, n_overflow))
        }

        # ---- Inside renderer: centered geom_text with auto-contrast color ----
        if (nrow(inside_dom) > 0L) {
          inside_dom$.text_col <- .pick_contrast_color(inside_dom$fill)
          p <- p + ggplot2::geom_text(
                     data = inside_dom,
                     ggplot2::aes(x = xmid, y = 0, label = .lbl),
                     color    = inside_dom$.text_col,
                     size     = .text_size,
                     fontface = if (.dlm == "number") "bold" else "plain",
                     vjust    = 0.5, hjust = 0.5,
                     inherit.aes = FALSE)
        }

        # ---- Below renderer: ggrepel with leader lines (legacy path) ----
        if (nrow(below_dom) > 0L) {
          if (.dlm == "number") {
            if (requireNamespace("ggrepel", quietly = TRUE) && nrow(below_dom) >= 4L) {
              p <- p + ggrepel::geom_text_repel(
                         data = below_dom,
                         ggplot2::aes(x = xmid, y = 0, label = .lbl),
                         size               = base_size * 0.28,
                         color              = "black",
                         fontface           = "bold",
                         direction          = "y",
                         nudge_y            = bar_half * 0.8,
                         ylim               = c(bar_half * 0.6, Inf),
                         segment.size       = 0.2,
                         segment.color      = "grey60",
                         box.padding        = 0.10,
                         point.padding      = 0.02,
                         min.segment.length = 0.05,
                         max.overlaps       = Inf,
                         seed               = 42L)
            } else {
              p <- p + ggplot2::geom_text(data = below_dom,
                                          ggplot2::aes(x = xmid, y = 0, label = .lbl),
                                          size  = base_size * 0.30,
                                          color = "black",
                                          fontface = "bold",
                                          vjust = 0.5, hjust = 0.5)
            }
          } else if (requireNamespace("ggrepel", quietly = TRUE)) {
            .y_start <- -bar_half - max(0.05, y_max_stack * 0.04)
            p <- p + ggrepel::geom_text_repel(
                       data = below_dom,
                       ggplot2::aes(x = xmid, y = .y_start, label = .lbl),
                       size               = base_size * 0.23,
                       color              = "black",
                       direction          = "y",
                       nudge_y            = -max(0.10, y_max_stack * 0.06),
                       segment.size       = 0.25,
                       segment.color      = "grey50",
                       box.padding        = 0.30,
                       point.padding      = 0.05,
                       min.segment.length = 0,
                       max.overlaps       = Inf,
                       seed               = 42L)
          } else {
            p <- p + ggplot2::geom_text(data = below_dom,
                                        ggplot2::aes(x = xmid, y = 0, label = .lbl),
                                        size  = base_size * 0.23,
                                        color = "black",
                                        vjust = 0.5, hjust = 0.5)
          }
        }
      }
      # In 'number' mode, also draw a side caption mapping number -> name
      # (rendered as a labelled fill scale legend via a sentinel aesthetic).
      if (.dlm == "number") {
        # rebuild fill so legend can show numbered legend entries
        legend_lvl <- sprintf("%d. %s", seq_len(nrow(domains_df)),
                              ifelse(nzchar(domains_df$name), domains_df$name,
                                     ifelse(is.na(domains_df$accession), "(unnamed)",
                                            domains_df$accession)))
        # keep order via factor levels
        domains_df$.legend <- factor(legend_lvl, levels = legend_lvl)
        # overlay an invisible rect just to inject the legend (the real fill
        # is the geom_rect above using scale_fill_identity()). Use a second
        # ggplot layer with discrete fill on top.
        if (requireNamespace("ggnewscale", quietly = TRUE)) {
          # Add a second fill scale that maps each domain to its colour with
          # a numbered legend. The in-plot rect is invisible (alpha = 0) but
          # the legend key glyph is forced to full opacity via override.aes
          # so the swatches actually show the domain colours.
          p <- p + ggnewscale::new_scale_fill() +
               ggplot2::geom_rect(data = domains_df,
                                  ggplot2::aes(xmin = xmin, xmax = xmax,
                                               ymin = ymin, ymax = ymax,
                                               fill = .legend),
                                  color = NA, alpha = 0) +
               ggplot2::scale_fill_manual(
                 name   = "Domains",
                 values = setNames(domains_df$fill,
                                   as.character(domains_df$.legend))) +
               ggplot2::guides(fill = ggplot2::guide_legend(
                 override.aes = list(alpha = 1, color = "grey25")))
        } else if (isTRUE(verbose)) {
          message("gvr_lollipop: install 'ggnewscale' to get a numbered domain legend; rendering numbers only.")
        }
      }
    }
  }

  p <- p +
    # stems (from top of bar up to stack top)
    ggplot2::geom_segment(data = stem_df,
                          ggplot2::aes(x = aa_pos, xend = aa_pos,
                                       y = bar_half, yend = top),
                          color = stem_color,
                          alpha = as.numeric(stem_alpha),
                          linewidth = 0.5) +
    # dots
    ggplot2::geom_point(data = dot_df,
                        ggplot2::aes(x = aa_pos, y = y, color = vc),
                        size = point_size) +
    ggplot2::scale_color_manual(values = col_map, name = "Variant class") +
    ggplot2::scale_x_continuous(limits = c(0, protein_length),
                                breaks = x_breaks, expand = c(0.01, 0.01)) +
    ggplot2::scale_y_continuous(
      limits = c(y_lower, y_upper),
      breaks = function(lim) {
        ub <- max(1, ceiling(lim[2]))
        seq(0, ub, by = max(1, ceiling(ub / 5)))
      },
      expand = c(0, 0)) +
    ggplot2::labs(title = gene,
                  subtitle = sprintf("Protein length: %d aa  -  %d variant%s in %d sample%s",
                                     protein_length,
                                     nrow(dt),
                                     if (nrow(dt) == 1L) "" else "s",
                                     length(unique(dt$.__sample__)),
                                     if (length(unique(dt$.__sample__)) == 1L) "" else "s"),
                  x = "Amino-acid position",
                  y = "Number of sample-variants",
                  color = "Variant class") +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(face = "bold",
                                                 size = base_size * 1.25,
                                                 hjust = 0,
                                                 margin = ggplot2::margin(b = 2)),
      plot.subtitle      = ggplot2::element_text(size = base_size * 0.85,
                                                 color = "grey35",
                                                 hjust = 0,
                                                 margin = ggplot2::margin(b = 8)),
      axis.title         = ggplot2::element_text(size = axis_title_size, color = axis_title_color),
      axis.text          = ggplot2::element_text(size = axis_text_size, color = axis_text_color),
      axis.line          = ggplot2::element_line(color = axis_line_color, linewidth = axis_line_width),
      axis.ticks         = ggplot2::element_line(color = axis_line_color, linewidth = axis_line_width),
      legend.title       = ggplot2::element_text(face = "bold", size = base_size * 0.85),
      legend.text        = ggplot2::element_text(size = base_size * 0.80),
      legend.key.size    = ggplot2::unit(0.7, "lines"),
      legend.position    = "right",
      legend.justification = c(0, 0.5),
      legend.box.spacing = ggplot2::unit(0.3, "lines"),
      plot.margin        = ggplot2::margin(t = 8, r = 12, b = 6, l = 6)
    )

  if (!is.null(label_df) && nrow(label_df) > 0L) {
    label_plot_df <- data.frame(aa_pos = label_df$.__aa_pos__,
                                top    = label_df$.__top__,
                                label  = label_df$.__label__,
                                stringsAsFactors = FALSE)
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(data = label_plot_df,
                                        ggplot2::aes(x = aa_pos, y = top,
                                                     label = label),
                                        size = base_size * 0.3,
                                        nudge_y = 0.5,
                                        segment.size = 0.2,
                                        max.overlaps = Inf)
    } else {
      p <- p + ggplot2::geom_text(data = label_plot_df,
                                  ggplot2::aes(x = aa_pos, y = top, label = label),
                                  size = base_size * 0.3,
                                  vjust = -0.6, hjust = 0.5)
    }
  }

  # ---- File output (default: auto-save, matching gvr_plot behaviour) ----
  if (!is.null(out_dir) && !is.null(out_prefix)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    svg_path <- file.path(out_dir, sprintf("%s_lollipop.svg", out_prefix))
    png_path <- file.path(out_dir, sprintf("%s_lollipop.png", out_prefix))

    # announce overwrite (matching gvr_plot convention)
    if (file.exists(png_path) && isTRUE(verbose))
      message(sprintf("gvr_lollipop: overwriting existing %s", png_path))

    svg_written <- .fuse_save(svg_path, function(tmp) {
      grDevices::svg(tmp, width = width, height = height)
      print(p)
      grDevices::dev.off()
    })
    png_written <- .fuse_save(png_path, function(tmp) {
      ggplot2::ggsave(tmp, plot = p, width = width, height = height,
                      units = "in", dpi = dpi, device = "png")
    })

    if (isTRUE(verbose)) {
      if (!is.na(svg_written)) message(sprintf("gvr_lollipop: wrote %s", svg_written))
      if (!is.na(png_written)) message(sprintf("gvr_lollipop: wrote %s", png_written))
    }
  }

  p
}
