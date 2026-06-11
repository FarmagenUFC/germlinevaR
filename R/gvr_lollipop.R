#' Per-gene amino-acid lollipop plot for a germline MAF
#'
#' @description
#' Builds a ggplot2 lollipop plot of every protein-altering variant of a single
#' gene in a `read.gvr()` / `gvr_filter()` / `gvr_novel()` MAF. Each amino-acid
#' position is a stem; each sample-variant carried at that position is one dot
#' stacked on the stem; dots are coloured by `Variant_Classification` using the
#' package's colourblind-safe palette. Optionally writes both an `.svg` and a
#' `.png` to disk.
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
#' equals the count of sample-variants at that position. Stems run from `y = 0`
#' up to the top of the stack.
#'
#' Labels: the top-`label_top` positions by count are annotated with their
#' `HGVSp_Short` (e.g. `p.R175H`). When multiple distinct `HGVSp_Short` share
#' a position, the most-common one is shown followed by `(+N more)`. Set
#' `label_top = 0` to disable labels, `label_top = Inf` to label every
#' position. If \pkg{ggrepel} is installed, labels are placed via
#' `ggrepel::geom_text_repel()`; otherwise a plain `geom_text()` is used.
#'
#' Empty-gene behaviour: if no row matches the gene, or none survives the
#' `vc_keep` filter, or none has a parseable amino-acid position, the function
#' issues a warning and returns a ggplot2 object with a single centered
#' "No protein-altering variants for <GENE>" annotation. The user can still
#' `ggsave()` it normally.
#'
#' Optional file output: when BOTH `out_dir` and `out_prefix` are non-NULL,
#' the function writes:
#' \itemize{
#'   \item `<out_dir>/<out_prefix>_lollipop.svg` - vector (lossless)
#'   \item `<out_dir>/<out_prefix>_lollipop.png` - raster at `dpi`
#' }
#' Both files are first rendered to `tempdir()` then `cp`'d to `out_dir`, so
#' S3-backed FUSE mounts (e.g. `/mnt/results/`) work without 0-byte issues.
#' The function always returns the ggplot2 object regardless of whether files
#' were written.
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
#' @param label_top Integer(1). Number of top-counted positions to label.
#'   `0` disables labels, `Inf` labels every position. Default `5L`.
#' @param point_size Numeric(1). Dot size. Default `3`.
#' @param stem_color Character(1). Stem (vertical line) colour. Default
#'   `"grey50"`.
#' @param base_size Numeric(1). ggplot2 base font size. Default `12`.
#' @param out_dir Character(1) or `NULL`. If both `out_dir` and `out_prefix`
#'   are non-NULL, the plot is written to disk (see Details). Default `NULL`
#'   (no file output).
#' @param out_prefix Character(1) or `NULL`. Filename prefix; output files are
#'   `<out_prefix>_lollipop.svg` / `<out_prefix>_lollipop.png`. Default `NULL`.
#' @param width Numeric(1). Plot width in inches. Default `10`.
#' @param height Numeric(1). Plot height in inches. Default `4`.
#' @param dpi Numeric(1). PNG resolution. Default `300`.
#' @param verbose Logical(1). If `TRUE` (default), emit progress messages
#'   (counts, dropped rows, file paths).
#'
#' @return A ggplot2 object. If `out_dir` + `out_prefix` are passed, the SVG
#'   and PNG files are also written to disk as a side effect.
#'
#' @seealso [read.gvr()], [gvr_filter()], [gvr_novel()], [gvr_summary()],
#'   [gvr_plot()]
#'
#' @examples
#' \dontrun{
#'   maf <- read.gvr("vcf_dir/", pattern = "\\.vep\\.vcf\\.gz$")
#'   f   <- gvr_filter(maf)
#'
#'   ## inline plot
#'   p <- gvr_lollipop(f, "MUC16")
#'   print(p)
#'
#'   ## save SVG + PNG
#'   gvr_lollipop(f, "MUC16",
#'                out_dir    = "results/lollipops",
#'                out_prefix = "MUC16")
#'
#'   ## only missense
#'   gvr_lollipop(f, "TP53", vc_keep = "Missense_Mutation")
#' }
#'
#' @importFrom data.table as.data.table copy is.data.table setDT :=
#' @importFrom ggplot2 ggplot aes geom_segment geom_point geom_text annotate
#'   scale_color_manual scale_x_continuous scale_y_continuous labs theme_classic
#'   theme element_text element_blank ggsave
#' @importFrom grDevices svg dev.off
#' @export
gvr_lollipop <- function(maf, gene,
                         vc_keep        = NULL,
                         protein_length = NULL,
                         label_top      = 5L,
                         point_size     = 3,
                         stem_color     = "grey50",
                         base_size      = 12,
                         out_dir        = NULL,
                         out_prefix     = NULL,
                         width          = 10,
                         height         = 4,
                         dpi            = 300,
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
  if (length(miss_req) > 0) {
    stop(sprintf("gvr_lollipop: required column(s) not found: %s",
                 paste(miss_req, collapse = ", ")))
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

  # ---- Derive protein length from Protein_position ----
  if (is.null(protein_length)) {
    totals <- suppressWarnings(as.integer(sub("^[^/]*/", "", dt$Protein_position)))
    totals <- totals[!is.na(totals) & totals > 0]
    if (length(totals) > 0L) {
      tab <- sort(table(totals), decreasing = TRUE)
      protein_length <- as.integer(names(tab)[1])
    } else {
      protein_length <- as.integer(ceiling(max(dt$.__aa_pos__, na.rm = TRUE) * 1.1))
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

  # ---- Class palette: present classes mapped via GVR_CLASS_COLORS ----
  classes_present <- sort(unique(dt$.__vc__))
  col_map <- GVR_CLASS_COLORS[classes_present]
  col_map[is.na(col_map)] <- GVR_CLASS_COLORS[["Other"]]
  names(col_map) <- classes_present

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
    ggplot2::geom_segment(data = stem_df,
                          ggplot2::aes(x = aa_pos, xend = aa_pos,
                                       y = 0, yend = top),
                          color = stem_color, linewidth = 0.5) +
    ggplot2::geom_point(data = dot_df,
                        ggplot2::aes(x = aa_pos, y = y, color = vc),
                        size = point_size) +
    ggplot2::scale_color_manual(values = col_map, name = "Variant_Classification") +
    ggplot2::scale_x_continuous(limits = c(0, protein_length),
                                breaks = x_breaks, expand = c(0.01, 0.01)) +
    ggplot2::scale_y_continuous(limits = c(0, max(pos_height$.__top__) + 1),
                                breaks = function(lim) seq(0, ceiling(lim[2]), by = max(1, ceiling(lim[2] / 5))),
                                expand = c(0, 0)) +
    ggplot2::labs(title = sprintf("%s (protein length: %d aa)", gene, protein_length),
                  x = "Amino-acid position",
                  y = "Number of sample-variants") +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0),
                   legend.position = "right")

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

  # ---- Optional file output ----
  if (!is.null(out_dir) && !is.null(out_prefix)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    svg_path <- file.path(out_dir, sprintf("%s_lollipop.svg", out_prefix))
    png_path <- file.path(out_dir, sprintf("%s_lollipop.png", out_prefix))

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
