#' Cohort oncoplot from a germline MAF (read.gvr / gvr_filter output)
#'
#' @description
#' Draws a maftools-style oncoplot from a MAF-style table - the output of
#' [read.gvr()], or of [gvr_filter()] - and writes it to a PNG file. Rows are the
#' top-`top_n` genes (ranked by number of distinct samples mutated, then by variant
#' count); columns are samples. Each gene x sample cell shows the single MOST-SEVERE
#' `Variant_Classification` observed for that gene in that sample. The plot is drawn
#' with \pkg{ComplexHeatmap}.
#'
#' @details
#' This is the standalone oncoplot previously produced inside [gvr_summary()]. It
#' needs only the `Hugo_Symbol` and `Variant_Classification` columns (plus the
#' per-sample column).
#'
#' Cell collapse: when a gene has several variant classes in one sample, the cell is
#' coloured by the most severe class, using this ranking (high to low):
#' Translation_Start_Site, Nonsense_Mutation, Nonstop_Mutation, Splice_Site,
#' Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins, Missense_Mutation,
#' Splice_Region, Protein_altering_variant, Silent, 5'UTR, 3'UTR, 5'Flank, 3'Flank,
#' RNA, Intron, IGR, Targeted_Region. Any class outside this list ranks last and is
#' coloured grey. Colours follow a colourblind-safe (Okabe-Ito) palette.
#'
#' Annotations: a right-side bar shows each gene's total variant burden; a top bar
#' shows each sample's total variant burden (axis labelled in thousands). Empty cells
#' (gene not mutated in that sample) are light grey.
#'
#' Data conventions:
#' \itemize{
#'   \item "Missing" means `NA` OR empty string `""`.
#'   \item Unknown/blank gene symbols are `Hugo_Symbol` in `c(".", "", "Unknown")`;
#'     these are excluded from the oncoplot.
#'   \item Works on ANY MAF-shaped table; it makes no assumption about prior
#'     filtering. It is commonly run on [gvr_filter()] output.
#' }
#'
#' @param maf A `data.table`/`data.frame` MAF from [read.gvr()] or [gvr_filter()].
#'   Required columns: `Hugo_Symbol`, `Variant_Classification`.
#' @param top_n Integer; number of genes (rows) shown, ranked by number of distinct
#'   samples mutated then by variant count. Default `20`.
#' @param sample_col Name of the per-sample column. Default `"Tumor_Sample_Barcode"`.
#'   If absent, all rows are pooled into a single sample `"All"` (with a warning).
#' @param out_dir Output directory for the PNG. Created if it does not exist.
#'   Default `"."` (current working directory).
#' @param file_prefix Base filename for the written PNG. Default `"gvr_plot"`;
#'   the file is written as `<file_prefix>.png` (fixed name, no timestamp), e.g.
#'   `gvr_plot.png`. An existing file at that path is overwritten (a message is
#'   emitted when `verbose = TRUE`).
#' @param verbose Logical; if `TRUE` (default) print the path of the file written.
#'
#' @return Invisibly, the path of the written PNG (character), or `NA_character_` if
#'   the plot was skipped (\pkg{ComplexHeatmap} not installed, or no known-gene
#'   variants present).
#'
#' @section Dependencies:
#' Requires \pkg{ComplexHeatmap} (a \pkg{Bioconductor} package, listed in `Suggests`).
#' If it is not installed the oncoplot is skipped with a warning and `NA_character_`
#' is returned.
#'
#' @seealso [gvr_summary()] for the tabular summary, [read.gvr()] to build the MAF,
#'   [gvr_filter()] to filter it before plotting.
#' @family germlinevaR
#' @author germlinevaR authors
#'
#' @examples
#' \dontrun{
#' maf <- read.gvr("/path/to/vcf_folder")
#'
#' ## write a top-20 oncoplot to the current directory
#' p <- gvr_plot(maf)
#' p                                  # path to the PNG
#'
#' ## smaller oncoplot of filtered hits, into a results folder
#' gvr_plot(gvr_filter(maf), top_n = 15, out_dir = "results/plots")
#' }
#'
#' @importFrom data.table as.data.table data.table setorder uniqueN :=
#' @importFrom grDevices png dev.off
#' @importFrom grid gpar grid.rect unit
#' @importFrom utils head
#' @export
gvr_plot <- function(maf,
                         top_n       = 20,
                         sample_col  = "Tumor_Sample_Barcode",
                         out_dir     = ".",
                         file_prefix = "gvr_plot",
                         verbose     = TRUE) {

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("gvr_plot requires the 'data.table' package.")
  }
  dt <- data.table::as.data.table(maf)

  # --- Soft guard for IMPACT column (used by top annotation) ----------------
  has_impact <- "IMPACT" %in% names(dt)
  if (!has_impact)
    warning("gvr_plot: 'IMPACT' column not found; falling back to total-burden bar.")

  .is_missing <- function(v) is.na(v) | v == ""
  UNKNOWN_GENE <- c(".", "", "Unknown")

  # --- Variant_Classification severity order (maftools-style, high -> low). ----
  #     Used to collapse multi-class gene x sample cells to a single most-severe
  #     class. Any class not listed sorts LAST and renders grey.
  GVR_SEVERITY <- c("Translation_Start_Site", "Nonsense_Mutation", "Nonstop_Mutation",
                    "Splice_Site", "Frame_Shift_Del", "Frame_Shift_Ins",
                    "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation",
                    "Splice_Region", "Protein_altering_variant", "Silent",
                    "5'UTR", "3'UTR", "5'Flank", "3'Flank", "RNA", "Intron",
                    "IGR", "Targeted_Region")
  .sev_rank <- function(x) { r <- match(x, GVR_SEVERITY); r[is.na(r)] <- length(GVR_SEVERITY) + 1L; r }
  # Colorblind-safe palette (Okabe-Ito + extensions), keyed by class; "Other" = grey.
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

  # --- Resolve sample column ---------------------------------------------------
  if (!sample_col %in% names(dt)) {
    warning(sprintf("gvr_plot: sample column '%s' not found; pooling all rows into 'All'.",
                    sample_col))
    dt[, .__sample__ := "All"]
  } else {
    dt[, .__sample__ := as.character(get(sample_col))]
    dt[.is_missing(.__sample__), .__sample__ := "NA_sample"]
  }
  samples <- sort(unique(dt$.__sample__))

  # --- Column-existence guard --------------------------------------------------
  req <- c("Hugo_Symbol", "Variant_Classification")
  miss_req <- req[!req %in% names(dt)]
  if (length(miss_req) > 0) {
    stop(sprintf("gvr_plot: required column(s) not found: %s",
                 paste(miss_req, collapse = ", ")))
  }

  # --- FUSE-safe PNG writer: render via `draw_fun(tmp)` to a tempdir path, then ---
  #     shell-cp to the final out_dir path (S3-backed mounts can 0-byte direct
  #     random-access writes). Returns the final path on success, NA on failure.
  .fuse_save_png <- function(final_path, draw_fun) {
    tmp <- file.path(tempdir(), basename(final_path))
    ok <- tryCatch({ draw_fun(tmp); file.exists(tmp) && file.info(tmp)$size > 0 },
                   error = function(e) {
                     warning(sprintf("gvr_plot: plot render failed (%s): %s",
                                     basename(final_path), conditionMessage(e))); FALSE })
    if (!ok) return(NA_character_)
    system2("cp", c("-f", shQuote(tmp), shQuote(final_path)))
    if (!file.exists(final_path) || file.info(final_path)$size == 0) {
      warning(sprintf("gvr_plot: copy to '%s' failed; left at '%s'.", final_path, tmp))
      return(tmp)
    }
    final_path
  }

  # --- Oncoplot builder (ComplexHeatmap): top-N genes x samples, each cell the ---
  #     single MOST-SEVERE Variant_Classification. Returns final path or NA.
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    warning("gvr_plot: 'ComplexHeatmap' not installed; skipping oncoplot."); return(invisible(NA_character_))
  }
  m <- dt[!(Hugo_Symbol %in% UNKNOWN_GENE)]
  if (nrow(m) == 0L) { warning("gvr_plot: no known-gene variants; skipping oncoplot."); return(invisible(NA_character_)) }
  gstat <- m[, .(n_var = .N, n_samp = data.table::uniqueN(.__sample__)), by = Hugo_Symbol]
  data.table::setorder(gstat, -n_samp, -n_var)
  top_g <- utils::head(gstat$Hugo_Symbol, top_n)
  sub <- m[Hugo_Symbol %in% top_g]
  sub[, .sr := .sev_rank(Variant_Classification)]
  cell <- sub[, .SD[which.min(.sr)], by = .(Hugo_Symbol, .__sample__)]
  mat <- matrix("", nrow = length(top_g), ncol = length(samples),
                dimnames = list(top_g, samples))
  for (i in seq_len(nrow(cell)))
    mat[cell$Hugo_Symbol[i], cell$.__sample__[i]] <- cell$Variant_Classification[i]
  # gene order: most-mutated at top (already ordered in top_g)
  mat <- mat[rev(top_g), , drop = FALSE]
  classes_present <- setdiff(unique(as.vector(mat)), "")
  col_map <- GVR_CLASS_COLORS[classes_present]
  col_map[is.na(col_map)] <- GVR_CLASS_COLORS[["Other"]]
  names(col_map) <- classes_present
  # per-gene total variant burden (right annotation) and per-sample burden (top)
  gene_burden <- m[Hugo_Symbol %in% top_g, .N, by = Hugo_Symbol]
  gb <- gene_burden$N[match(rownames(mat), gene_burden$Hugo_Symbol)]; gb[is.na(gb)] <- 0
  samp_burden <- vapply(samples, function(s) sum(m$.__sample__ == s), integer(1))

  # IMPACT palette and per-sample counts matrix (rows = levels, cols = samples)
  IMPACT_LEVELS  <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
  IMPACT_COLORS  <- c(HIGH = "#D55E00", MODERATE = "#E69F00",
                       LOW  = "#009E73", MODIFIER  = "#BBBBBB")
  if (has_impact) {
    imp_mat <- do.call(rbind, lapply(IMPACT_LEVELS, function(lv) {
      vapply(samples, function(s)
        sum(m$.__sample__ == s & !is.na(m$IMPACT) & m$IMPACT == lv),
        integer(1))
    }))
    rownames(imp_mat) <- IMPACT_LEVELS
  }

  # Top annotation: stacked IMPACT bar (or fallback total-burden bar)
  ta <- if (has_impact) {
    ComplexHeatmap::HeatmapAnnotation(
      `Variant impact` = ComplexHeatmap::anno_barplot(
        t(imp_mat), border = FALSE, beside = FALSE,
        gp = grid::gpar(fill = IMPACT_COLORS[IMPACT_LEVELS], col = NA),
        axis_param = list(
          at     = pretty(c(0, colSums(imp_mat)), n = 3),
          labels = paste0(round(pretty(c(0, colSums(imp_mat)), n = 3) / 1000), "k"),
          gp     = grid::gpar(fontsize = 7))),
      height = grid::unit(1.6, "cm"),
      annotation_name_gp = grid::gpar(fontsize = 9))
  } else {
    # fallback: original total-burden bar
    ComplexHeatmap::HeatmapAnnotation(
      `Burden` = ComplexHeatmap::anno_barplot(
        samp_burden, border = FALSE,
        gp = grid::gpar(fill = "#0279EE", col = NA),
        axis_param = list(
          at     = pretty(c(0, samp_burden), n = 3),
          labels = paste0(round(pretty(c(0, samp_burden), n = 3) / 1000), "k"),
          gp     = grid::gpar(fontsize = 7))),
      height = grid::unit(1.6, "cm"),
      annotation_name_gp = grid::gpar(fontsize = 9))
  }
  cell_fun <- function(j, i, x, y, width, height, fill) {
    v <- mat[i, j]
    grid::grid.rect(x, y, width = width * 0.95, height = height * 0.95,
                    gp = grid::gpar(fill = if (v == "") "#F2F2F2" else col_map[[v]], col = "white", lwd = 1))
  }
  ht <- ComplexHeatmap::Heatmap(
    mat, name = "Most severe\nclass", col = col_map, rect_gp = grid::gpar(type = "none"),
    cell_fun = cell_fun, na_col = "#F2F2F2",
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_heatmap_legend = TRUE, row_names_side = "left", column_names_side = "top",
    column_names_rot = 45, right_annotation = ra, top_annotation = ta,
    column_title = sprintf("Top %d genes \u00d7 %d sample(s) \u2014 cells show most-severe class",
                           length(top_g), length(samples)),
    column_title_gp = grid::gpar(fontsize = 11))

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # Fixed output name (no timestamp), mirroring gvr_summary; re-runs overwrite in
  # place and announce the overwrite when verbose.
  final_path <- file.path(out_dir, sprintf("%s.png", file_prefix))
  if (file.exists(final_path) && isTRUE(verbose))
    message(sprintf("gvr_plot: overwriting existing %s", final_path))
  # IMPACT legend (anno_barplot does not auto-generate one)
  impact_lgd <- if (has_impact)
    ComplexHeatmap::Legend(
      labels    = IMPACT_LEVELS,
      title     = "IMPACT",
      legend_gp = grid::gpar(fill = IMPACT_COLORS[IMPACT_LEVELS]))
  else NULL

  path <- .fuse_save_png(final_path, function(tmp) {
    grDevices::png(tmp, width = max(1100, 360 + 150 * length(samples)),
                   height = max(720, 110 + 34 * length(top_g)), res = 150)
    ComplexHeatmap::draw(ht, heatmap_legend_side = "right", merge_legend = TRUE,
                         annotation_legend_list = if (!is.null(impact_lgd)) list(impact_lgd) else list(),
                         padding = grid::unit(c(2, 6, 2, 2), "mm"))
    grDevices::dev.off()
  })
  if (!is.na(path) && isTRUE(verbose)) message(sprintf("gvr_plot: written %s", path))
  invisible(path)
}

# NOTE: globalVariables() declarations for this package are consolidated in
# R/globals.R (one package-scoped block covering all functions).
