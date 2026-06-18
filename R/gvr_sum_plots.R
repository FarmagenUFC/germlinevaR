# gvr_sum_plots.R
# ============================================================================
# germlinevaR: produce gvr_summary() plots as individual image files + panels
#
# Companion function to gvr_summary(). Recomputes the same cohort summaries
# gvr_summary() builds internally, then writes each chart as a standalone
# image plus a multi-panel composite, all in a new folder.
#
# Public function: gvr_sum_plots(maf, ...) -- same MAF requirements as
# gvr_summary(). Returns the output folder path invisibly.
#
# Engine: ggplot2 + ggsave(); panels assembled with patchwork when available,
# fallback to gridExtra::grid.arrange. Default format is "png"; any device
# name accepted by ggsave() is permitted.
# ----------------------------------------------------------------------------

#' Export gvr_summary() plots as standalone image files
#'
#' Produces the same cohort-level and per-sample plots that [gvr_summary()]
#' renders in its HTML/PDF dashboard, but as individual image files and as
#' multi-panel composites, written into a new folder.
#'
#' @details
#' The function recomputes the summary sections from the MAF (it does not
#' depend on a prior call to [gvr_summary()]), then writes:
#'
#' Cohort-level plots (always produced):
#' \itemize{
#'   \item `top_genes.<ext>`: top genes by variant count
#'   \item `variant_classification.<ext>`: top variant classifications
#'   \item `impact.<ext>`: VEP IMPACT severity
#'   \item `top_variants.<ext>`: top variants by `dbSNP_RS` recurrence
#'     (omitted if the column is missing or all empty)
#' }
#'
#' Per-sample plots (when `per_sample = TRUE`, one file per sample):
#' \itemize{
#'   \item `per_sample/top_genes__<sample>.<ext>`
#' }
#'
#' Panels (when `panel = TRUE`):
#' \itemize{
#'   \item `panel_cohort.<ext>`: 2x2 grid of the four cohort plots
#'   \item `panel_per_sample.<ext>`: grid of every per-sample plot
#'     (omitted if `n_samples == 1` since it would duplicate the
#'     standalone per-sample image)
#' }
#'
#' All files are written under `out_dir/folder_name/`. Existing files of the
#' same name are overwritten; other files in the folder are left alone.
#'
#' @param maf data.table or data.frame produced by [read.gvr()] /
#'   [read.gvr.snpeff()] / [read.gvr.dual()]. Required columns:
#'   `Hugo_Symbol`, `Variant_Classification`, `Variant_Type`, `IMPACT`,
#'   `CLIN_SIG`; `dbSNP_RS` is optional.
#' @param out_dir parent directory under which the output folder is created.
#'   Default `"."` (current working directory).
#' @param folder_name name of the new folder created under `out_dir`.
#'   Default `"gvr_sum_plots"`.
#' @param format one of `"png"` (default), `"pdf"`, `"svg"`, `"jpeg"`,
#'   `"tiff"`, `"bmp"`, `"eps"`, `"ps"`, `"tex"`, `"wmf"`. Mirrors the device
#'   list of [ggplot2::ggsave()].
#' @param width plot width in inches (passed to [ggplot2::ggsave()]).
#'   Default 7.
#' @param height plot height in inches. Default 5.
#' @param dpi resolution for raster formats (ignored for vector). Default 300.
#' @param sample_col column holding the sample identifier. Default
#'   `"Tumor_Sample_Barcode"`.
#' @param top_n_genes number of top genes to keep in `top_genes` and per-sample
#'   plots. Default 20.
#' @param top_n_variants number of top variants to keep in `top_variants`.
#'   Default 20.
#' @param per_sample logical; if `FALSE`, skip per-sample plots and the
#'   per-sample panel. Default `TRUE`.
#' @param panel logical; if `FALSE`, skip both panel images. Default `TRUE`.
#' @param verbose logical; print progress messages. Default `TRUE`.
#'
#' @return The output folder path, invisibly.
#'
#' @seealso [gvr_summary()] for the dashboard view, [read.gvr()] for the
#'   reader.
#'
#' @family germlinevaR plotting
#'
#' @author Thiago Loreto Matos
#'
#' @examples
#' \dontrun{
#'   maf <- read.gvr("/path/to/vcf/")
#'   gvr_sum_plots(maf, out_dir = "results", folder_name = "S6_plots",
#'                 format = "pdf")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_col coord_flip scale_fill_manual labs theme_minimal theme element_text element_blank facet_wrap as_labeller scale_y_continuous position_dodge ggsave
#' @importFrom data.table as.data.table data.table setorder
#' @export
gvr_sum_plots <- function(maf,
                          out_dir        = ".",
                          folder_name    = "gvr_sum_plots",
                          format         = "png",
                          width          = 7,
                          height         = 5,
                          dpi            = 300,
                          sample_col     = "Tumor_Sample_Barcode",
                          top_n_genes    = 20,
                          top_n_variants = 20,
                          per_sample     = TRUE,
                          panel          = TRUE,
                          verbose        = TRUE) {

  # ---------- arg validation ----------------------------------------------
  valid_fmts <- c("png", "pdf", "svg", "jpeg", "tiff", "bmp",
                  "eps", "ps", "tex", "wmf")
  format <- tolower(format[[1L]])
  if (!format %in% valid_fmts) {
    stop(sprintf(
      "gvr_sum_plots: 'format' must be one of: %s. Got '%s'.",
      paste(valid_fmts, collapse = ", "), format), call. = FALSE)
  }
  if (!requireNamespace("ggplot2",   quietly = TRUE))
    stop("gvr_sum_plots requires 'ggplot2'.", call. = FALSE)
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("gvr_sum_plots requires 'data.table'.", call. = FALSE)
  # Format-specific package requirements (ggsave defers these to dev() time;
  # we surface them up-front so the user fixes them before any plot is built).
  fmt_pkg <- switch(format,
                    svg  = "svglite",
                    wmf  = "ragg",
                    NULL)
  if (!is.null(fmt_pkg) && !requireNamespace(fmt_pkg, quietly = TRUE)) {
    stop(sprintf(
      "gvr_sum_plots: format '%s' requires the '%s' package. Install it with install.packages('%s'), or pick another format (e.g. 'pdf' for vector or 'png' for raster).",
      format, fmt_pkg, fmt_pkg), call. = FALSE)
  }
  if (!is.numeric(width)  || width  <= 0) stop("gvr_sum_plots: 'width' must be a positive number.",  call. = FALSE)
  if (!is.numeric(height) || height <= 0) stop("gvr_sum_plots: 'height' must be a positive number.", call. = FALSE)

  dt <- data.table::as.data.table(maf)
  if (nrow(dt) == 0L)
    stop("gvr_sum_plots: 'maf' has zero rows; nothing to plot.", call. = FALSE)

  # ---------- output folder ------------------------------------------------
  out_path <- file.path(out_dir, folder_name)
  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
  if (per_sample) dir.create(file.path(out_path, "per_sample"),
                             recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(verbose))
    message(sprintf("gvr_sum_plots: writing to '%s' (format: %s)", out_path, format))

  # ---------- build the summary sections (verbatim of gvr_summary block) ---
  built <- .gvr_sp_build_sections(dt,
                                  sample_col     = sample_col,
                                  top_n_genes    = top_n_genes,
                                  top_n_variants = top_n_variants,
                                  verbose        = verbose)
  sections <- built$sections
  samples  <- built$samples

  # ---------- color palette mirrors gvr_summary ---------------------------
  PHYLO_BLUE  <- "#0279EE"
  PHYLO_GREEN <- "#75A025"
  STONE       <- "#ECE9E2"
  CREAM       <- "#FAF9F3"
  INK         <- "#000000"
  ORANGE      <- "#E69F00"
  VERMIL      <- "#D55E00"
  YELLOW      <- "#E9ED4C"
  PINK        <- "#FD9BED"
  PHYLO_PALETTE <- c(PHYLO_BLUE, PHYLO_GREEN, ORANGE, VERMIL, YELLOW, PINK)
  fill_vals <- stats::setNames(
    PHYLO_PALETTE[((seq_along(samples) - 1L) %% length(PHYLO_PALETTE)) + 1L],
    samples)
  FACET_THRESHOLD <- 6L

  # Compact facet labels (strip common >= 3-char prefix)
  .lab_fun <- local({
    sn <- as.character(samples); pre <- ""
    if (length(sn) > 1L) {
      mn <- min(nchar(sn)); i <- 0L
      while (i < mn && length(unique(substr(sn, 1L, i + 1L))) == 1L) i <- i + 1L
      if (i >= 3L) pre <- substr(sn[1], 1L, i)
    }
    function(x) if (nzchar(pre)) sub(paste0("^", pre), "", x) else x
  })

  # ---------- build the four cohort ggplots --------------------------------
  cohort_specs <- list(
    top_genes = list(
      sec    = sections$top_genes,
      catcol = "Hugo_Symbol",
      top    = top_n_genes,
      title  = sprintf("Top genes by variant count (top %d)", top_n_genes),
      file   = "top_genes"),
    variant_classification = list(
      sec    = sections$variant_classification,
      catcol = "Variant_Classification",
      top    = 10L,
      title  = "Variant classification (top 10)",
      file   = "variant_classification"),
    impact = list(
      sec    = sections$impact,
      catcol = "IMPACT",
      top    = NULL,
      title  = "Functional impact (VEP IMPACT)",
      file   = "impact"))

  if (!is.null(sections$top_variants) && nrow(sections$top_variants) > 0L) {
    cohort_specs$top_variants <- list(
      sec    = sections$top_variants,
      catcol = "dbSNP_RS",
      top    = top_n_variants,
      title  = sprintf("Top variants by recurrence (top %d)", top_n_variants),
      file   = "top_variants")
  } else if (isTRUE(verbose)) {
    message("  Note: top_variants section empty/missing -- skipping that plot.")
  }

  cohort_plots <- lapply(cohort_specs, function(s)
    .gvr_sp_bar_gg(s$sec, s$catcol, samples = samples, top = s$top,
                   title = s$title, fill_vals = fill_vals,
                   facet_threshold = FACET_THRESHOLD,
                   lab_fun = .lab_fun, phylo_blue = PHYLO_BLUE))

  # ---------- save cohort plots --------------------------------------------
  for (nm in names(cohort_specs)) {
    fp <- file.path(out_path,
                    sprintf("%s.%s", cohort_specs[[nm]]$file, format))
    .gvr_sp_save_one(cohort_plots[[nm]], fp, format, width, height, dpi)
    if (isTRUE(verbose)) message("  wrote ", basename(fp))
  }

  # ---------- per-sample plots ---------------------------------------------
  ps_plots <- list()
  if (per_sample) {
    tgps <- sections$top_genes_per_sample
    for (sm in names(tgps)) {
      sm_dt <- tgps[[sm]]
      if (!nrow(sm_dt)) next
      gg <- .gvr_sp_bar_gg_single(sm_dt, cat_col = "Hugo_Symbol",
                                  val_col = sm, top = 10L,
                                  title = sprintf("Top genes (top 10) -- %s", sm),
                                  fill_color = PHYLO_BLUE,
                                  phylo_blue = PHYLO_BLUE)
      ps_plots[[sm]] <- gg
      safe_sm <- .gvr_sp_safe_filename(sm)
      fp <- file.path(out_path, "per_sample",
                      sprintf("top_genes__%s.%s", safe_sm, format))
      .gvr_sp_save_one(gg, fp, format, width, height, dpi)
    }
    if (isTRUE(verbose))
      message(sprintf("  wrote %d per-sample plot(s) to '%s/per_sample/'",
                      length(ps_plots), folder_name))
  }

  # ---------- panels -------------------------------------------------------
  if (panel) {
    # Cohort panel
    panel_w <- width  * 2
    panel_h <- height * ceiling(length(cohort_plots) / 2L)
    panel_cohort <- .gvr_sp_panel(cohort_plots, ncol = 2L)
    fp <- file.path(out_path, sprintf("panel_cohort.%s", format))
    .gvr_sp_save_one(panel_cohort, fp, format, panel_w, panel_h, dpi)
    if (isTRUE(verbose)) message("  wrote ", basename(fp))

    # Per-sample panel (skip if only 1 sample: it would duplicate the single
    # per-sample plot)
    if (per_sample && length(ps_plots) >= 2L) {
      ncol_ps <- ceiling(sqrt(length(ps_plots)))
      nrow_ps <- ceiling(length(ps_plots) / ncol_ps)
      panel_ps <- .gvr_sp_panel(ps_plots, ncol = ncol_ps)
      psp_w <- width  * ncol_ps
      psp_h <- height * nrow_ps
      fp <- file.path(out_path, sprintf("panel_per_sample.%s", format))
      .gvr_sp_save_one(panel_ps, fp, format, psp_w, psp_h, dpi)
      if (isTRUE(verbose)) message("  wrote ", basename(fp))
    }
  }

  if (isTRUE(verbose)) message("gvr_sum_plots: done.")
  invisible(out_path)
}

# ============================================================================
# Private helpers
# ============================================================================

# ---------------------------------------------------------------------------
# .gvr_sp_build_sections()
#   Builds the same `sections` list gvr_summary() builds internally. Copied
#   verbatim from gvr_summary.R (lines ~219-399) with the verbose digest
#   trimmed (gvr_sum_plots() owns its own progress messages).
# ---------------------------------------------------------------------------
.gvr_sp_build_sections <- function(dt, sample_col, top_n_genes, top_n_variants,
                                   verbose = FALSE) {
  .is_missing <- function(v) is.na(v) | v == ""
  UNKNOWN_GENE <- c(".", "", "Unknown")

  if (!sample_col %in% names(dt)) {
    warning(sprintf("gvr_sum_plots: sample column '%s' not found; pooling all rows into 'All'.",
                    sample_col))
    dt[, .__sample__ := "All"]
  } else {
    dt[, .__sample__ := as.character(get(sample_col))]
    dt[.is_missing(.__sample__), .__sample__ := "NA_sample"]
  }
  samples <- sort(unique(dt$.__sample__))

  req <- c("Hugo_Symbol", "Variant_Classification", "Variant_Type", "IMPACT", "CLIN_SIG")
  miss_req <- req[!req %in% names(dt)]
  if (length(miss_req) > 0) {
    stop(sprintf("gvr_sum_plots: required column(s) not found: %s",
                 paste(miss_req, collapse = ", ")), call. = FALSE)
  }

  .count_by_sample <- function(valuevec, samplevec, category_name,
                               order_levels = NULL, drop_values = NULL) {
    keep <- rep(TRUE, length(valuevec))
    if (!is.null(drop_values)) keep <- !(valuevec %in% drop_values)
    v <- valuevec[keep]; s <- samplevec[keep]
    if (length(v) == 0L) {
      out <- data.table::data.table(X = character(0))
      data.table::setnames(out, "X", category_name)
      for (sm in samples) out[, (sm) := integer(0)]
      out[, Total := integer(0)]
      return(out[])
    }
    tab <- table(factor(v), factor(s, levels = samples))
    m <- data.table::as.data.table(unclass(tab), keep.rownames = TRUE)
    data.table::setnames(m, "rn", category_name)
    for (sm in samples) if (!sm %in% names(m)) m[, (sm) := 0L]
    data.table::setcolorder(m, c(category_name, samples))
    m[, Total := rowSums(as.matrix(.SD)), .SDcols = samples]
    if (!is.null(order_levels)) {
      m <- m[match(order_levels, m[[category_name]])]
      m <- m[!is.na(m[[category_name]])]
    } else {
      data.table::setorder(m, -Total)
    }
    m[]
  }

  # Top genes
  gene_tab <- .count_by_sample(dt$Hugo_Symbol, dt$.__sample__, "Hugo_Symbol",
                               drop_values = UNKNOWN_GENE)
  top_genes <- utils::head(gene_tab, top_n_genes)

  # Per-sample top genes
  top_genes_per_sample <- list()
  for (sm in samples) {
    dt_sm <- dt[.__sample__ == sm & !(Hugo_Symbol %in% UNKNOWN_GENE)]
    if (nrow(dt_sm) > 0L) {
      sm_gene_tab <- dt_sm[, .(N = .N), by = Hugo_Symbol]
      data.table::setorder(sm_gene_tab, -N)
      sm_gene_tab <- utils::head(sm_gene_tab, top_n_genes)
      data.table::setnames(sm_gene_tab, "N", sm)
      top_genes_per_sample[[sm]] <- sm_gene_tab
    }
  }

  # Variant classification
  variant_classification <- .count_by_sample(dt$Variant_Classification, dt$.__sample__,
                                              "Variant_Classification")

  # Top variants by dbSNP_RS
  top_variants <- NULL
  has_dbsnp <- "dbSNP_RS" %in% names(dt)
  if (has_dbsnp) {
    rs <- dt$dbSNP_RS
    rs_miss <- .is_missing(rs) | rs == "novel" | tolower(rs) == "unknown"
    dt_rs <- dt[!rs_miss]
    if (nrow(dt_rs) > 0L) {
      dt_rs[, .__rs__ := dbSNP_RS]
      rs_tab <- dt_rs[, .(Hugo_Symbol = {
        v <- Hugo_Symbol[!(Hugo_Symbol %in% UNKNOWN_GENE) & !.is_missing(Hugo_Symbol)]
        if (length(v)) v[1L] else ""
      }, Total = .N), by = .__rs__]
      dt_rs_cast <- dt_rs[, .(.__rs__, .__sample__, .__n__ = 1L)]
      samp_counts <- data.table::dcast(
        dt_rs_cast, .__rs__ ~ .__sample__,
        fun.aggregate = sum, value.var = ".__n__", fill = 0L)
      for (sm in samples) if (!sm %in% names(samp_counts)) samp_counts[, (sm) := 0L]
      rs_tab <- samp_counts[rs_tab, on = ".__rs__"]
      data.table::setcolorder(rs_tab, c(".__rs__", "Hugo_Symbol", samples, "Total"))
      data.table::setorder(rs_tab, -Total)
      data.table::setnames(rs_tab, ".__rs__", "dbSNP_RS")
      top_variants <- utils::head(rs_tab, top_n_variants)
    }
  }

  # Impact
  impact_order <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
  present_impact <- c(intersect(impact_order, unique(dt$IMPACT)),
                      setdiff(unique(dt$IMPACT[!.is_missing(dt$IMPACT)]), impact_order))
  impact <- .count_by_sample(dt$IMPACT, dt$.__sample__, "IMPACT",
                             order_levels = present_impact)

  sections <- list(top_genes = top_genes,
                   top_genes_per_sample = top_genes_per_sample,
                   variant_classification = variant_classification,
                   impact = impact)
  if (!is.null(top_variants)) sections$top_variants <- top_variants

  list(sections = sections, samples = samples)
}

# ---------------------------------------------------------------------------
# .gvr_sp_bar_gg()
#   Returns a ggplot object (not a grob) for a cohort-level bar chart.
#   Switches between grouped (<=6 samples) and faceted (>6) layout, same as
#   .bar_grob() in gvr_summary.R but minus the trailing ggplotGrob().
# ---------------------------------------------------------------------------
.gvr_sp_bar_gg <- function(dt, cat_col, samples, top, title, fill_vals,
                           facet_threshold, lab_fun, phylo_blue) {
  d <- as.data.frame(dt, stringsAsFactors = FALSE)
  d <- d[d[[cat_col]] != "" & !is.na(d[[cat_col]]), , drop = FALSE]
  if (!is.null(top) && nrow(d) > top)
    d <- d[order(-d$Total)[seq_len(top)], , drop = FALSE]
  long <- do.call(rbind, lapply(samples, function(s)
    data.frame(Category = d[[cat_col]], Sample = s, n = d[[s]],
               stringsAsFactors = FALSE)))
  long$Category <- factor(long$Category, levels = d[[cat_col]][order(d$Total)])
  many <- length(samples) > facet_threshold
  gg <- ggplot2::ggplot(long, ggplot2::aes(x = Category, y = n, fill = Sample)) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = fill_vals) +
    ggplot2::labs(title = title, x = NULL, y = NULL, fill = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                       colour = phylo_blue,
                                                       size = 13),
                   panel.grid.major.y = ggplot2::element_blank())
  if (many) {
    gg <- gg + ggplot2::geom_col(show.legend = FALSE) +
      ggplot2::facet_wrap(~ Sample, ncol = ceiling(sqrt(length(samples))),
                          labeller = ggplot2::as_labeller(lab_fun)) +
      ggplot2::scale_y_continuous(
        breaks = scales::breaks_extended(3),
        labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
      ggplot2::theme(legend.position = "none",
                     strip.text = ggplot2::element_text(size = 7, face = "bold"),
                     axis.text  = ggplot2::element_text(size = 5.5),
                     panel.spacing = grid::unit(4, "pt"))
  } else {
    gg <- gg + ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75),
                                  width = 0.68) +
      ggplot2::scale_y_continuous(
        breaks = scales::breaks_extended(4),
        labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
      ggplot2::theme(legend.position = "top", legend.justification = "left",
                     legend.key.size = grid::unit(10, "pt"),
                     legend.text = ggplot2::element_text(size = 8),
                     axis.text = ggplot2::element_text(size = 9))
  }
  gg
}

# ---------------------------------------------------------------------------
# .gvr_sp_bar_gg_single()
#   Single-sample horizontal bar chart (per-sample top genes). ggplot2
#   translation of gvr_summary()'s .plt_bar_single() plotly helper.
# ---------------------------------------------------------------------------
.gvr_sp_bar_gg_single <- function(dt, cat_col, val_col, top = 10L, title = NULL,
                                  fill_color, phylo_blue) {
  d <- as.data.frame(dt, stringsAsFactors = FALSE)
  d <- d[d[[cat_col]] != "" & !is.na(d[[cat_col]]), , drop = FALSE]
  if (!is.null(top) && nrow(d) > top)
    d <- d[order(-d[[val_col]])[seq_len(top)], , drop = FALSE]
  d[[cat_col]] <- factor(d[[cat_col]], levels = d[[cat_col]][order(d[[val_col]])])
  ggplot2::ggplot(d, ggplot2::aes(x = .data[[cat_col]], y = .data[[val_col]])) +
    ggplot2::geom_col(fill = fill_color, width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::scale_y_continuous(
      breaks = scales::breaks_extended(4),
      labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                       colour = phylo_blue,
                                                       size = 12),
                   panel.grid.major.y = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size = 9))
}

# ---------------------------------------------------------------------------
# .gvr_sp_save_one()
#   Thin wrapper around ggsave(). Passes dpi only for raster formats; for
#   vector formats ggsave ignores it anyway but we keep the call clean.
# ---------------------------------------------------------------------------
.gvr_sp_save_one <- function(gg, path, format, width, height, dpi) {
  raster <- format %in% c("png", "jpeg", "tiff", "bmp")
  args <- list(filename = path, plot = gg, device = format,
               width = width, height = height, units = "in")
  if (raster) args$dpi <- dpi
  do.call(ggplot2::ggsave, args)
  invisible(path)
}

# ---------------------------------------------------------------------------
# .gvr_sp_panel()
#   Combine a list of ggplot objects into one composite plot. Uses
#   patchwork::wrap_plots when available, otherwise gridExtra::grid.arrange.
#   Returns either a ggplot/patchwork object (saveable directly via ggsave)
#   or a grob (saveable via ggsave when wrapped).
# ---------------------------------------------------------------------------
.gvr_sp_panel <- function(plots, ncol) {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(patchwork::wrap_plots(plots, ncol = ncol))
  }
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    g <- gridExtra::arrangeGrob(grobs = plots, ncol = ncol)
    # Wrap the grob so ggsave can serialize it as a single plot
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::annotation_custom(g))
  }
  stop("gvr_sum_plots: panel assembly needs 'patchwork' or 'gridExtra'.",
       call. = FALSE)
}

# ---------------------------------------------------------------------------
# .gvr_sp_safe_filename()
#   Sanitize sample names for use inside filenames. Anything outside
#   [A-Za-z0-9._-] becomes '_'.
# ---------------------------------------------------------------------------
.gvr_sp_safe_filename <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)
