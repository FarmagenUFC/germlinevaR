#' Multi-section summary of a germline MAF (read.gvr / gvr_filter output)
#'
#' @description
#' Produces a multi-section overview of a MAF-style table - the output of
#' [read.gvr()], or of [gvr_filter()] - covering variant burden, affected genes,
#' functional classes, clinical significance and predicted impact. Every section is
#' returned as a tidy `data.table` with one column per sample plus a `Total` column.
#' Optionally writes a multi-sheet Excel workbook and/or a multi-page PDF report
#' (both into a `gvr_summary/` subfolder of `out_dir`). For a cohort oncoplot, see
#' [gvr_oncoplot()].
#'
#' @details
#' Sections returned (as a named list of `data.table`s):
#' \itemize{
#'   \item `overview` - cohort-level counts: total variants, variants per sample,
#'     distinct genes affected, and variants with no gene symbol.
#'   \item `top_genes` - the `top_n_genes` genes with the most variants
#'     (per-sample + `Total`); unknown/blank genes excluded.
#'   \item `variant_classification` - counts per `Variant_Classification` (maftools
#'     functional class), per-sample + `Total`, sorted by `Total` descending.
#'   \item `variant_type` - counts per `Variant_Type` (SNP/DEL/INS/ONP/DNP/TNP).
#'   \item `clin_sig` - counts per `CLIN_SIG` token (ClinVar categories). `CLIN_SIG`
#'     is split on `&` and `/`, so a variant annotated `"pathogenic&benign"`
#'     increments BOTH categories; category counts therefore sum to \eqn{\ge} the
#'     number of variants. A `missing/unclassified` row counts NA/"" `CLIN_SIG`.
#'   \item `impact` - counts per VEP `IMPACT` (HIGH/MODERATE/LOW/MODIFIER), in
#'     severity order rather than count order.
#' }
#'
#' The section tables are the core return value. By default the function also writes
#' an Excel workbook (`save_excel = TRUE`) and a PDF report (`save_pdf = TRUE`) into
#' `out_dir/gvr_summary/`; set either to `FALSE` to skip it. The cohort oncoplot lives
#' in [gvr_oncoplot()].
#'
#' Data conventions:
#' \itemize{
#'   \item "Missing" means `NA` OR empty string `""`.
#'   \item Unknown/blank gene symbols are `Hugo_Symbol` in `c(".", "", "Unknown")`;
#'     these are excluded from the distinct-gene tally and from `top_genes`, but their
#'     variants are still counted in the totals (and reported as "variants with no
#'     gene symbol").
#'   \item Works on ANY MAF-shaped table; it makes no assumption about prior
#'     filtering. It is commonly run on [gvr_filter()] output to summarise the
#'     retained hits.
#' }
#'
#' @param maf A `data.table`/`data.frame` MAF from [read.gvr()] or [gvr_filter()].
#'   Required columns: `Hugo_Symbol`, `Variant_Classification`, `Variant_Type`,
#'   `IMPACT`, `CLIN_SIG`.
#' @param sample_col Name of the per-sample column. Default `"Tumor_Sample_Barcode"`.
#'   If absent, all rows are pooled into a single sample `"All"` (with a warning).
#' @param top_n_genes Integer; number of genes to report in `top_genes` (by total
#'   variant count). Default `20`.
#' @param save_excel Logical; if `TRUE` (default), write a multi-sheet `.xlsx`.
#'   The workbook is written into the `gvr_summary/` subfolder of `out_dir`
#'   (see `out_dir`). Pass `FALSE` for a compute-only run.
#' @param save_pdf Logical; if `TRUE` (default), write a multi-page PDF report into
#'   the `gvr_summary/` subfolder of `out_dir`. The report has a title/metadata page,
#'   the charted sections (top genes, variant classification, predicted impact) each
#'   shown as a table together with its bar chart, and the remaining tables grouped on
#'   their own page(s). Requires \pkg{gridExtra} and \pkg{ggplot2}; if unavailable, the
#'   PDF is skipped with a warning and the sections are still returned. Pass `FALSE`
#'   for a compute-only run.
#' @param out_dir Parent output directory. All written outputs (Excel and/or PDF) are
#'   placed in a `gvr_summary/` subfolder created inside `out_dir`. The subfolder is
#'   created only when `save_excel` or `save_pdf` is `TRUE`. Default `"."` (current
#'   working directory), i.e. outputs go to `./gvr_summary/`.
#' @param file_prefix Base filename for written outputs. Default `"gvr_summary"`, giving
#'   `<file_prefix>.xlsx` and `<file_prefix>_report.pdf` (no timestamp). Filenames are
#'   fixed, so re-running into the same `out_dir` overwrites the previous files (a
#'   message is printed when `verbose = TRUE`).
#' @param verbose Logical; if `TRUE` (default) print a compact console digest and the
#'   path(s) of any file(s) written.
#'
#' @return Invisibly, a named list of `data.table`s: `overview`, `top_genes`,
#'   `variant_classification`, `variant_type`, `clin_sig`, `impact`. The return value
#'   is identical regardless of whether the Excel/PDF files are written.
#'
#' @section Dependencies:
#' Core summary uses \pkg{data.table}. The optional Excel export uses \pkg{openxlsx};
#' the optional PDF report uses \pkg{gridExtra} + \pkg{ggplot2} (rendered via the base
#' \code{grDevices::pdf} device). Each optional output degrades gracefully: if its
#' package(s) are unavailable, that output is skipped with a warning and the section
#' tables are still returned. PDF text is ASCII-only (the base `pdf()` device does not
#' embed glyphs for non-ASCII punctuation).
#'
#' @seealso [read.gvr()] to build the MAF, [gvr_filter()] to filter it before
#'   summarising, [gvr_oncoplot()] for a cohort oncoplot.
#' @family germlinevaR
#' @author germlinevaR authors
#'
#' @examples
#' \dontrun{
#' maf <- read.gvr("/path/to/vcf_folder")
#'
#' ## default: writes gvr_summary.xlsx + gvr_summary_report.pdf into ./gvr_summary/
#' s <- gvr_summary(maf)
#' s$variant_classification          # inspect a section
#' s$impact                          # HIGH -> MODIFIER, severity order
#'
#' ## compute only (no files written); print a console digest
#' s <- gvr_summary(maf, save_excel = FALSE, save_pdf = FALSE)
#'
#' ## summarise filtered hits, writing into results/summary/gvr_summary/
#' gvr_summary(gvr_filter(maf), out_dir = "results/summary")
#' }
#'
#' @importFrom data.table as.data.table data.table setnames setcolorder setorder
#'   uniqueN copy melt :=
#' @importFrom utils head
#' @importFrom openxlsx createWorkbook
#' @importFrom grDevices pdf dev.off
#' @export


gvr_summary <- function(maf,
                        sample_col     = "Tumor_Sample_Barcode",
                        top_n_genes    = 20,
                        save_excel     = TRUE,
                        save_pdf       = TRUE,
                        out_dir        = ".",
                        file_prefix    = "gvr_summary",
                        verbose        = TRUE) {

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("gvr_summary requires the 'data.table' package.")
  }
  dt <- data.table::as.data.table(maf)
  n_total <- nrow(dt)

  .is_missing <- function(v) is.na(v) | v == ""
  UNKNOWN_GENE <- c(".", "", "Unknown")
  # --- Resolve sample column ---------------------------------------------------
  if (!sample_col %in% names(dt)) {
    warning(sprintf("gvr_summary: sample column '%s' not found; pooling all rows into 'All'.",
                    sample_col))
    dt[, .__sample__ := "All"]
  } else {
    dt[, .__sample__ := as.character(get(sample_col))]
    dt[.is_missing(.__sample__), .__sample__ := "NA_sample"]
  }
  samples <- sort(unique(dt$.__sample__))

  # --- Column-existence guard for the analytic columns -------------------------
  req <- c("Hugo_Symbol", "Variant_Classification", "Variant_Type", "IMPACT", "CLIN_SIG")
  miss_req <- req[!req %in% names(dt)]
  if (length(miss_req) > 0) {
    stop(sprintf("gvr_summary: required column(s) not found: %s",
                 paste(miss_req, collapse = ", ")))
  }

  # --- Helper: counts of `valuevec` x sample as a wide data.table --------------
  #     rows = categories, cols = samples + Total. `order_levels` optionally fixes
  #     row order; otherwise sorted by Total desc.
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
    # ensure all sample columns exist (table drops empty factor combos only if level absent)
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
  # ============================================================================
  # SECTION 1: overview (variants & genes)
  # ============================================================================
  is_known_gene <- !(dt$Hugo_Symbol %in% UNKNOWN_GENE)
  ov_rows <- list()
  # total variants
  per_sample_var <- vapply(samples, function(sm) sum(dt$.__sample__ == sm), integer(1))
  ov_rows[["Total variants"]] <- c(per_sample_var, Total = n_total)
  # distinct genes (known only)
  per_sample_genes <- vapply(samples, function(sm)
    data.table::uniqueN(dt$Hugo_Symbol[dt$.__sample__ == sm & is_known_gene]), integer(1))
  ov_rows[["Distinct genes (known)"]] <- c(per_sample_genes,
                                           Total = data.table::uniqueN(dt$Hugo_Symbol[is_known_gene]))
  # variants with no gene symbol
  per_sample_nogene <- vapply(samples, function(sm)
    sum(dt$.__sample__ == sm & !is_known_gene), integer(1))
  ov_rows[["Variants with no gene symbol"]] <- c(per_sample_nogene, Total = sum(!is_known_gene))

  overview <- data.table::data.table(Metric = names(ov_rows))
  for (sm in c(samples, "Total")) {
    overview[, (sm) := vapply(ov_rows, function(r) r[[sm]], numeric(1))]
  }

  # ============================================================================
  # SECTION 1b: top genes (known only), by Total variant count
  # ============================================================================
  gene_tab <- .count_by_sample(dt$Hugo_Symbol, dt$.__sample__, "Hugo_Symbol",
                               drop_values = UNKNOWN_GENE)
  top_genes <- utils::head(gene_tab, top_n_genes)

  # ============================================================================
  # SECTION 2: functional classes
  # ============================================================================
  variant_classification <- .count_by_sample(dt$Variant_Classification, dt$.__sample__,
                                              "Variant_Classification")
  variant_type <- .count_by_sample(dt$Variant_Type, dt$.__sample__, "Variant_Type")

  # ============================================================================
  # SECTION 3: clinical categories (token-split on & and /)
  # ============================================================================
  cs   <- dt$CLIN_SIG
  miss <- .is_missing(cs)
  # expand non-missing CLIN_SIG into (token, sample) pairs
  idx_nm <- which(!miss)
  if (length(idx_nm) > 0) {
    toks_list <- strsplit(cs[idx_nm], "[&/]")
    ntok <- lengths(toks_list)
    tok_vec <- trimws(unlist(toks_list, use.names = FALSE))
    samp_vec <- rep(dt$.__sample__[idx_nm], ntok)
    # drop any empty tokens produced by stray delimiters
    keep_tok <- tok_vec != ""
    clin_counts <- .count_by_sample(tok_vec[keep_tok], samp_vec[keep_tok], "CLIN_SIG")
  } else {
    clin_counts <- .count_by_sample(character(0), character(0), "CLIN_SIG")
  }
  # append missing/unclassified row
  per_sample_miss <- vapply(samples, function(sm) sum(miss & dt$.__sample__ == sm), integer(1))
  miss_row <- data.table::as.data.table(c(list(CLIN_SIG = "missing/unclassified"),
                                          as.list(per_sample_miss),
                                          list(Total = sum(miss))))
  clin_sig <- rbind(clin_counts, miss_row, fill = TRUE)

  # ============================================================================
  # SECTION 4: impact severity (fixed severity order)
  # ============================================================================
  impact_order <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
  present_impact <- c(intersect(impact_order, unique(dt$IMPACT)),
                      setdiff(unique(dt$IMPACT[!.is_missing(dt$IMPACT)]), impact_order))
  impact <- .count_by_sample(dt$IMPACT, dt$.__sample__, "IMPACT",
                             order_levels = present_impact)

  sections <- list(overview = overview,
                   top_genes = top_genes,
                   variant_classification = variant_classification,
                   variant_type = variant_type,
                   clin_sig = clin_sig,
                   impact = impact)

  # ============================================================================
  # Verbose console digest
  # ============================================================================
  if (isTRUE(verbose)) {
    message(sprintf("gvr_summary: %d variants across %d sample(s): %s",
                    n_total, length(samples), paste(samples, collapse = ", ")))
    message(sprintf("  Distinct genes (known): %d   |  variants with no gene symbol: %d",
                    data.table::uniqueN(dt$Hugo_Symbol[is_known_gene]), sum(!is_known_gene)))
    topfc <- variant_classification[1:min(3, nrow(variant_classification))]
    message(sprintf("  Top functional classes: %s",
                    paste(sprintf("%s=%d", topfc$Variant_Classification, topfc$Total), collapse = ", ")))
    message(sprintf("  IMPACT: %s",
                    paste(sprintf("%s=%d", impact$IMPACT, impact$Total), collapse = ", ")))
    nrare <- clin_sig[CLIN_SIG != "missing/unclassified"]
    if (nrow(nrare) > 0) {
      message(sprintf("  CLIN_SIG (top tokens): %s",
                      paste(sprintf("%s=%d", utils::head(nrare$CLIN_SIG, 4),
                                    utils::head(nrare$Total, 4)), collapse = ", ")))
    }
  }

  # ============================================================================
  # Output folder: ALL written summary outputs go into <out_dir>/gvr_summary/.
  # The subfolder is created only when something is actually written.
  # ============================================================================
  out_subdir <- file.path(out_dir, "gvr_summary")
  if (isTRUE(save_excel) || isTRUE(save_pdf)) {
    if (!dir.exists(out_subdir))
      dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)
  }

  # ============================================================================
  # Optional Excel export  ->  <out_dir>/gvr_summary/<file_prefix>.xlsx
  # Fixed filename (no timestamp): re-running overwrites the previous workbook.
  # ============================================================================
  if (isTRUE(save_excel)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      warning("gvr_summary: 'openxlsx' not installed; skipping Excel export.")
    } else {
      xlsx_name  <- sprintf("%s.xlsx", file_prefix)
      final_xlsx <- file.path(out_subdir, xlsx_name)
      if (file.exists(final_xlsx) && isTRUE(verbose))
        message(sprintf("  Overwriting existing Excel: %s", final_xlsx))

      sheet_map <- c(overview = "Overview", top_genes = "Top_genes",
                     variant_classification = "Variant_classification",
                     variant_type = "Variant_type", clin_sig = "Clinical", impact = "Impact")
      wb <- openxlsx::createWorkbook()
      hs <- openxlsx::createStyle(textDecoration = "bold", halign = "center")
      for (nm in names(sections)) {
        sh <- sheet_map[[nm]]
        openxlsx::addWorksheet(wb, sh)
        openxlsx::writeData(wb, sh, sections[[nm]], headerStyle = hs)
        openxlsx::freezePane(wb, sh, firstRow = TRUE)
        openxlsx::setColWidths(wb, sh, cols = seq_len(ncol(sections[[nm]])), widths = "auto")
      }
      # Write to a local temp file first, then shell-cp to out_subdir (FUSE-safe:
      # openxlsx uses zip random-access writes that can fail / 0-byte on S3-backed mounts).
      tmp_xlsx <- file.path(tempdir(), xlsx_name)
      wrote_ok <- tryCatch({ openxlsx::saveWorkbook(wb, tmp_xlsx, overwrite = TRUE); TRUE },
                           error = function(e) { warning(sprintf("gvr_summary: Excel write failed: %s", conditionMessage(e))); FALSE })
      if (wrote_ok) {
        cp <- system2("cp", c(shQuote(tmp_xlsx), shQuote(final_xlsx)))
        if (!file.exists(final_xlsx) || file.info(final_xlsx)$size == 0) {
          warning(sprintf("gvr_summary: copy to '%s' may have failed; Excel left at '%s'.",
                          final_xlsx, tmp_xlsx))
          final_xlsx <- tmp_xlsx
        }
        if (isTRUE(verbose)) message(sprintf("  Excel written: %s", final_xlsx))
      }
    }
  }

  # ============================================================================
  # Optional PDF report -> <out_dir>/gvr_summary/<file_prefix>_report.pdf (fixed name).
  # Layout: title/metadata page; each CHARTED section (top genes, variant
  # classification, predicted impact) shown as its table together with its bar
  # chart (side-by-side when the table is short, table full-width + chart below
  # when the table is tall); then the remaining chart-less tables (overview,
  # variant type, clinical significance) grouped/stacked on their own page(s).
  # Requires gridExtra + ggplot2; if unavailable the PDF is skipped with a warning
  # (sections are still returned).
  # ============================================================================
  if (isTRUE(save_pdf)) {
    have_pkgs <- requireNamespace("gridExtra", quietly = TRUE) &&
                 requireNamespace("ggplot2",   quietly = TRUE) &&
                 requireNamespace("grid",      quietly = TRUE) &&
                 requireNamespace("scales",    quietly = TRUE)
    if (!have_pkgs) {
      warning("gvr_summary: 'gridExtra'/'ggplot2'/'scales' not installed; skipping PDF report.")
    } else {
      # Internal renderer for the multi-page PDF (nested: single-use, not exported).
      # ASCII-only text (base pdf() substitutes non-ASCII glyphs); writes to a POSIX
      # temp path then shell-cp to final_pdf (R file.copy() can 0-byte on S3 FUSE).
      .gvr_summary_pdf <- function(sections, samples, meta, final_pdf, file_prefix = "gvr_summary") {
        PHYLO_BLUE <- "#0279EE"; PHYLO_GREEN <- "#75A025"; CREAM <- "#FAF9F3"; STONE <- "#ECE9E2"
        INK <- "#000000"; ORANGE <- "#E69F00"; VERMIL <- "#D55E00"

        W <- 8.27; H <- 11.69                       # A4 portrait inches
        USABLE_W_IN <- 7.3                          # printable body width (matches shipped)
        CEX_CORE_FLOOR <- 0.50                      # min table body font scale before col-paginate
        CEX_HEAD_FLOOR <- 0.55
        FACET_THRESHOLD <- 6L                       # > this many samples -> facet charts

        pal <- c("#0279EE", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#75A025",
                 "#D55E00", "#332288", "#AA4499", "#888888")
        fill_vals <- stats::setNames(pal[((seq_along(samples) - 1L) %% length(pal)) + 1L], samples)

        fmt <- function(x) format(x, big.mark = ",", trim = TRUE)

        # ---- measured grob size in inches (drive fit by measurement, not guesses) --------
        .grob_w_in <- function(g) grid::convertWidth(sum(g$widths),  "in", valueOnly = TRUE)
        .grob_h_in <- function(g) grid::convertHeight(sum(g$heights), "in", valueOnly = TRUE)

        # ---- dashboard table theme at a given core font scale ----------------------------
        .mk_theme <- function(cex_core = 0.72) {
          cex_head <- max(CEX_HEAD_FLOOR, cex_core + 0.06)
          gridExtra::ttheme_minimal(
            core    = list(fg_params = list(cex = cex_core, hjust = 1, x = 0.95),
                           bg_params = list(fill = c(CREAM, STONE), col = NA)),
            colhead = list(fg_params = list(col = "white", fontface = "bold", cex = cex_head),
                           bg_params = list(fill = PHYLO_BLUE, col = NA)))
        }

        # ---- robust table-grob factory (PORTED VERBATIM from shipped renderer): -----------
        #   row pagination -> font shrink -> column pagination. Returns flat list of grobs,
        #   each with attr "colspan" (e.g. "samples 1-6 of 20" or NULL) and "n_rows".
        .mk_table_grobs <- function(dt, rows_per_page = 34L) {
          df <- as.data.frame(dt, stringsAsFactors = FALSE)
          for (j in seq_len(ncol(df)))
            if (is.numeric(df[[j]])) df[[j]] <- format(df[[j]], big.mark = ",", trim = TRUE)
          # Dashboard tables are DENSE by default (start at 0.60, not gridExtra's larger
          # default) so the common few-sample report packs tightly; the ladder still steps
          # DOWN toward the floor when a wide many-sample table overflows the page width
          # (after which .mk_table_grobs column-paginates).
          cex_try <- c(0.60, 0.55, CEX_CORE_FLOOR)
          cex_use <- CEX_CORE_FLOOR
          for (cx in cex_try) {
            g0 <- gridExtra::tableGrob(df[1, , drop = FALSE], rows = NULL, theme = .mk_theme(cx))
            if (.grob_w_in(g0) <= USABLE_W_IN) { cex_use <- cx; break }
            cex_use <- cx
          }
          th <- .mk_theme(cex_use)
          ncol_df <- ncol(df)
          probe   <- gridExtra::tableGrob(df[1, , drop = FALSE], rows = NULL, theme = th)
          col_groups <- list(seq_len(ncol_df))
          if (ncol_df >= 3L && .grob_w_in(probe) > USABLE_W_IN) {
            cat_i <- 1L; total_i <- ncol_df
            mid   <- setdiff(seq_len(ncol_df), c(cat_i, total_i))
            wcol  <- vapply(seq_len(ncol_df), function(j)
              .grob_w_in(gridExtra::tableGrob(df[1, j, drop = FALSE], rows = NULL, theme = th)),
              numeric(1))
            fixed_w <- wcol[cat_i] + wcol[total_i]
            budget  <- USABLE_W_IN - fixed_w
            col_groups <- list(); cur <- integer(0); used <- 0
            for (j in mid) {
              if (length(cur) > 0L && used + wcol[j] > budget) {
                col_groups[[length(col_groups) + 1L]] <- c(cat_i, cur, total_i)
                cur <- integer(0); used <- 0
              }
              cur <- c(cur, j); used <- used + wcol[j]
            }
            if (length(cur)) col_groups[[length(col_groups) + 1L]] <- c(cat_i, cur, total_i)
          }
          n_mid_total <- if (length(col_groups) > 1L) ncol_df - 2L else NA_integer_
          npg_rows <- max(1L, ceiling(nrow(df) / rows_per_page))
          out <- list()
          for (cg_idx in seq_along(col_groups)) {
            cg <- col_groups[[cg_idx]]
            colspan_lbl <- if (length(col_groups) > 1L) {
              mids <- setdiff(cg, c(1L, ncol_df))
              sprintf("samples %d-%d of %d", mids[1] - 1L, mids[length(mids)] - 1L, n_mid_total)
            } else NULL
            for (pg in seq_len(npg_rows)) {
              rs <- ((pg - 1L) * rows_per_page + 1L):min(pg * rows_per_page, nrow(df))
              g  <- gridExtra::tableGrob(df[rs, cg, drop = FALSE], rows = NULL, theme = th)
              attr(g, "colspan") <- colspan_lbl
              attr(g, "n_rows")  <- length(rs)
              out[[length(out) + 1L]] <- g
            }
          }
          out
        }

        # ---- bar chart grob: grouped (<=6 samples) or faceted (>6) — PORTED -------------
        .bar_grob <- function(dt, cat_col, top = NULL, title = NULL) {
          d <- as.data.frame(dt, stringsAsFactors = FALSE)
          d <- d[d[[cat_col]] != "" & !is.na(d[[cat_col]]), , drop = FALSE]
          if (!is.null(top) && nrow(d) > top) d <- d[order(-d$Total)[seq_len(top)], , drop = FALSE]
          long <- do.call(rbind, lapply(samples, function(s)
            data.frame(Category = d[[cat_col]], Sample = s, n = d[[s]], stringsAsFactors = FALSE)))
          long$Category <- factor(long$Category, levels = d[[cat_col]][order(d$Total)])
          many <- length(samples) > FACET_THRESHOLD
          gg <- ggplot2::ggplot(long, ggplot2::aes(x = Category, y = n, fill = Sample)) +
            ggplot2::coord_flip() +
            ggplot2::scale_fill_manual(values = fill_vals) +
            ggplot2::labs(title = title, x = NULL, y = NULL, fill = NULL) +
            ggplot2::theme_minimal(base_size = 10) +
            ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", colour = PHYLO_BLUE, size = 13),
                           panel.grid.major.y = ggplot2::element_blank())
          if (many) {
            sn <- as.character(samples); pre <- ""
            if (length(sn) > 1L) {
              mn <- min(nchar(sn)); i <- 0L
              while (i < mn && length(unique(substr(sn, 1L, i + 1L))) == 1L) i <- i + 1L
              if (i >= 3L) pre <- substr(sn[1], 1L, i)
            }
            lab_fun <- function(x) if (nzchar(pre)) sub(paste0("^", pre), "", x) else x
            gg <- gg + ggplot2::geom_col(show.legend = FALSE) +
              ggplot2::facet_wrap(~ Sample, ncol = ceiling(sqrt(length(samples))),
                                  labeller = ggplot2::as_labeller(lab_fun)) +
              ggplot2::scale_y_continuous(breaks = scales::breaks_extended(3),
                                          labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
              ggplot2::theme(legend.position = "none",
                             strip.text = ggplot2::element_text(size = 7, face = "bold"),
                             axis.text  = ggplot2::element_text(size = 5.5),
                             panel.spacing = grid::unit(4, "pt"))
          } else {
            gg <- gg + ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.68) +
              ggplot2::scale_y_continuous(breaks = scales::breaks_extended(4),
                                          labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
              ggplot2::theme(legend.position = "top", legend.justification = "left",
                             legend.key.size = grid::unit(10, "pt"),
                             legend.text = ggplot2::element_text(size = 8),
                             axis.text = ggplot2::element_text(size = 9))
          }
          ggplot2::ggplotGrob(gg)
        }

        # ---- left-justify a table grob within a full-width row ---------------------------
        .left_just <- function(g) {
          w_in <- .grob_w_in(g)
          gridExtra::arrangeGrob(g, grid::nullGrob(), ncol = 2,
                                 widths = grid::unit.c(grid::unit(w_in, "in"), grid::unit(1, "null")))
        }

        # ---- KPI card. Big-number font auto-shrinks for long values (e.g. 7-digit cohort
        #      totals) so the number never kisses the card edge. ----------------------------
        .mk_kpi <- function(value, label, fill = PHYLO_BLUE, fg = "white") {
          nchar_v <- nchar(value)
          num_fs  <- if (nchar_v <= 6L) 30 else if (nchar_v <= 8L) 24 else 20
          grid::grobTree(
            grid::roundrectGrob(gp = grid::gpar(fill = fill, col = NA), r = grid::unit(6, "pt")),
            grid::textGrob(value, y = 0.60, gp = grid::gpar(col = fg, fontface = "bold", fontsize = num_fs)),
            grid::textGrob(label, y = 0.22, gp = grid::gpar(col = fg, fontsize = 11)))
        }

        # ---- cairo_pdf blank-first-page guard --------------------------------------------
        .make_new_page <- function() {
          first <- TRUE
          function() { if (first) first <<- FALSE else grid::grid.newpage() }
        }

        # ============================ HERO PAGE ===========================================
        .render_hero <- function(np) {
          np()
          title_grob <- grid::grobTree(
            grid::textGrob("germlinevaR \u2013 Cohort Summary", x = 0.02, y = 0.70, hjust = 0,
                           gp = grid::gpar(fontface = "bold", fontsize = 24, col = PHYLO_BLUE)),
            grid::textGrob(sprintf("%d sample%s  \u00b7  %s total variants  \u00b7  %s distinct genes",
                                   meta$n_samples, if (meta$n_samples == 1L) "" else "s",
                                   fmt(meta$n_total), fmt(meta$n_genes)),
                           x = 0.02, y = 0.26, hjust = 0, gp = grid::gpar(fontsize = 11, col = INK)))
          hi <- sections$impact$Total[sections$impact$IMPACT == "HIGH"]; if (!length(hi)) hi <- 0L
          kpis <- list(
            .mk_kpi(fmt(meta$n_total),    "Total variants",         PHYLO_BLUE),
            .mk_kpi(fmt(meta$n_samples),  "Samples",                PHYLO_GREEN),
            .mk_kpi(fmt(meta$n_genes),    "Distinct genes (known)", ORANGE, fg = INK),
            .mk_kpi(fmt(hi),              "HIGH-impact variants",   VERMIL))
          cards <- lapply(kpis, function(g)
            gridExtra::arrangeGrob(g, vp = grid::viewport(width = 0.94, height = 0.86)))
          cards_row <- gridExtra::arrangeGrob(grobs = cards, nrow = 1)
          chart_top <- .bar_grob(sections$top_genes, "Hugo_Symbol",
                                 title = "Top genes by variant count (top 10)", top = 10L)
          chart_cls <- .bar_grob(sections$variant_classification, "Variant_Classification",
                                 title = "Variant classification (top 10)", top = 10L)
          hero <- gridExtra::arrangeGrob(
            grobs = list(title_grob, cards_row, chart_top, chart_cls), ncol = 1,
            heights = grid::unit.c(grid::unit(0.95, "in"), grid::unit(1.5, "in"),
                                   grid::unit(1, "null"), grid::unit(1, "null")))
          grid::pushViewport(grid::viewport(width = grid::unit(W - 1, "in"),
                                            height = grid::unit(H - 1, "in")))
          grid::grid.draw(hero); grid::popViewport()
        }

        # ===================== AUTO SIDE-BY-SIDE REFERENCE SECTION ========================
        # Each section -> .mk_table_grobs (already font-shrunk / column-paginated). A
        # section's grob is SIDE-BY-SIDE eligible only when it produced a SINGLE grob whose
        # measured width fits half the content width (so two fit with a gap). Otherwise it
        # is FULL-WIDTH and stacks. This degrades safely for many-sample (wide / paginated)
        # tables, matching the user's "side-by-side if it fits, else stack" rule.
        .render_reference <- function(np, target_dev,
                                      title_in = 0.30, tgap_in = 0.10, bgap_in = 0.18,
                                      col_gap_in = 0.30, draw_frac = 0.97) {
          content_w <- W - 1
          budget    <- (H - 1) * draw_frac
          tbl_specs <- list(
            list(s = "overview",               t = "Overview"),
            list(s = "top_genes",              t = "Top genes"),
            list(s = "variant_classification", t = "Variant classification"),
            list(s = "variant_type",           t = "Variant type"),
            list(s = "clin_sig",               t = "Clinical significance"),
            list(s = "impact",                 t = "Functional impact (table)"))
          # MEASUREMENT must happen on a throwaway device that is opened AND closed HERE,
          # so the cairo report device is current during all subsequent drawing. (Leaving
          # pdf(NULL) open across the draw loop sends every page to the throwaway device
          # and the report comes out empty/0-byte.)
          items <- local({
            grDevices::pdf(NULL, width = W, height = H)
            on.exit(grDevices::dev.off())
            grid::pushViewport(grid::viewport(width = grid::unit(content_w, "in"),
                                              height = grid::unit(budget, "in")))
            it <- list()
            for (sp in tbl_specs) {
              gl <- .mk_table_grobs(sections[[sp$s]])
              multi <- length(gl) > 1L
              for (gi in seq_along(gl)) {
                g  <- gl[[gi]]
                cs <- attr(g, "colspan")
                lbl <- if (multi && !is.null(cs)) sprintf("%s (%s)", sp$t, cs) else sp$t
                wi <- .grob_w_in(g); hi <- .grob_h_in(g)
                # Side-by-side eligible iff this section produced a SINGLE grob (not
                # column-paginated). Whether two eligible tables actually fit a row is
                # decided later by the real pairwise width sum (a$w + gap + b$w <= content_w),
                # which is more permissive than a rigid half-width gate (a wide table can
                # still pair with a narrow one). Column-paginated (multi) tables stay
                # full-width and stack -- the many-sample safe fallback.
                eligible <- !multi
                it[[length(it) + 1L]] <- list(type = "table", lbl = lbl, g = g,
                                              w = wi, h = hi, pair_ok = eligible)
              }
            }
            grid::popViewport()
            it
          })
          # Closing the throwaway measurement device above can leave a DIFFERENT device
          # current (R falls back to the most-recent surviving device, not necessarily our
          # report device). Force the report device current before any drawing.
          if (target_dev %in% grDevices::dev.list()) grDevices::dev.set(target_dev)
          imp_chart <- list(type = "chart", lbl = "Functional impact (VEP IMPACT)",
                            g = .bar_grob(sections$impact, "IMPACT", title = NULL),
                            w = content_w, h = 2.4, pair_ok = FALSE)

          # ---- height-balanced pairing among pair_ok items; others stay single --------
          ord  <- order(vapply(items, function(u) u$h, numeric(1)), decreasing = TRUE)
          done <- rep(FALSE, length(items)); rows <- list()
          for (ii in seq_along(ord)) {
            i <- ord[ii]; if (done[i]) next
            a <- items[[i]]
            if (!isTRUE(a$pair_ok)) { done[i] <- TRUE
              rows[[length(rows)+1]] <- list(kind = "single", a = a, h = a$h); next }
            partner <- NA_integer_
            for (jj in seq_along(ord)) {
              j <- ord[jj]
              if (j == i || done[j] || !isTRUE(items[[j]]$pair_ok)) next
              if (a$w + col_gap_in + items[[j]]$w <= content_w) { partner <- j; break }
            }
            if (!is.na(partner)) {
              b <- items[[partner]]; done[i] <- TRUE; done[partner] <- TRUE
              rows[[length(rows)+1]] <- list(kind = "pair", a = a, b = b, h = max(a$h, b$h))
            } else {
              done[i] <- TRUE
              rows[[length(rows)+1]] <- list(kind = "single", a = a, h = a$h)
            }
          }
          rows[[length(rows)+1]] <- list(kind = "single", a = imp_chart, h = imp_chart$h)

          per_row_overhead <- title_in + tgap_in + bgap_in
          pages <- list(); cur <- list(); used <- 0
          for (r in rows) {
            cost <- r$h + per_row_overhead
            if (length(cur) > 0 && used + cost > budget) { pages[[length(pages)+1]] <- cur; cur <- list(); used <- 0 }
            cur[[length(cur)+1]] <- r; used <- used + cost
          }
          if (length(cur)) pages[[length(pages)+1]] <- cur

          block_for <- function(item, full_width, row_h) {
            body <- if (item$type == "table") { if (full_width) .left_just(item$g) else item$g } else item$g
            # Top-anchor the title within its band (vjust = 1 at y = 1) so its
            # baseline sits high and there is full clearance to the table body that
            # follows the tgap spacer -- prevents the body riding up under the title.
            ttl <- grid::textGrob(item$lbl, x = 0.02, y = 1, hjust = 0, vjust = 1,
                                  gp = grid::gpar(fontsize = 10.5, fontface = "bold", col = PHYLO_BLUE))
            pad <- max(row_h - item$h, 0)
            gridExtra::arrangeGrob(ttl, grid::nullGrob(), body, grid::nullGrob(), ncol = 1,
                                   heights = grid::unit.c(grid::unit(title_in, "in"), grid::unit(tgap_in, "in"),
                                                          grid::unit(item$h, "in"), grid::unit(pad, "in")))
          }

          for (pg in pages) {
            np()
            grobs <- list(grid::textGrob("Reference tables", x = 0.02, hjust = 0, vjust = 1, y = 1,
                                         gp = grid::gpar(fontface = "bold", fontsize = 10, col = PHYLO_GREEN)))
            rel <- c(0.18)
            for (r in pg) {
              if (r$kind == "pair") {
                ga <- block_for(r$a, FALSE, r$h); gb <- block_for(r$b, FALSE, r$h)
                row_grob <- gridExtra::arrangeGrob(ga, grid::nullGrob(), gb, ncol = 3,
                              widths = grid::unit.c(grid::unit(1, "null"), grid::unit(col_gap_in, "in"),
                                                    grid::unit(1, "null")))
              } else {
                row_grob <- block_for(r$a, TRUE, r$h)
              }
              grobs[[length(grobs)+1]] <- row_grob
              rel <- c(rel, r$h + title_in + tgap_in)
              grobs[[length(grobs)+1]] <- grid::nullGrob(); rel <- c(rel, bgap_in)
            }
            if (sum(rel) < budget) { grobs[[length(grobs)+1]] <- grid::nullGrob(); rel <- c(rel, budget - sum(rel)) }
            grid::pushViewport(grid::viewport(y = 0.5, height = grid::unit(budget, "in"),
                                              width = grid::unit(content_w, "in")))
            grid::grid.draw(gridExtra::arrangeGrob(grobs = grobs, ncol = 1, heights = grid::unit(rel, "in")))
            grid::popViewport()
          }
        }

        # ============================ DRIVE THE DEVICE ====================================
        tmp_pdf <- file.path(tempdir(), basename(final_pdf))
        grDevices::cairo_pdf(tmp_pdf, width = W, height = H, onefile = TRUE)
        cairo_dev <- grDevices::dev.cur()
        on.exit(if (cairo_dev %in% grDevices::dev.list()) grDevices::dev.off(cairo_dev), add = TRUE)
        np <- .make_new_page()
        .render_hero(np)
        .render_reference(np, target_dev = cairo_dev)
        if (cairo_dev %in% grDevices::dev.list()) grDevices::dev.off(cairo_dev)
        on.exit()
        if (!file.exists(tmp_pdf) || file.info(tmp_pdf)$size == 0)
          stop("PDF device produced no/zero-byte file.")
        system2("cp", c("-f", shQuote(tmp_pdf), shQuote(final_pdf)))
        invisible(final_pdf)
      }

      pdf_name  <- sprintf("%s_report.pdf", file_prefix)
      final_pdf <- file.path(out_subdir, pdf_name)
      if (file.exists(final_pdf) && isTRUE(verbose))
        message(sprintf("  Overwriting existing PDF report: %s", final_pdf))
      meta <- list(
        out_dir   = normalizePath(out_dir, mustWork = FALSE),
        n_samples = length(samples), samples = samples,
        n_total   = n_total,
        n_genes   = data.table::uniqueN(dt$Hugo_Symbol[is_known_gene]),
        n_nogene  = sum(!is_known_gene),
        generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
      ok_pdf <- tryCatch({
        .gvr_summary_pdf(sections, samples, meta, final_pdf, file_prefix)
        file.exists(final_pdf) && file.info(final_pdf)$size > 0
      }, error = function(e) {
        warning(sprintf("gvr_summary: PDF report failed: %s", conditionMessage(e))); FALSE })
      if (isTRUE(ok_pdf) && isTRUE(verbose)) message(sprintf("  PDF report written: %s", final_pdf))
    }
  }

  invisible(sections)
}

# NOTE: globalVariables() declarations for this package are consolidated in
# R/globals.R (one package-scoped block covering all functions).
