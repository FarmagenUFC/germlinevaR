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
                 requireNamespace("grid",      quietly = TRUE)
    if (!have_pkgs) {
      warning("gvr_summary: 'gridExtra'/'ggplot2' not installed; skipping PDF report.")
    } else {
      # Internal renderer for the multi-page PDF (nested: single-use, not exported).
      # ASCII-only text (base pdf() substitutes non-ASCII glyphs); writes to a POSIX
      # temp path then shell-cp to final_pdf (R file.copy() can 0-byte on S3 FUSE).
      .gvr_summary_pdf <- function(sections, samples, meta, final_pdf, file_prefix = "gvr_summary") {
        PHYLO_BLUE <- "#0279EE"; PHYLO_GREEN <- "#75A025"; CREAM <- "#FAF9F3"; STONE <- "#ECE9E2"
        # bar-chart fill palette: Okabe-Ito-style, colorblind-safe ordering. Blue + orange
        # lead (the safest 2-class pair) so the common 2-sample case is maximally legible.
        pal <- c("#0279EE", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#75A025",
                 "#D55E00", "#332288", "#AA4499", "#888888")
        fill_vals <- stats::setNames(pal[((seq_along(samples) - 1L) %% length(pal)) + 1L], samples)

        # ---- Layout constants (A4 portrait), tunable in one place ---------------
        USABLE_W_IN    <- 7.3    # printable body width  (A4 8.27in - margins)
        USABLE_H_IN    <- 9.0    # printable body height (A4 11.69in - header/footer)
        CEX_CORE_FLOOR <- 0.50   # smallest table body font scale before paginating cols
        CEX_HEAD_FLOOR <- 0.55   # smallest table header font scale
        FACET_THRESHOLD <- 6L    # > this many samples -> facet charts instead of dodge

        # table theme at a given font scale (core cex). header scales in proportion.
        .mk_theme <- function(cex_core = 0.72) {
          cex_head <- max(CEX_HEAD_FLOOR, cex_core + 0.06)
          gridExtra::ttheme_minimal(
            core    = list(fg_params = list(cex = cex_core, hjust = 1, x = 0.95),
                           bg_params = list(fill = c(CREAM, STONE), col = NA)),
            colhead = list(fg_params = list(col = "white", fontface = "bold", cex = cex_head),
                           bg_params = list(fill = PHYLO_BLUE, col = NA)))
        }
        .tt <- .mk_theme(0.72)

        # measured grob size in inches (used to drive fit decisions, not guesses)
        .grob_w_in <- function(g) grid::convertWidth(sum(g$widths),  "in", valueOnly = TRUE)
        .grob_h_in <- function(g) grid::convertHeight(sum(g$heights), "in", valueOnly = TRUE)
        # Draw a (short) table grob anchored to the TOP of the body viewport instead of
        # vertically centred, so short tables don't leave a large empty band above them.
        # The page body is ~0.86 of an 11.69in A4 page (~10.05in); we cap the grob's
        # natural height at the body and pin it to the top edge.
        .body_in <- 0.86 * 11.69
        .draw_table_top <- function(g) {
          frac <- min(1, .grob_h_in(g) / .body_in)
          vp <- grid::viewport(y = 1, just = "top", height = frac)
          grid::pushViewport(vp); grid::grid.draw(g); grid::popViewport()
        }

        # Build one or more tableGrobs for a section. Robust to MANY samples:
        #   1) ROW pagination (rows_per_page) as before;
        #   2) FONT auto-scale: shrink core cex toward CEX_CORE_FLOOR while the grob is
        #      wider than USABLE_W_IN;
        #   3) COLUMN pagination: if still too wide at the font floor, split the
        #      per-sample columns into groups, repeating the FIRST column (category) and
        #      the LAST column ("Total") on every column-page; pages are labelled via the
        #      "colspan" attribute (e.g. "samples 1-6 of 20").
        # Returns a flat list of tableGrobs; each carries attr(,"colspan") (or NULL).
        .mk_table_grobs <- function(dt, rows_per_page = 28L) {
          df <- as.data.frame(dt, stringsAsFactors = FALSE)
          for (j in seq_len(ncol(df)))
            if (is.numeric(df[[j]])) df[[j]] <- format(df[[j]], big.mark = ",", trim = TRUE)

          # choose a font scale that fits width (using the full df width as the probe);
          # step down from 0.72 to the floor.
          cex_try <- c(0.72, 0.66, 0.60, 0.55, CEX_CORE_FLOOR)
          cex_use <- CEX_CORE_FLOOR
          for (cx in cex_try) {
            g0 <- gridExtra::tableGrob(df[1, , drop = FALSE], rows = NULL, theme = .mk_theme(cx))
            if (.grob_w_in(g0) <= USABLE_W_IN) { cex_use <- cx; break }
            cex_use <- cx
          }
          th <- .mk_theme(cex_use)

          # decide column groups. ncol layout: [1]=category ... [ncol]=Total.
          ncol_df <- ncol(df)
          probe   <- gridExtra::tableGrob(df[1, , drop = FALSE], rows = NULL, theme = th)
          col_groups <- list(seq_len(ncol_df))   # default: all columns in one group
          if (ncol_df >= 3L && .grob_w_in(probe) > USABLE_W_IN) {
            cat_i   <- 1L
            total_i <- ncol_df
            mid     <- setdiff(seq_len(ncol_df), c(cat_i, total_i))   # sample columns
            # per-column width estimate (inches) from a 1-row probe; pack greedily.
            wcol <- vapply(seq_len(ncol_df), function(j)
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
              attr(g, "n_rows")  <- length(rs)   # data rows in this grob (for packing)
              out[[length(out) + 1L]] <- g
            }
          }
          out
        }
        .section_header <- function(title, subtitle = NULL) {
          grid::grid.text(title, x = 0.06, y = 0.95, just = c("left", "top"),
                          gp = grid::gpar(fontsize = 15, fontface = "bold", col = PHYLO_BLUE))
          if (!is.null(subtitle))
            grid::grid.text(subtitle, x = 0.06, y = 0.905, just = c("left", "top"),
                            gp = grid::gpar(fontsize = 9, col = "grey40"))
        }
        # Build the per-sample dodged bar chart as a GROB (so it can be composed on
        # the same page as a table via gridExtra::arrangeGrob, rather than printed to
        # its own page). The ggplot spec is unchanged from the previous renderer.
        .bar_grob <- function(title, dt, cat_col, top = NULL) {
          d <- data.table::as.data.table(data.table::copy(dt))
          if (!is.null(top) && nrow(d) > top) d <- utils::head(d[order(-d[["Total"]])], top)
          long <- data.table::melt(d, id.vars = cat_col, measure.vars = samples,
                                   variable.name = "Sample", value.name = "n")
          data.table::setnames(long, cat_col, "Category")
          long[, Category := factor(Category, levels = rev(d[[cat_col]]))]
          many <- length(samples) > FACET_THRESHOLD
          gg <- ggplot2::ggplot(long, ggplot2::aes(x = Category, y = n, fill = Sample)) +
            ggplot2::coord_flip() +
            ggplot2::scale_fill_manual(values = fill_vals) +
            ggplot2::labs(title = NULL, x = NULL, y = "Variants", fill = NULL) +
            ggplot2::theme_minimal(base_size = 9) +
            ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", colour = PHYLO_BLUE))
          if (many) {
            # MANY samples: one small panel per sample (no dodge, no legend). To keep the
            # grid legible we (a) abbreviate strip labels by dropping the longest common
            # prefix shared by all sample names, (b) use only a few count-axis breaks
            # with SI-style labels so ticks don't collide in narrow panels, and (c)
            # shrink category + axis text. Panels wrap into ~sqrt(n) columns.
            sn  <- as.character(samples)
            pre <- ""
            if (length(sn) > 1L) {
              # longest common prefix across sample names
              mn <- min(nchar(sn)); i <- 0L
              while (i < mn && length(unique(substr(sn, 1L, i + 1L))) == 1L) i <- i + 1L
              if (i >= 3L) pre <- substr(sn[1], 1L, i)   # only strip a meaningful prefix
            }
            lab_fun <- function(x) if (nzchar(pre)) sub(paste0("^", pre), "", x) else x
            gg <- gg + ggplot2::geom_col(show.legend = FALSE) +
              ggplot2::facet_wrap(~ Sample, ncol = ceiling(sqrt(length(samples))),
                                  labeller = ggplot2::as_labeller(lab_fun)) +
              ggplot2::scale_y_continuous(
                breaks = scales::breaks_extended(3),
                labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
              ggplot2::theme(legend.position = "none",
                             strip.text = ggplot2::element_text(size = 7, face = "bold"),
                             axis.text.y = ggplot2::element_text(size = 5.5),
                             axis.text.x = ggplot2::element_text(size = 5.5),
                             panel.spacing = grid::unit(4, "pt"))
          } else {
            # FEW samples: grouped (dodged) bars with a top legend (the legible default).
            gg <- gg + ggplot2::geom_col(position = "dodge") +
              ggplot2::theme(legend.position = "top")
          }
          ggplot2::ggplotGrob(gg)
        }

        # Render ONE charted section: its table together with its bar chart.
        #  - SHORT table (<= tall_threshold rows): table LEFT, chart RIGHT (side-by-side).
        #  - TALL  table (>  tall_threshold rows): table FULL-WIDTH on top, chart BELOW.
        # Tables wider than one page paginate via .mk_table_grobs (28 rows/page); the
        # chart is drawn on the LAST page of the section (beside or below the final
        # table chunk per the short/tall rule applied to that chunk's row count).
        .charted_block <- function(title, dt, cat_col, top = NULL,
                                   chart_title = NULL, tall_threshold = 10L) {
          tgrobs <- .mk_table_grobs(dt)
          cgrob  <- .bar_grob(chart_title, dt, cat_col, top = top)
          npg    <- length(tgrobs)
          facet  <- length(samples) > FACET_THRESHOLD

          if (facet) {
            # MANY samples: the faceted chart is a full multi-panel grid and the table is
            # wide (column-paginated). Cramming both on a page clips. So render the
            # table page(s) TABLE-ONLY, then the chart on its OWN full page.
            for (k in seq_len(npg)) {
              grid::grid.newpage()
              ttl <- if (npg > 1L) sprintf("%s (table %d/%d)", title, k, npg) else title
              .section_header(ttl, subtitle = attr(tgrobs[[k]], "colspan"))
              body <- grid::viewport(y = 0.46, height = 0.86)
              grid::pushViewport(body); .draw_table_top(tgrobs[[k]]); grid::popViewport()
            }
            grid::grid.newpage()
            .section_header(title, subtitle = "per-sample chart")
            cbody <- grid::viewport(y = 0.46, height = 0.86)
            grid::pushViewport(cbody); grid::grid.draw(cgrob); grid::popViewport()
            return(invisible(NULL))
          }

          for (k in seq_len(npg)) {
            grid::grid.newpage()
            ttl <- if (npg > 1L) sprintf("%s (page %d/%d)", title, k, npg) else title
            cs  <- attr(tgrobs[[k]], "colspan")
            .section_header(ttl, subtitle = cs)
            body <- grid::viewport(y = 0.46, height = 0.86)
            grid::pushViewport(body)
            if (k < npg) {
              # intermediate table page (only happens for very long tables): table only
              .draw_table_top(tgrobs[[k]])
            } else {
              tall <- nrow(dt) > tall_threshold
              if (tall) {
                # full-width table on top, chart below (uses the whole body height)
                g <- gridExtra::arrangeGrob(tgrobs[[k]], cgrob, ncol = 1,
                                            heights = grid::unit(c(0.55, 0.45), "npc"))
                grid::grid.draw(g)
              } else {
                # side-by-side: table left, chart right. Draw it into a TOP-ANCHORED
                # band whose height scales with the row count, so a short chart (e.g.
                # the 4-category IMPACT chart) is not stretched over the full page
                # height. band_h is a fraction of the body viewport.
                g <- gridExtra::arrangeGrob(tgrobs[[k]], cgrob, ncol = 2,
                                            widths = grid::unit(c(0.46, 0.54), "npc"))
                # band_h is a fraction of the body viewport. Keep it small so a short
                # chart (e.g. 4-category IMPACT) is compact, not stretched: ~0.30 of
                # the body for a 4-row table, growing modestly toward the tall_threshold.
                band_h <- max(0.22, min(0.42, 0.06 + 0.06 * nrow(dt)))
                band   <- grid::viewport(y = 1 - band_h / 2, height = band_h)
                grid::pushViewport(band)
                grid::grid.draw(g)
                grid::popViewport()
              }
            }
            grid::popViewport()
          }
        }

        # Render the chart-less tables GROUPED: stack several section tables on one
        # page (each with a small sub-header), top to bottom, packing as many as fit
        # and spilling to a new page at a section boundary. `items` is a list of
        # lists(title=, dt=, subtitle=). A single table taller than one page falls back
        # to .mk_table_grobs pagination on its own page(s).
        .draw_grouped_tables <- function(items, page_title = NULL,
                                        rows_per_page = 30L) {
          # Pre-expand every item into its rendered table GROB(s). A wide table is
          # column-paginated by .mk_table_grobs into several grobs that stack
          # VERTICALLY (one per sample-column group), each a self-contained table with
          # its own header. We therefore pack at the GROB level: each grob is an atomic
          # unit whose vertical cost is its own row count (+2 for sub-title+header).
          # Packing at the item level (the old behaviour) undercounted wide tables and
          # overflowed the page (e.g. a 19-row clin_sig split into 3 groups = 3x the
          # height but counted once), clipping rows and the repeated sub-title.
          units <- list()
          for (it in items) {
            gl <- .mk_table_grobs(it$dt)
            multi <- length(gl) > 1L
            for (gi in seq_along(gl)) {
              cs  <- attr(gl[[gi]], "colspan")
              lbl <- if (multi && !is.null(cs)) sprintf("%s (%s)", it$title, cs) else it$title
              nr  <- attr(gl[[gi]], "n_rows"); if (is.null(nr)) nr <- nrow(it$dt)
              units[[length(units) + 1L]] <- list(lbl = lbl, g = gl[[gi]],
                                                  cost = nr + 2L, rows = nr)
            }
          }
          # greedy bin-pack the grob units by vertical cost
          bins <- list(); cur <- list(); used <- 0L
          for (u in units) {
            if (used > 0L && used + u$cost > rows_per_page) {
              bins[[length(bins) + 1L]] <- cur; cur <- list(); used <- 0L
            }
            cur[[length(cur) + 1L]] <- u; used <- used + u$cost
          }
          if (length(cur)) bins[[length(bins) + 1L]] <- cur
          np <- length(bins)
          for (bi in seq_len(np)) {
            grid::grid.newpage()
            hdr <- if (!is.null(page_title)) {
              if (np > 1L) sprintf("%s (page %d/%d)", page_title, bi, np) else page_title
            } else NULL
            if (!is.null(hdr)) .section_header(hdr)
            grp <- bins[[bi]]
            sub_grobs <- list(); rel_h <- numeric(0)
            for (u in grp) {
              sub_grobs[[length(sub_grobs) + 1L]] <-
                grid::textGrob(u$lbl, x = 0.02, just = "left",
                               gp = grid::gpar(fontsize = 11, fontface = "bold",
                                               col = PHYLO_BLUE))
              rel_h <- c(rel_h, 0.6)
              sub_grobs[[length(sub_grobs) + 1L]] <- u$g
              rel_h <- c(rel_h, max(1.2, u$rows * 0.5))
            }
            # Append a bottom spacer so a lightly-filled page (e.g. just Overview +
            # Variant type) packs toward the TOP instead of spreading across the whole
            # page. Spacer weight = remaining budget up to the page row capacity.
            content_w <- sum(rel_h)
            cap_w     <- rows_per_page * 0.5            # same 0.5/row scale as tables
            if (content_w < cap_w) {
              sub_grobs[[length(sub_grobs) + 1L]] <- grid::nullGrob()
              rel_h <- c(rel_h, cap_w - content_w)
            }
            stacked <- gridExtra::arrangeGrob(grobs = sub_grobs, ncol = 1,
                                              heights = grid::unit(rel_h, "null"))
            vp <- grid::viewport(y = 0.46, height = 0.86)
            grid::pushViewport(vp); grid::grid.draw(stacked); grid::popViewport()
          }
        }

        # Render section 3 (Predicted impact, charted) TOGETHER with sections 4
        # (Overview) and 5 (Variant type) on ONE page when they fit; otherwise split.
        # FIT is decided by MEASURED grob heights vs USABLE_H_IN (not row counts), so it
        # auto-splits if any of the three grows (more rows or many-sample facet chart).
        .impact_combined_block <- function(impact, overview, vartype) {
          imp_tg <- .mk_table_grobs(impact$dt)         # impact table grob(s)
          imp_cg <- .bar_grob(impact$chart_title, impact$dt, impact$cat_col)
          ov_tg  <- .mk_table_grobs(overview$dt)
          vt_tg  <- .mk_table_grobs(vartype$dt)

          # Conditions that force the split (each piece then laid out the standard way):
          #  - impact table needed column pagination (wide -> many samples),
          #  - the chart is in facet (many-sample) mode (tall, needs a full page),
          #  - any companion table is itself column-paginated,
          #  - measured stacked height exceeds the usable page height.
          facet_mode <- length(samples) > FACET_THRESHOLD
          multi_imp  <- length(imp_tg) > 1L
          multi_comp <- length(ov_tg) > 1L || length(vt_tg) > 1L

          # impact band height: table vs a compact chart height, whichever is taller.
          imp_band_in <- max(.grob_h_in(imp_tg[[1]]), 2.6)   # compact chart ~2.6in
          # companion stack: sub-title (~0.28in) + table height, per item.
          comp_in <- (0.28 + .grob_h_in(ov_tg[[1]])) + (0.28 + .grob_h_in(vt_tg[[1]]))
          total_in <- imp_band_in + 0.4 + comp_in           # +0.4 spacing

          if (facet_mode || multi_imp || multi_comp || total_in > USABLE_H_IN) {
            # ---- SPLIT: impact on its own page, then 4+5 grouped (existing paths) ----
            .charted_block("3. Predicted impact (severity order)", impact$dt,
                           cat_col = impact$cat_col, chart_title = impact$chart_title)
            .draw_grouped_tables(
              list(list(title = "4. Overview",      dt = overview$dt),
                   list(title = "5. Variant type",  dt = vartype$dt)),
              page_title = "Summary tables")
            return(invisible(FALSE))
          }

          # ---- ONE PAGE: impact table+chart band on top, 4 & 5 stacked below -------
          grid::grid.newpage()
          .section_header("3. Predicted impact (severity order)")
          body <- grid::viewport(y = 0.46, height = 0.86)
          grid::pushViewport(body)

          # proportion of body for the impact band vs the companion tables, from inches.
          band_frac <- max(0.30, min(0.55, imp_band_in / (imp_band_in + comp_in + 0.4)))
          # impact band (top): table left, compact chart right. A small right gutter
          # column keeps the chart's largest axis label clear of the page margin.
          imp_g <- gridExtra::arrangeGrob(imp_tg[[1]], imp_cg, grid::nullGrob(), ncol = 3,
                                          widths = grid::unit(c(0.45, 0.52, 0.03), "npc"))
          band  <- grid::viewport(y = 1 - band_frac / 2, height = band_frac)
          grid::pushViewport(band); grid::grid.draw(imp_g); grid::popViewport()

          # companion tables (below the band): Overview then Variant type, sub-titled
          sub_grobs <- list(); rel_h <- numeric(0)
          for (it in list(list(title = "4. Overview", g = ov_tg[[1]], n = nrow(overview$dt)),
                          list(title = "5. Variant type", g = vt_tg[[1]], n = nrow(vartype$dt)))) {
            sub_grobs[[length(sub_grobs) + 1L]] <-
              grid::textGrob(it$title, x = 0.02, just = "left",
                             gp = grid::gpar(fontsize = 11, fontface = "bold", col = PHYLO_BLUE))
            rel_h <- c(rel_h, 0.5)
            sub_grobs[[length(sub_grobs) + 1L]] <- it$g
            rel_h <- c(rel_h, max(1.0, it$n * 0.5))
          }
          sub_grobs[[length(sub_grobs) + 1L]] <- grid::nullGrob()   # bottom spacer
          rel_h <- c(rel_h, sum(rel_h) * 0.15)
          stacked <- gridExtra::arrangeGrob(grobs = sub_grobs, ncol = 1,
                                            heights = grid::unit(rel_h, "null"))
          lower <- grid::viewport(y = (1 - band_frac) / 2, height = 1 - band_frac - 0.02)
          grid::pushViewport(lower); grid::grid.draw(stacked); grid::popViewport()

          grid::popViewport()   # body
          invisible(TRUE)
        }

        tmp_pdf <- file.path(tempdir(), basename(final_pdf))
        grDevices::pdf(tmp_pdf, width = 8.27, height = 11.69, onefile = TRUE)  # A4 portrait
        on.exit(grDevices::dev.off(), add = TRUE)

        # Page 1: title / metadata (ASCII only)
        grid::grid.newpage()
        grid::grid.text("germlinevaR - gvr_summary report", y = 0.78,
                        gp = grid::gpar(fontsize = 22, fontface = "bold", col = PHYLO_BLUE))
        meta_lines <- c(
          sprintf("Output folder: %s", file.path(meta$out_dir, "gvr_summary")),
          sprintf("Generated: %s", meta$generated),
          "",
          sprintf("Samples (%d): %s", meta$n_samples, paste(meta$samples, collapse = ", ")),
          sprintf("Total variants: %s", format(meta$n_total, big.mark = ",")),
          sprintf("Distinct genes (known): %s", format(meta$n_genes, big.mark = ",")),
          sprintf("Variants with no gene symbol: %s", format(meta$n_nogene, big.mark = ",")))
        grid::grid.text(paste(meta_lines, collapse = "\n"), x = 0.5, y = 0.5,
                        gp = grid::gpar(fontsize = 12), just = "center")
        grid::grid.text("germlinevaR :: VEP germline VCF -> MAF toolkit", y = 0.07,
                        gp = grid::gpar(fontsize = 9, col = "grey50"))

        # --- Charted sections: each table shown together with its bar chart ---
        # (tall tables -> full-width with chart below; short tables -> side-by-side).
        .charted_block("1. Top genes (by total variants)", sections$top_genes,
                       cat_col = "Hugo_Symbol", top = 15L,
                       chart_title = "Top 15 genes (per sample)")
        .charted_block("2. Variant classification", sections$variant_classification,
                       cat_col = "Variant_Classification",
                       chart_title = "Variant classification (per sample)")

        # Sections 3 (Predicted impact, charted) + 4 (Overview) + 5 (Variant type)
        # share ONE page when they fit (measured), else auto-split (see fn).
        .impact_combined_block(
          impact   = list(dt = sections$impact, cat_col = "IMPACT",
                          chart_title = "Predicted impact (per sample)"),
          overview = list(dt = sections$overview),
          vartype  = list(dt = sections$variant_type))

        # Section 6 (Clinical significance) on its own grouped page (19 tokens).
        .draw_grouped_tables(
          list(list(title = "6. Clinical significance (CLIN_SIG tokens)",
                    dt = sections$clin_sig)),
          page_title = "Summary tables")

        grDevices::dev.off(); on.exit()   # close device now so the file is complete before copy
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
