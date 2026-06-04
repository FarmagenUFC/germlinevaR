#' Multi-section summary of a germline MAF (read.gvr / gvr_filter output)
#'
#' @description
#' Produces a multi-section overview of a MAF-style table - the output of
#' [read.gvr()], or of [gvr_filter()] - covering variant burden, affected genes,
#' functional classes, clinical significance and predicted impact. Every section is
#' returned as a tidy `data.table` with one column per sample plus a `Total` column.
#' Optionally writes a multi-sheet Excel workbook. For a cohort oncoplot, see
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
#' This function is purely tabular: its only side effect is the optional Excel
#' workbook (written when `save_excel = TRUE`). It produces no plots; the cohort
#' oncoplot now lives in [gvr_oncoplot()].
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
#' @param save_excel Logical; if `TRUE`, write a multi-sheet `.xlsx` to `out_dir`.
#'   Default `FALSE`.
#' @param out_dir Output directory for the Excel workbook. Created if it does not
#'   exist. Default `"."` (current working directory).
#' @param file_prefix Filename prefix for the written workbook. Default
#'   `"gvr_summary"`; a `YYYYMMDD_HHMMSS` timestamp is appended.
#' @param verbose Logical; if `TRUE` (default) print a compact console digest and the
#'   path of the Excel file when written.
#'
#' @return Invisibly, a named list of `data.table`s: `overview`, `top_genes`,
#'   `variant_classification`, `variant_type`, `clin_sig`, `impact`. The return value
#'   is identical whether or not the Excel file is written.
#'
#' @section Dependencies:
#' Core summary uses \pkg{data.table}; the optional Excel export uses \pkg{openxlsx}.
#' If \pkg{openxlsx} is unavailable, the Excel export is skipped with a warning and the
#' sections are still returned.
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
#' ## compute only (no files written); print a console digest
#' s <- gvr_summary(maf)
#' s$variant_classification          # inspect a section
#' s$impact                          # HIGH -> MODIFIER, severity order
#'
#' ## summarise filtered hits and write the Excel workbook
#' gvr_summary(gvr_filter(maf), save_excel = TRUE, out_dir = "results/summary")
#' }
#'
#' @importFrom data.table as.data.table data.table setnames setcolorder setorder
#'   uniqueN copy :=
#' @importFrom utils head
#' @importFrom openxlsx createWorkbook
#' @export
gvr_summary <- function(maf,
                        sample_col     = "Tumor_Sample_Barcode",
                        top_n_genes    = 20,
                        save_excel     = FALSE,
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
    m <- as.data.table(unclass(tab), keep.rownames = TRUE)
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
  # Optional Excel export
  # ============================================================================
  if (isTRUE(save_excel)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      warning("gvr_summary: 'openxlsx' not installed; skipping Excel export.")
    } else {
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
      xlsx_name <- sprintf("%s_%s.xlsx", file_prefix, ts)
      final_xlsx <- file.path(out_dir, xlsx_name)

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
      # Write to a local temp file first, then shell-cp to out_dir (FUSE-safe: openxlsx
      # uses zip random-access writes that can fail / 0-byte on S3-backed mounts).
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

  invisible(sections)
}

# Silence R CMD check NOTEs for data.table non-standard-evaluation symbols
# (column references and special symbols used inside `dt[...]`). These are not
# undefined globals; this is the idiomatic data.table remedy.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".", ".N", ".SD", ".__sample__",
    "Hugo_Symbol", "Total", "Variant_Classification", "CLIN_SIG"
  ))
}
