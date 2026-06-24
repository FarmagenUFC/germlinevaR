#' Modular, individually-toggleable filtering of a read.gvr table
#'
#' @description
#' Applies a set of independent variant filters to the MAF-like table produced by
#' [read.gvr()]. Each distinct filter is its own argument; setting an argument to
#' `NULL` disables that filter entirely (no rows removed by it). With all defaults,
#' `gvr_filter(gvr)` reproduces the canonical rare / clinically-relevant /
#' called-genotype pipeline (AF filters + CLIN_SIG + GT exclusion).
#'
#' @details
#' Filters are applied in a fixed order; each step operates on the survivors of the
#' previous one:
#' \enumerate{
#'   \item gnomAD exome AF        - `gnomADe_AF`        (+ `gnomADe_AF_keep_missing`)
#'   \item gnomAD genome / VCF AF - `AF`                (+ `AF_keep_missing`)
#'   \item ABraOM AF             - `ABraOM_AF`         (+ `ABraOM_AF_keep_missing`)
#'   \item Clinical significance  - `clin_sig_terms`    (+ `clin_sig_keep_missing`)
#'   \item Remove benign          - `remove_benign`
#'   \item Biotype                - `biotype_keep`
#'   \item Genotype exclusion     - `gt_exclude`
#'   \item Variant classification - `vc_nonSyn`
#'   \item Gene subset            - `genes`
#' }
#'
#' Important data notes (true of [read.gvr()] output):
#' \itemize{
#'   \item AF columns are CHARACTER (e.g. `"0.8781"`), so they are coerced with
#'     `as.numeric()` before comparison.
#'   \item "Missing" means EITHER `NA` OR empty string `""` ([read.gvr()] uses `""`
#'     for absent values). Both are treated as missing everywhere in this function.
#'   \item `CLIN_SIG` matching is SUBSTRING + case-insensitive. A compound annotation
#'     such as `"pathogenic&benign"` or
#'     `"uncertain_significance&likely_benign&benign"` is KEPT because it CONTAINS a
#'     wanted term. This matches the dplyr `str_detect()` / `grepl()` convention. Use
#'     exact-token matching only if you split `CLIN_SIG` yourself.
#'   \item The default `gt_exclude = c("0", "0/0")` is a no-op on data whose `GT`
#'     column only contains called alt genotypes (e.g. `0/1`, `1/1`, `1/2`); it is
#'     retained for portability to data that does carry `"0"`/`"0/0"`.
#' }
#'
#' The input is never modified: filtering operates on an internal `data.table` copy.
#'
#' @param gvr data.table / data.frame from [read.gvr()] (or compatible). Filtered in a
#'   copy; the input object is not modified.
#'
#' @param gnomADe_AF Numeric upper threshold for the gnomAD exome AF column
#'   `gnomADe_AF` (keep rows with AF < threshold). `NULL` disables this filter.
#'   Default 0.01.
#' @param AF Numeric upper threshold for the `AF` column (gnomAD genome / VCF allele
#'   frequency). `NULL` disables this filter. Default 0.01.
#' @param ABraOM_AF Numeric upper threshold for the Brazilian-cohort (ABraOM SABE 609)
#'   column `ABraOM_AF`. `NULL` disables this filter. Default 0.01.
#' @param gnomADe_AF_keep_missing Logical; if TRUE (default), keep rows whose
#'   `gnomADe_AF` is missing (NA or ""); if FALSE drop them. Ignored when
#'   `gnomADe_AF` is NULL.
#' @param AF_keep_missing Logical; missing-value handling for the `AF` filter.
#'   TRUE (default) keeps missing. Ignored when `AF` is NULL.
#' @param ABraOM_AF_keep_missing Logical; missing-value handling for the ABraOM
#'   filter. TRUE (default) retains variants absent from the Brazilian cohort
#'   (where absence often means "not catalogued", not "common"). Ignored when
#'   `ABraOM_AF` is NULL.
#'
#' @param clin_sig_terms Character vector of clinical-significance terms to keep
#'   (substring, case-insensitive, OR-combined). `NULL` disables the CLIN_SIG filter.
#'   Default: c("likely_pathogenic","pathogenic","uncertain_significance").
#' @param clin_sig_keep_missing Logical; if TRUE (default) rows with missing CLIN_SIG
#'   (NA/"") are kept. Only relevant when `clin_sig_terms` is non-NULL.
#' @param remove_benign Logical; if TRUE, remove rows whose `CLIN_SIG` contains
#'   "benign" (substring, case-insensitive). This catches `benign`,
#'   `likely_benign`, and compound annotations like
#'   `"uncertain_significance&likely_benign"`. Applied AFTER the `clin_sig_terms`
#'   keep-filter, so a row that matched a wanted term but also contains "benign"
#'   is still removed. `FALSE` (default) does not remove benign rows.
#' @param biotype_keep Character vector of BIOTYPE values to keep (exact match via %in%).
#'   `NULL` (default) disables the biotype filter — all biotypes are kept.
#'   Pass e.g. `c("protein_coding", "protein_coding_LoF")` to restrict to
#'   protein-coding transcripts.
#' @param gt_exclude Character vector of GT values to remove (exact match). `NULL`
#'   disables the genotype filter. Default: c("0","0/0").
#' @param vc_nonSyn Logical or character vector. Controls which
#'   `Variant_Classification` values are retained. `FALSE` (default) keeps all.
#'   `TRUE` keeps only the 9 protein-altering classes (Frame_Shift_Del,
#'   Frame_Shift_Ins, Splice_Site, Translation_Start_Site, Nonsense_Mutation,
#'   Nonstop_Mutation, In_Frame_Del, In_Frame_Ins, Missense_Mutation).
#'   A custom character vector keeps only those classifications. Rows with
#'   missing/blank `Variant_Classification` are removed when this filter is active.
#' @param missense_only Logical; if `TRUE`, keep only rows whose
#'   `Variant_Classification` equals `"Missense_Mutation"` (added in vN+5).
#'   Default `FALSE` preserves prior behaviour byte-for-byte. Combines
#'   non-contradictorily with `vc_nonSyn`: `vc_nonSyn` runs first (keeping
#'   9 protein-altering classes), then `missense_only` narrows to the missense
#'   subset. Errors with a clear message if `Variant_Classification` is missing.
#' @param genes Character vector of `Hugo_Symbol`s to keep (exact, case-insensitive),
#'   or `NULL` (default) to keep all genes.
#' @param save_excel Logical; if TRUE, also write the FILTERED table to an `.xlsx`
#'   workbook (single `"Filtered"` sheet) at `<out_dir>/<file_prefix>.xlsx`. Requires
#'   the \pkg{openxlsx} package (a `Suggests` dependency); if it is not installed the
#'   export is skipped with a warning. Default FALSE. The write is a side effect only:
#'   the returned `data.table` is identical whether or not `save_excel` is TRUE.
#' @param out_dir Output directory for the Excel file. `NULL` (default) uses the
#'   current working directory. Created if it does not exist. Only used when
#'   `save_excel = TRUE`.
#' @param file_prefix Filename prefix (without extension) for the Excel file. Default
#'   `"gvr_filter"` -> `gvr_filter.xlsx`. Only used when `save_excel = TRUE`.
#' @param verbose Logical; if TRUE (default) print a per-filter breakdown (rows in -> out
#'   and rows removed by each active step).
#'
#' @return A `data.table` of the surviving rows, with the same columns as the input.
#'   A plain `data.frame` input is returned as a `data.table`. The input object is not
#'   modified. With `verbose = TRUE`, a per-filter breakdown (rows in -> out, and rows
#'   removed by each active step) is printed as it runs.
#'
#' @seealso [read.gvr()] to build the table, [gvr_summary()] to summarise the filtered
#'   variants.
#' @family germlinevaR
#' @author germlinevaR authors
#'
#' @examples
#' \dontrun{
#' gvr <- read.gvr("/path/to/vcf_folder")
#'
#' ## Default pipeline: rare variants + clinically relevant + called genotypes:
#' gvr_clean <- gvr_filter(gvr)
#'
#' ## Add protein-coding biotype filter:
#' gvr_filter(gvr, biotype_keep = c("protein_coding", "protein_coding_LoF"))
#'
#' ## Only the rarity filter on gnomAD exome AF, nothing else:
#' gvr_filter(gvr, gnomADe_AF = 0.001, AF = NULL, ABraOM_AF = NULL,
#'            clin_sig_terms = NULL, gt_exclude = NULL,
#'            vc_nonSyn = FALSE, genes = NULL)
#'
#' ## Pathogenic-only, protein-coding:
#' gvr_filter(gvr, clin_sig_terms = c("pathogenic", "likely_pathogenic"),
#'            biotype_keep = "protein_coding")
#'
#' ## Remove benign annotations (including likely_benign and compound entries):
#' gvr_filter(gvr, remove_benign = TRUE)
#'
#' ## Keep only protein-altering variants and a gene panel:
#' gvr_filter(gvr, vc_nonSyn = TRUE, genes = c("TP53", "BRCA1", "BRCA2"))
#' }
#'
#' @importFrom data.table as.data.table
#' @importFrom openxlsx createWorkbook
#' @export
gvr_filter <- function(gvr,
                       gnomADe_AF = 0.01,
                       AF = 0.01,
                       ABraOM_AF = 0.01,
                       gnomADe_AF_keep_missing = TRUE,
                       AF_keep_missing = TRUE,
                       ABraOM_AF_keep_missing = TRUE,
                       clin_sig_terms = c("likely_pathogenic",
                                          "pathogenic",
                                          "uncertain_significance"),
                       clin_sig_keep_missing = TRUE,
                       remove_benign = FALSE,
                       biotype_keep = NULL,
                       gt_exclude = c("0", "0/0"),
                       vc_nonSyn = FALSE,
                       missense_only = FALSE,  # vN+5: strict missense-only filter (Variant_Classification == "Missense_Mutation")
                       genes = NULL,
                       save_excel = FALSE,
                       out_dir = NULL,
                       file_prefix = "gvr_filter",
                       verbose = TRUE) {

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("gvr_filter requires the 'data.table' package.")
  }
  `%notin%` <- function(x, table) !(x %in% table)

  # --- Work on a copy; never mutate the caller's object -----------------------
  dt <- data.table::as.data.table(gvr)   # copies if input is a data.frame/data.table
  n_in_total <- nrow(dt)

  # --- Helper: detect missing (NA OR empty string) ----------------------------
  .is_missing <- function(v) is.na(v) | v == ""

  # --- Helper: verbose per-step logger ----------------------------------------
  .log_step <- function(label, before, after) {
    if (isTRUE(verbose)) {
      removed <- before - after
      pct <- if (before > 0) 100 * removed / before else 0
      message(sprintf("  [%-26s] %8d -> %8d   (removed %7d, %5.1f%%)",
                      label, before, after, removed, pct))
    }
  }

  # --- Standard vc_nonSyn classes (same as read.gvr) --------------------------
  .vc_nonSyn_default <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site",
                           "Translation_Start_Site", "Nonsense_Mutation",
                           "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins",
                           "Missense_Mutation")

  # --- AF filter specs: each is an independent, NULL-disabled argument --------
  #     (kept in one table so all three columns share identical logic) ---------
  af_specs <- list(
    list(col = "gnomADe_AF",        thr = gnomADe_AF,        keep_miss = isTRUE(gnomADe_AF_keep_missing)),
    list(col = "AF",                thr = AF,                keep_miss = isTRUE(AF_keep_missing)),
    list(col = "ABraOM_AF",         thr = ABraOM_AF,         keep_miss = isTRUE(ABraOM_AF_keep_missing))
  )
  af_active <- Filter(function(s) !is.null(s$thr), af_specs)

  # --- Column-existence guard: only require columns for ACTIVE filters ---------
  # For AF columns (gnomADe_AF, AF, ABraOM_AF), if the column is missing we
  # auto-disable that filter with a warning (the column may be absent because
  # e.g. add_abraom = FALSE was passed to read.gvr). For all other filter
  # columns, a missing column is a hard error (the user explicitly asked for
  # that filter but the data cannot support it).
  needed <- character(0)
  if (!is.null(clin_sig_terms) && length(clin_sig_terms) > 0) needed <- c(needed, "CLIN_SIG")
  if (isTRUE(remove_benign))                                  needed <- c(needed, "CLIN_SIG")
  if (!is.null(biotype_keep)   && length(biotype_keep)   > 0) needed <- c(needed, "BIOTYPE")
  if (!is.null(gt_exclude)     && length(gt_exclude)     > 0) needed <- c(needed, "GT")
  if (!identical(vc_nonSyn, FALSE))                            needed <- c(needed, "Variant_Classification")
  if (isTRUE(missense_only))                                   needed <- c(needed, "Variant_Classification")
  if (!is.null(genes)         && length(genes)          > 0) needed <- c(needed, "Hugo_Symbol")
  needed <- unique(needed)
  missing_cols <- needed[needed %notin% names(dt)]
  if (length(missing_cols) > 0) {
    stop(sprintf("gvr_filter: required column(s) not found for the requested filters: %s",
                 paste(missing_cols, collapse = ", ")))
  }

  # Auto-disable AF filters whose columns are absent (e.g. ABraOM_AF when
  # add_abraom = FALSE was used in read.gvr). Warn so the user knows.
  af_active <- Filter(function(s) {
    if (s$col %notin% names(dt)) {
      warning(sprintf("gvr_filter: column '%s' not found in data; skipping %s filter. (Set %s = NULL to silence this warning.)",
                       s$col, s$col, s$col))
      FALSE
    } else TRUE
  }, af_active)

  if (isTRUE(verbose)) {
    message(sprintf("gvr_filter: %d rows in", n_in_total))
  }

  # ============================================================================
  # 1-3. Allele-frequency (rare-variant) filters (gnomADe_AF, AF, ABraOM)
  # ============================================================================
  for (s in af_active) {
    thr <- s$thr
    if (!is.numeric(thr) || length(thr) != 1L) {
      stop(sprintf("gvr_filter: '%s' must be a single numeric threshold (or NULL to disable).", s$col))
    }
    before <- nrow(dt)
    raw  <- dt[[s$col]]
    miss <- .is_missing(raw)
    x    <- suppressWarnings(as.numeric(raw))
    below <- !is.na(x) & x < thr
    if (s$keep_miss) {
      keep <- miss | below
    } else {
      keep <- !miss & below
    }
    dt <- dt[keep]
    .log_step(sprintf("AF %s<%g%s", s$col, thr, if (s$keep_miss) " (keep NA)" else ""),
              before, nrow(dt))
  }

  # ============================================================================
  # 4. Clinical significance (substring, case-insensitive, optional keep-missing)
  # ============================================================================
  if (!is.null(clin_sig_terms) && length(clin_sig_terms) > 0) {
    before <- nrow(dt)
    cs   <- dt[["CLIN_SIG"]]
    miss <- .is_missing(cs)
    pat  <- paste(clin_sig_terms, collapse = "|")
    hit  <- grepl(pat, cs, ignore.case = TRUE) & !miss
    keep <- hit | (isTRUE(clin_sig_keep_missing) & miss)
    dt <- dt[keep]
    .log_step(sprintf("CLIN_SIG match%s", if (isTRUE(clin_sig_keep_missing)) "|NA" else ""),
              before, nrow(dt))
  }

  # ============================================================================
  # 5. Remove benign (substring match on CLIN_SIG, after clin_sig_terms filter)
  # ============================================================================
  if (isTRUE(remove_benign)) {
    before <- nrow(dt)
    cs <- dt[["CLIN_SIG"]]
    # Remove any row whose CLIN_SIG contains "benign" (case-insensitive).
    # This catches: benign, likely_benign, uncertain_significance&likely_benign, etc.
    # Missing/blank CLIN_SIG is NOT removed (benign is absent, so keep).
    has_benign <- grepl("benign", cs, ignore.case = TRUE) & !.is_missing(cs)
    dt <- dt[!has_benign]
    .log_step("CLIN_SIG remove benign", before, nrow(dt))
  }

  # ============================================================================
  # 6. Biotype keep-set
  # ============================================================================
  if (!is.null(biotype_keep) && length(biotype_keep) > 0) {
    before <- nrow(dt)
    keep <- dt[["BIOTYPE"]] %in% biotype_keep
    dt <- dt[keep]
    .log_step("BIOTYPE keep-set", before, nrow(dt))
  }

  # ============================================================================
  # 7. Genotype exclusion
  # ============================================================================
  if (!is.null(gt_exclude) && length(gt_exclude) > 0) {
    before <- nrow(dt)
    keep <- dt[["GT"]] %notin% gt_exclude
    dt <- dt[keep]
    .log_step("GT exclude", before, nrow(dt))
  }

  # ============================================================================
  # 8. Variant classification (vc_nonSyn)
  # ============================================================================
  if (!identical(vc_nonSyn, FALSE)) {
    vc_keep <- if (isTRUE(vc_nonSyn)) .vc_nonSyn_default else as.character(vc_nonSyn)
    vc_keep <- vc_keep[!is.na(vc_keep) & nzchar(vc_keep)]
    if (length(vc_keep) > 0L) {
      before <- nrow(dt)
      vc <- dt[["Variant_Classification"]]
      # Remove rows with missing/blank Variant_Classification when filter is active
      keep <- !.is_missing(vc) & vc %in% vc_keep
      dt <- dt[keep]
      .log_step("Variant_Classification keep", before, nrow(dt))
    }
  }

  # ============================================================================
  # 8b. vN+5: strict missense-only filter
  # ============================================================================
  # Independent of `vc_nonSyn`. When TRUE, keeps ONLY rows whose
  # Variant_Classification is "Missense_Mutation" (a strict subset of the 9
  # protein-altering classes that `vc_nonSyn = TRUE` keeps). The two filters
  # compose: TRUE+TRUE yields the same result as missense_only alone.
  if (isTRUE(missense_only)) {
    before <- nrow(dt)
    vc <- dt[["Variant_Classification"]]
    if (is.null(vc)) {
      stop("gvr_filter: missense_only=TRUE but 'Variant_Classification' ",
           "column is missing.", call. = FALSE)
    }
    keep <- !.is_missing(vc) & vc == "Missense_Mutation"
    dt <- dt[keep]
    .log_step("missense_only", before, nrow(dt))
  }

  # ============================================================================
  # 9. Gene subset (exact, case-insensitive)
  # ============================================================================
  if (!is.null(genes) && length(genes) > 0L) {
    genes_chr <- as.character(genes)
    genes_chr <- genes_chr[!is.na(genes_chr) & nzchar(genes_chr)]
    if (length(genes_chr) > 0L) {
      before <- nrow(dt)
      keep <- toupper(trimws(as.character(dt[["Hugo_Symbol"]]))) %in% toupper(genes_chr)
      dt <- dt[keep]
      .log_step("Hugo_Symbol gene subset", before, nrow(dt))
    }
  }

  if (isTRUE(verbose)) {
    kept <- nrow(dt)
    message(sprintf("gvr_filter: %d rows out (%.1f%% of input, %d removed total)",
                    kept, if (n_in_total > 0) 100 * kept / n_in_total else 0,
                    n_in_total - kept))
  }

  # ============================================================================
  # 10. Optional Excel export of the FILTERED table  ->  <out_dir>/<file_prefix>.xlsx
  # ============================================================================
  if (isTRUE(save_excel)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      warning("gvr_filter: 'openxlsx' not installed; skipping Excel export.")
    } else {
      if (is.null(out_dir)) out_dir <- "."
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      xlsx_name  <- sprintf("%s.xlsx", file_prefix)
      final_xlsx <- file.path(out_dir, xlsx_name)
      if (file.exists(final_xlsx) && isTRUE(verbose))
        message(sprintf("  Overwriting existing Excel: %s", final_xlsx))
      if (nrow(dt) > 1000000L)
        warning(sprintf("gvr_filter: filtered table has %d rows; Excel's per-sheet limit is 1,048,576 rows.",
                        nrow(dt)))
      wb <- openxlsx::createWorkbook()
      hs <- openxlsx::createStyle(textDecoration = "bold", halign = "center")
      openxlsx::addWorksheet(wb, "Filtered")
      openxlsx::writeData(wb, "Filtered", as.data.frame(dt), headerStyle = hs)
      openxlsx::freezePane(wb, "Filtered", firstRow = TRUE)
      openxlsx::setColWidths(wb, "Filtered", cols = seq_len(ncol(dt)), widths = "auto")
      tmp_xlsx <- file.path(tempdir(), xlsx_name)
      wrote_ok <- tryCatch({ openxlsx::saveWorkbook(wb, tmp_xlsx, overwrite = TRUE); TRUE },
                           error = function(e) {
                             warning(sprintf("gvr_filter: Excel write failed: %s", conditionMessage(e))); FALSE })
      if (wrote_ok) {
        system2("cp", c(shQuote(tmp_xlsx), shQuote(final_xlsx)))
        sz <- suppressWarnings(file.info(final_xlsx)$size)
        if (is.na(sz) || sz == 0) {
          warning(sprintf("gvr_filter: copy to '%s' may have failed; Excel left at '%s'.",
                          final_xlsx, tmp_xlsx))
          final_xlsx <- tmp_xlsx
        }
        if (isTRUE(verbose)) message(sprintf("  Excel written: %s", final_xlsx))
      }
    }
  }

  dt[]
}
