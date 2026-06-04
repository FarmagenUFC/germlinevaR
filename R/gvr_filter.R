#' Modular, individually-toggleable filtering of a read.gvr MAF
#'
#' @description
#' Applies a set of independent variant filters to the MAF-style table produced by
#' [read.gvr()]. Each distinct filter is its own argument; setting an argument to
#' `NULL` disables that filter entirely (no rows removed by it). With all defaults,
#' `gvr_filter(maf)` reproduces the canonical rare / clinically-relevant /
#' protein-coding / called-genotype germline pipeline.
#'
#' @details
#' Filters are applied in a fixed order; each step operates on the survivors of the
#' previous one:
#' \enumerate{
#'   \item gnomAD exome AF        - `gnomADe_AF`        (+ `gnomADe_AF_keep_missing`)
#'   \item gnomAD genome / VCF AF - `AF`                (+ `AF_keep_missing`)
#'   \item ABraOM (SABE 609) AF   - `ABraOM_SABE609_AF` (+ `ABraOM_SABE609_AF_keep_missing`)
#'   \item Clinical significance  - `clin_sig_terms`    (+ `clin_sig_keep_missing`)
#'   \item Biotype                - `biotype_keep`
#'   \item Genotype exclusion     - `gt_exclude`
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
#' @param maf data.table / data.frame from [read.gvr()] (or compatible). Filtered in a
#'   copy; the input object is not modified.
#'
#' @param gnomADe_AF Numeric upper threshold for the gnomAD exome AF column
#'   `gnomADe_AF` (keep rows with AF < threshold). `NULL` disables this filter.
#'   Default 0.01.
#' @param AF Numeric upper threshold for the `AF` column (gnomAD genome / VCF allele
#'   frequency). `NULL` disables this filter. Default 0.01.
#' @param ABraOM_SABE609_AF Numeric upper threshold for the Brazilian-cohort column
#'   `ABraOM_SABE609_AF`. `NULL` disables this filter. Default 0.01.
#' @param gnomADe_AF_keep_missing Logical; if TRUE, keep rows whose `gnomADe_AF` is
#'   missing (NA or ""); if FALSE (default) drop them. Matches the literal dplyr example
#'   (missing -> dropped). Ignored when `gnomADe_AF` is NULL.
#' @param AF_keep_missing Logical; missing-value handling for the `AF` filter. FALSE
#'   (default) drops missing. Ignored when `AF` is NULL.
#' @param ABraOM_SABE609_AF_keep_missing Logical; missing-value handling for the ABraOM
#'   filter. FALSE (default) drops missing. Set TRUE to retain variants absent from the
#'   Brazilian cohort (where absence often means "not catalogued", not "common").
#'   Ignored when `ABraOM_SABE609_AF` is NULL.
#'
#' @param clin_sig_terms Character vector of clinical-significance terms to keep
#'   (substring, case-insensitive, OR-combined). `NULL` disables the CLIN_SIG filter.
#'   Default: c("likely_pathogenic","pathogenic","uncertain_significance").
#' @param clin_sig_keep_missing Logical; if TRUE (default) rows with missing CLIN_SIG
#'   (NA/"") are kept. Only relevant when `clin_sig_terms` is non-NULL.
#' @param biotype_keep Character vector of BIOTYPE values to keep (exact match via %in%).
#'   `NULL` disables the biotype filter.
#'   Default: c("protein_coding","protein_coding_LoF").
#' @param gt_exclude Character vector of GT values to remove (exact match). `NULL`
#'   disables the genotype filter. Default: c("0","0/0").
#' @param verbose Logical; if TRUE (default) print a per-filter breakdown (rows in -> out
#'   and rows removed by each active step).
#'
#' @return A `data.table` of the surviving rows, with the same columns as the input.
#'   A plain `data.frame` input is returned as a `data.table`. The input object is not
#'   modified. With `verbose = TRUE`, a per-filter breakdown (rows in -> out, and rows
#'   removed by each active step) is printed as it runs.
#'
#' @seealso [read.gvr()] to build the MAF, [gvr_summary()] to summarise the filtered
#'   variants.
#' @family germlinevaR
#' @author germlinevaR authors
#'
#' @examples
#' \dontrun{
#' maf <- read.gvr("/path/to/vcf_folder")
#'
#' ## Reproduce the default rare / clinically-relevant / protein-coding pipeline:
#' maf_clean <- gvr_filter(maf)
#'
#' ## Strict gnomAD, but keep variants absent from ABraOM:
#' gvr_filter(maf, ABraOM_SABE609_AF_keep_missing = TRUE)
#'
#' ## Only the rarity filter on gnomAD exome AF, nothing else:
#' gvr_filter(maf, gnomADe_AF = 0.001, AF = NULL, ABraOM_SABE609_AF = NULL,
#'            clin_sig_terms = NULL, biotype_keep = NULL, gt_exclude = NULL)
#'
#' ## Pathogenic-only (drop uncertain_significance), exact protein_coding:
#' gvr_filter(maf, clin_sig_terms = c("pathogenic", "likely_pathogenic"),
#'            biotype_keep = "protein_coding")
#' }
#'
#' @importFrom data.table as.data.table
#' @export
gvr_filter <- function(maf,
                       gnomADe_AF = 0.01,
                       AF = 0.01,
                       ABraOM_SABE609_AF = 0.01,
                       gnomADe_AF_keep_missing = FALSE,
                       AF_keep_missing = FALSE,
                       ABraOM_SABE609_AF_keep_missing = FALSE,
                       clin_sig_terms = c("likely_pathogenic",
                                          "pathogenic",
                                          "uncertain_significance"),
                       clin_sig_keep_missing = TRUE,
                       biotype_keep = c("protein_coding", "protein_coding_LoF"),
                       gt_exclude = c("0", "0/0"),
                       verbose = TRUE) {

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("gvr_filter requires the 'data.table' package.")
  }
  `%notin%` <- function(x, table) !(x %in% table)

  # --- Work on a copy; never mutate the caller's object -----------------------
  dt <- data.table::as.data.table(maf)   # copies if input is a data.frame/data.table
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

  # --- AF filter specs: each is an independent, NULL-disabled argument --------
  #     (kept in one table so all three columns share identical logic) ---------
  af_specs <- list(
    list(col = "gnomADe_AF",        thr = gnomADe_AF,        keep_miss = isTRUE(gnomADe_AF_keep_missing)),
    list(col = "AF",                thr = AF,                keep_miss = isTRUE(AF_keep_missing)),
    list(col = "ABraOM_SABE609_AF", thr = ABraOM_SABE609_AF, keep_miss = isTRUE(ABraOM_SABE609_AF_keep_missing))
  )
  af_active <- Filter(function(s) !is.null(s$thr), af_specs)

  # --- Column-existence guard: only require columns for ACTIVE filters ---------
  needed <- character(0)
  if (length(af_active) > 0) needed <- c(needed, vapply(af_active, `[[`, "", "col"))
  if (!is.null(clin_sig_terms) && length(clin_sig_terms) > 0) needed <- c(needed, "CLIN_SIG")
  if (!is.null(biotype_keep)   && length(biotype_keep)   > 0) needed <- c(needed, "BIOTYPE")
  if (!is.null(gt_exclude)     && length(gt_exclude)     > 0) needed <- c(needed, "GT")
  needed <- unique(needed)
  missing_cols <- needed[needed %notin% names(dt)]
  if (length(missing_cols) > 0) {
    stop(sprintf("gvr_filter: required column(s) not found for the requested filters: %s",
                 paste(missing_cols, collapse = ", ")))
  }

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
  # 5. Biotype keep-set
  # ============================================================================
  if (!is.null(biotype_keep) && length(biotype_keep) > 0) {
    before <- nrow(dt)
    keep <- dt[["BIOTYPE"]] %in% biotype_keep
    dt <- dt[keep]
    .log_step("BIOTYPE keep-set", before, nrow(dt))
  }

  # ============================================================================
  # 6. Genotype exclusion
  # ============================================================================
  if (!is.null(gt_exclude) && length(gt_exclude) > 0) {
    before <- nrow(dt)
    keep <- dt[["GT"]] %notin% gt_exclude
    dt <- dt[keep]
    .log_step("GT exclude", before, nrow(dt))
  }

  if (isTRUE(verbose)) {
    kept <- nrow(dt)
    message(sprintf("gvr_filter: %d rows out (%.1f%% of input, %d removed total)",
                    kept, if (n_in_total > 0) 100 * kept / n_in_total else 0,
                    n_in_total - kept))
  }

  dt[]
}
