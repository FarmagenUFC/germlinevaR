#' Filter a read.gvr table down to candidate novel variants
#'
#' @description
#' Returns only the rows of a `read.gvr()` / `gvr_filter()` table that have NO
#' database evidence in any of the standard catalogues: no rsID in `dbSNP_RS`,
#' AND no allele frequency in `gnomADe_AF`, `AF`, or `ABraOM_AF`. This is the
#' canonical "show me only the novel variants" subsetter and complements
#' [gvr_filter()] (which uses presence of AF as a continuous filter).
#'
#' @details
#' A column value counts as "missing/empty" iff `is.na(x) | x == ""`. This is
#' the same convention used throughout the package (see [gvr_filter()] and
#' [gvr_summary()]): `read.gvr()` writes `""` for absent values, but NA can
#' creep in if upstream code re-coerced columns.
#'
#' A row is kept iff ALL FOUR of the following columns are missing/empty:
#' \itemize{
#'   \item `dbSNP_RS`     - rsID from VEP `Existing_variation`
#'   \item `gnomADe_AF`   - gnomAD exome allele frequency
#'   \item `AF`           - gnomAD genome / VCF allele frequency
#'   \item `ABraOM_AF`    - ABraOM SABE-609 Brazilian-cohort allele frequency
#' }
#' If any of the four columns is absent from the input (e.g. ABraOM disabled
#' at read time), the function treats that column as "all missing" - i.e. it
#' does not exclude any row on the basis of the missing column. The audit log
#' will say `(column not present, skipped)` for transparency.
#'
#' The cascade order is fixed (dbSNP_RS -> gnomADe_AF -> AF -> ABraOM_AF) so
#' the verbose audit is reproducible across runs.
#'
#' The input is never modified: filtering operates on an internal `data.table`
#' copy and the output has the same column set and order as the input.
#'
#' @param gvr data.table / data.frame from [read.gvr()] (or compatible, e.g. the
#'   output of [gvr_filter()]). Filtered in a copy; the input object is not
#'   modified.
#' @param verbose logical(1). If `TRUE` (default), print a step-by-step audit
#'   line for each of the 4 cascaded filters showing rows-in -> rows-out and
#'   percentage retained, ending with a one-line summary.
#'
#' @return A `data.table` with the same columns as `gvr` but only the rows that
#'   pass all four "no database evidence" checks. If none of the 4 columns is
#'   present in `gvr`, the output equals the input and a warning is issued.
#'
#' @seealso [read.gvr()], [gvr_filter()], [gvr_summary()]
#'
#' @examples
#' ## Load the shipped example table and find candidate novel variants
#' gvr <- readRDS(system.file("extdata", "example_gvr.rds",
#'                            package = "germlinevaR"))
#' nov <- gvr_novel(gvr, verbose = FALSE)
#' dim(nov)
#'
#' ## Sanity-check: every kept row really is novel
#' stopifnot(all(is.na(nov$dbSNP_RS)   | nov$dbSNP_RS   == ""))
#' stopifnot(all(is.na(nov$gnomADe_AF) | nov$gnomADe_AF == ""))
#'
#' ## Combine with gvr_filter() to restrict to filtered novel variants
#' filt <- gvr_filter(gvr, ABraOM_AF = NULL, verbose = FALSE)
#' nov_filt <- gvr_novel(filt, verbose = FALSE)
#' dim(nov_filt)
#' @importFrom data.table as.data.table copy setDT is.data.table
#' @export
gvr_novel <- function(gvr, verbose = TRUE) {

  # ---- Nested helpers ----
  .is_missing <- function(x) is.na(x) | x == ""

  .log_step <- function(label, before, after, present) {
    if (!isTRUE(verbose)) return(invisible())
    if (!isTRUE(present)) {
      message(sprintf("  [%-16s] (column not present, skipped)", label))
      return(invisible())
    }
    pct <- if (before == 0L) 0 else 100 * after / before
    message(sprintf("  [%-16s] %6d -> %6d (kept %6d, %5.1f%%)",
                    label, before, after, after, pct))
  }

  # ---- Input validation ----
  if (!is.data.frame(gvr)) {
    stop("gvr_novel: 'gvr' must be a data.frame / data.table.")
  }
  if (nrow(gvr) == 0L) {
    if (isTRUE(verbose)) message("gvr_novel: input has 0 rows; returning empty table.")
    return(data.table::as.data.table(gvr))
  }

  # ---- Work on a data.table copy (do not modify input) ----
  dt <- if (data.table::is.data.table(gvr)) data.table::copy(gvr) else data.table::as.data.table(gvr)

  # ---- Cascade columns + audit-friendly labels (fixed order) ----
  cascade <- list(
    list(col = "dbSNP_RS",   label = "no dbSNP_RS"),
    list(col = "gnomADe_AF", label = "no gnomADe_AF"),
    list(col = "AF",         label = "no AF"),
    list(col = "ABraOM_AF",  label = "no ABraOM_AF")
  )

  present_cols <- vapply(cascade, function(s) s$col %in% names(dt), logical(1L))
  if (!any(present_cols)) {
    warning("gvr_novel: none of dbSNP_RS / gnomADe_AF / AF / ABraOM_AF found in input; ",
            "returning input unchanged.")
    if (isTRUE(verbose)) message(sprintf("gvr_novel: %d rows in, %d rows out (filter inactive).",
                                          nrow(dt), nrow(dt)))
    return(dt)
  }

  if (isTRUE(verbose)) message(sprintf("gvr_novel: %d rows in\n", nrow(dt)))

  in_rows <- nrow(dt)

  # ---- Apply cascade ----
  for (i in seq_along(cascade)) {
    spec   <- cascade[[i]]
    before <- nrow(dt)
    present <- spec$col %in% names(dt)
    if (present) {
      dt <- dt[.is_missing(get(spec$col))]
    }
    after <- nrow(dt)
    .log_step(spec$label, before, after, present)
  }

  # ---- Closing summary ----
  if (isTRUE(verbose)) {
    out_rows <- nrow(dt)
    pct_kept <- if (in_rows == 0L) 0 else 100 * out_rows / in_rows
    n_known_genes <- if ("Hugo_Symbol" %in% names(dt))
      length(unique(dt$Hugo_Symbol[!.is_missing(dt$Hugo_Symbol)])) else NA_integer_
    message(sprintf("\ngvr_novel: %d rows out (%.1f%% of input%s)",
                    out_rows, pct_kept,
                    if (!is.na(n_known_genes)) sprintf(", %d distinct known genes", n_known_genes) else ""))
  }

  dt
}
