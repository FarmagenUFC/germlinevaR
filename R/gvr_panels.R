# gvr_panels.R - Disease gene-panel registry for germlinevaR
#
# Provides two exported helpers:
#   gvr_list_panels()       - sorted character vector of available disease names
#   gvr_panel_genes(panel)  - Hugo_Symbol vector for one disease (case-insensitive,
#                             accepts underscores as aliases for spaces)
#
# These are consumed by read.gvr() via its `panel` argument, which unions the
# resolved gene list with the user-supplied `genes` vector before filtering.
#
# Adding a new disease: append one entry to .gvr_panel_registry() below.
# Canonical keys MUST be lowercase with spaces. Gene vectors should use
# current HGNC symbols (uppercased on storage, no duplicates).
# ----------------------------------------------------------------------------


# Internal: the registry itself. Returned as a list so callers can iterate,
# but normal users go through gvr_list_panels() / gvr_panel_genes().
# Canonical keys: lowercase, spaces. Gene vectors: HGNC symbols, uppercase,
# sorted alphabetically, unique.
.gvr_panel_registry <- function() {
  list(
    "breast cancer" = sort(unique(toupper(c(
      "AKT1", "BRCA1", "BRCA2", "ERBB2", "ESR1", "PIK3CA"
    )))),

    "hereditary prostate cancer" = sort(unique(toupper(c(
      "ATM", "ATR", "BRCA1", "BRCA2", "BRIP1", "CHEK2", "EPCAM",
      "HOXB13", "MLH1", "MRE11", "MSH2", "MSH6", "NBN", "PALB2",
      "PMS2", "PTEN", "RAD51C", "RAD51D", "TP53"
    ))))
  )
}


# Internal: canonicalise a user-supplied panel name. Lowercase, trim,
# collapse underscores to spaces, squeeze multi-space to single space.
# Returns "" if the input doesn't normalise to a non-empty string.
.gvr_panel_normalise <- function(x) {
  if (is.null(x)) return(character(0L))
  s <- as.character(x)
  s <- s[!is.na(s)]
  if (length(s) == 0L) return(character(0L))
  s <- trimws(s)
  s <- tolower(s)
  s <- gsub("_+", " ", s, perl = TRUE)
  s <- gsub("[[:space:]]+", " ", s, perl = TRUE)
  s
}


# Internal: a single helpful catalog string used in error messages.
# Looks like:
#   Available panels:
#     - breast cancer (6 genes)
#     - hereditary prostate cancer (19 genes)
.gvr_panel_catalog <- function() {
  reg <- .gvr_panel_registry()
  nm  <- sort(names(reg))
  lines <- vapply(nm, function(n) {
    sprintf("  - %s (%d genes)", n, length(reg[[n]]))
  }, character(1L), USE.NAMES = FALSE)
  paste(c("Available panels:", lines), collapse = "\n")
}


#' List Available Disease Gene Panels
#'
#' Returns the sorted character vector of every disease name accepted by
#' [gvr_panel_genes()] and by [read.gvr()]'s `panel` argument.
#'
#' @return Character vector of canonical panel names (lowercase, spaces).
#'
#' @details
#' Each canonical name corresponds to a curated Hugo_Symbol gene list.
#' Panel names are matched case-insensitively and underscores are treated
#' as spaces, so `"breast cancer"`, `"Breast Cancer"`, and
#' `"breast_cancer"` all resolve to the same entry.
#'
#' Panels currently shipped:
#'
#' \describe{
#'   \item{`"breast cancer"`}{6 genes — `AKT1`, `BRCA1`, `BRCA2`, `ERBB2`,
#'     `ESR1`, `PIK3CA`.}
#'   \item{`"hereditary prostate cancer"`}{19 genes — `ATM`, `ATR`, `BRCA1`,
#'     `BRCA2`, `BRIP1`, `CHEK2`, `EPCAM`, `HOXB13`, `MLH1`, `MRE11`, `MSH2`,
#'     `MSH6`, `NBN`, `PALB2`, `PMS2`, `PTEN`, `RAD51C`, `RAD51D`, `TP53`.}
#' }
#'
#' Use [gvr_panel_genes()] to retrieve the gene vector for a specific
#' disease, or pass the panel name (or vector of names) to
#' [read.gvr()] via its `panel` argument to filter a MAF.
#'
#' @examples
#' gvr_list_panels()
#'
#' ## Use in read.gvr() to keep only breast-cancer genes
#' \dontrun{
#' maf <- read.gvr("/path/to/folder", panel = "breast cancer")
#' }
#'
#' @seealso [gvr_panel_genes()], [read.gvr()].
#' @family germlinevaR
#' @author germlinevaR authors
#' @export
gvr_list_panels <- function() {
  sort(names(.gvr_panel_registry()))
}


#' Genes for a Disease Panel
#'
#' Returns the Hugo_Symbol vector associated with a given disease panel
#' name. The same registry is used by [read.gvr()] when the user passes a
#' `panel` argument.
#'
#' @param panel Single character string naming a disease panel
#'   (case-insensitive; underscores are treated as spaces, e.g.
#'   `"breast_cancer"` is equivalent to `"breast cancer"`). Must resolve
#'   to exactly one entry in the registry.
#'
#' @return Character vector of Hugo_Symbols (uppercase, sorted, unique).
#'
#' @details
#' If `panel` does not resolve to a known disease, an error is thrown
#' listing all available panels with their gene counts, so the caller
#' sees the catalog at the point of failure.
#'
#' @examples
#' gvr_panel_genes("breast cancer")
#' gvr_panel_genes("Breast_Cancer")            # case + underscore alias
#' gvr_panel_genes("hereditary prostate cancer")
#'
#' \dontrun{
#' ## Use in read.gvr()
#' maf <- read.gvr("/path/to/folder",
#'                 panel = c("breast cancer", "hereditary prostate cancer"))
#' }
#'
#' @seealso [gvr_list_panels()], [read.gvr()].
#' @family germlinevaR
#' @author germlinevaR authors
#' @export
gvr_panel_genes <- function(panel) {
  if (missing(panel) || is.null(panel))
    stop("gvr_panel_genes: 'panel' is required.\n", .gvr_panel_catalog(),
         call. = FALSE)
  if (length(panel) != 1L)
    stop("gvr_panel_genes: 'panel' must be a single name (length 1). ",
         "To union multiple panels, pass a vector to read.gvr(panel = ...).\n",
         .gvr_panel_catalog(),
         call. = FALSE)

  key <- .gvr_panel_normalise(panel)
  reg <- .gvr_panel_registry()
  if (length(key) == 0L || !nzchar(key) || !(key %in% names(reg))) {
    stop(sprintf("gvr_panel_genes: unknown panel '%s'.\n%s",
                 as.character(panel), .gvr_panel_catalog()),
         call. = FALSE)
  }
  reg[[key]]
}


# ---------------------------------------------------------------------------
# Internal: bulk resolver used by read.gvr() and read.gvr.snpeff().
# Accepts a (possibly NULL) character vector of panel names; returns the
# uppercased, deduplicated union of their gene lists. Empty input -> empty
# character vector. Unknown names trigger the same catalog error as
# gvr_panel_genes().
# ---------------------------------------------------------------------------
.gvr_resolve_panels <- function(panel) {
  if (is.null(panel)) return(character(0L))
  v <- as.character(panel)
  v <- v[!is.na(v) & nzchar(trimws(v))]
  if (length(v) == 0L) return(character(0L))
  unique(toupper(unlist(lapply(v, gvr_panel_genes), use.names = FALSE)))
}
