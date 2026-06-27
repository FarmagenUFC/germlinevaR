# gvr_panels.R - Disease gene-panel registry for germlinevaR
#
# Provides two exported helpers:
#   gvr_list_panels()       - sorted character vector of available disease names
#   gvr_panel_genes(panel)  - Hugo_Symbol vector for one disease (case-insensitive,
#                             accepts underscores as aliases for spaces, and a small
#                             panel-name alias table for legacy / typo / synonym keys)
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
    "men1" = sort(unique(toupper(c(
      "AIP", "AP2S1", "CASR", "CDC73", "CDKN1A", "CDKN1B",
      "CDKN2A", "CDKN2C", "GCM2", "GNA11", "MAX", "MEN1",
      "RET"
    )))),

    "acromegaly" = sort(unique(toupper(c(
      "AIP", "CDKN1A", "CDKN1B", "CDKN2A", "GNAS", "GPR101",
      "HESX1", "LHX3", "LHX4", "MEN1", "NF1", "PDE4A",
      "PDE4B", "PRKACB", "PROP1", "RB1", "SDHA", "SDHD",
      "TERT", "TP53", "TSC1", "USP48", "USP8"
    )))),

    "pheochromocytoma" = sort(unique(toupper(c(
      "ATM", "DLST", "EGLN1", "EGLN2", "EPAS1", "FH",
      "HRAS", "KIF1B", "MAX", "MDH2", "MEN1", "MERTK",
      "MET", "NF1", "RET", "SDHA", "SDHAF2", "SDHB",
      "SDHC", "SDHD", "SLC25A11", "TMEM127", "TP53", "VHL"
    )))),

    "hereditary cancer" = sort(unique(toupper(c(
      "ALK", "APC", "ATM", "AXIN2", "BAP1", "BARD1",
      "BLM", "BMPR1A", "BRCA1", "BRCA2", "BRIP1", "BUB1B",
      "CDC73", "CDH1", "CDK4", "CDKN1B", "CDKN2A", "CEP57",
      "CHEK2", "CTNNA1", "CYLD", "DICER1", "EGFR", "EGLN1",
      "EPCAM", "EXT1", "EXT2", "FAN1", "FH", "FLCN",
      "GALNT12", "GREM1", "HOXB13", "HRAS", "KIF1B", "KIT",
      "LZTR1", "MAX", "MC1R", "MEN1", "MET", "MITF",
      "MLH1", "MLH3", "MSH2", "MSH3", "MSH6", "MUTYH",
      "NF1", "NF2", "NTHL1", "PALB2", "PDGFRA", "PHOX2B",
      "PIK3CA", "PMS2", "POLD1", "POLE", "POT1", "PRKAR1A",
      "PTCH1", "PTCH2", "PTEN", "RAD51C", "RAD51D", "RB1",
      "RET", "RHBDF2", "RNF43", "SDHA", "SDHAF2", "SDHB",
      "SDHC", "SDHD", "SMAD4", "SMARCA4", "SMARCB1", "SMARCE1",
      "STK11", "SUFU", "TERT", "TMEM127", "TP53", "TSC1",
      "TSC2", "VHL", "WT1"
    )))),

    "gist" = sort(unique(toupper(c(
      "BRAF", "KIT", "NF1", "NTRK1", "NTRK2", "NTRK3",
      "PDGFRA", "SDHA", "SDHB", "SDHC", "SDHD"
    )))),

    "lynch syndrome" = sort(unique(toupper(c(
      "EPCAM", "MLH1", "MSH2", "MSH6", "PMS2"
    )))),

    "li-fraumeni syndrome" = sort(unique(toupper(c(
      "ARHGAP30", "CHEK2", "PIK3CA", "TP53"
    )))),

    "hereditary gastric cancer" = sort(unique(toupper(c(
      "APC", "ATM", "BLM", "BMPR1A", "BRCA1", "BRCA2",
      "CDH1", "CTNNA1", "EPCAM", "GREM1", "KIT", "MEN1",
      "MLH1", "MSH2", "MSH6", "MUTYH", "NF1", "PALB2",
      "PDGFRA", "PMS2", "POLD1", "POLE", "PTEN", "RNF43",
      "SMAD4", "STK11", "TP53", "VHL"
    )))),

    "hereditary colorectal cancer" = sort(unique(toupper(c(
      "APC", "ATM", "AXIN2", "BMPR1A", "CDH1", "CHEK2",
      "EPCAM", "GALNT12", "GREM1", "MLH1", "MSH2", "MSH3",
      "MSH6", "MUTYH", "NTHL1", "PMS2", "POLD1", "POLE",
      "PTEN", "RNF43", "SMAD4", "STK11", "TP53"
    )))),

    "familial adenomatous polyposis" = sort(unique(toupper(c(
      "APC", "AXIN2", "BMPR1A", "EPCAM", "MBD4", "MLH1",
      "MSH2", "MSH3", "MSH6", "MUTYH", "NTHL1", "PMS2",
      "POLD1", "POLE", "PTEN", "RNF43", "SMAD4", "STK11",
      "TP53"
    )))),

    "hereditary melanoma cancer" = sort(unique(toupper(c(
      "BAP1", "BRCA1", "BRCA2", "CDK4", "CDKN2A", "CHEK2",
      "EPCAM", "MC1R", "MITF", "MLH1", "MSH2", "MSH6",
      "PMS2", "POT1", "PTCH1", "PTEN", "RB1", "TERT",
      "TP53", "TYR"
    )))),

    "hereditary prostate cancer" = sort(unique(toupper(c(
      "ATM", "ATR", "BRCA1", "BRCA2", "BRIP1", "CHEK2",
      "EPCAM", "HOXB13", "MLH1", "MRE11", "MSH2", "MSH6",
      "NBN", "PALB2", "PMS2", "PTEN", "RAD51C", "RAD51D",
      "TP53"
    )))),

    "hereditary breast and ovarian cancer" = sort(unique(toupper(c(
      "ATM", "BARD1", "BRCA1", "BRCA2", "BRIP1", "CDH1",
      "CHEK2", "EPCAM", "MLH1", "MRE11", "MSH2", "MSH6",
      "NBN", "NF1", "PALB2", "PMS1", "PMS2", "PTEN",
      "RAD50", "RAD51C", "RAD51D", "SMARCA4", "STK11", "TP53",
      "XRCC2"
    )))),

    "breast cancer" = sort(unique(toupper(c(
      "AKT1", "ATM", "BARD1", "BRCA1", "BRCA2", "BRIP1",
      "CDH1", "CHEK2", "EPCAM", "ERBB2", "ESR1", "MLH1",
      "MSH2", "MSH6", "NBN", "NF1", "PALB2", "PIK3CA",
      "PMS2", "PTEN", "RAD51C", "RAD51D", "STK11", "TP53"
    )))),

    "breast cancer somatic" = sort(unique(toupper(c(
      "AKT1", "BRCA1", "BRCA2", "ERBB2", "ESR1", "PIK3CA"
    ))))
  )
}


# Internal: panel-name alias map. Maps legacy or alternative spellings
# (typos, long-form synonyms) to canonical keys. Applied AFTER
# .gvr_panel_normalise() basic normalisation. Returning an empty vector
# is fine; lookups just fall through.
.gvr_panel_alias_map <- function() {
  c(
    "gastrointestinal stromal tumor" = "gist",
    "pheocromocytoma"                = "pheochromocytoma",
    "breast cancer somatic panel"    = "breast cancer somatic"
  )
}


# Internal: canonicalise a user-supplied panel name. Lowercase, trim,
# collapse underscores to spaces, squeeze multi-space to single space,
# then apply the panel-name alias map.
# Returns "" if the input does not normalise to a non-empty string.
.gvr_panel_normalise <- function(x) {
  if (is.null(x)) return(character(0L))
  s <- as.character(x)
  s <- s[!is.na(s)]
  if (length(s) == 0L) return(character(0L))
  s <- trimws(s)
  s <- tolower(s)
  s <- gsub("_+", " ", s, perl = TRUE)
  s <- gsub("[[:space:]]+", " ", s, perl = TRUE)
  amap <- .gvr_panel_alias_map()
  hit  <- s %in% names(amap)
  s[hit] <- amap[s[hit]]
  s
}


# Internal: a single helpful catalog string used in error messages.
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
#' as spaces. A small alias table is also recognised: e.g.
#' `"gastrointestinal stromal tumor"` is accepted as a synonym for
#' `"gist"`, and the common typo `"pheocromocytoma"` resolves to
#' `"pheochromocytoma"`.
#'
#' Panels currently shipped:
#'
#' \describe{
#'   \item{`"men1"`}{13 genes.}
#'   \item{`"acromegaly"`}{23 genes.}
#'   \item{`"pheochromocytoma"`}{24 genes. (alias: `"pheocromocytoma"`)}
#'   \item{`"hereditary cancer"`}{87 genes.}
#'   \item{`"gist"`}{11 genes. (alias: `"gastrointestinal stromal tumor"`)}
#'   \item{`"lynch syndrome"`}{5 genes.}
#'   \item{`"li-fraumeni syndrome"`}{4 genes.}
#'   \item{`"hereditary gastric cancer"`}{28 genes.}
#'   \item{`"hereditary colorectal cancer"`}{23 genes.}
#'   \item{`"familial adenomatous polyposis"`}{19 genes.}
#'   \item{`"hereditary melanoma cancer"`}{20 genes.}
#'   \item{`"hereditary prostate cancer"`}{19 genes.}
#'   \item{`"hereditary breast and ovarian cancer"`}{25 genes.}
#'   \item{`"breast cancer"`}{24 genes.}
#'   \item{`"breast cancer somatic"`}{6 genes. (alias: `"breast cancer somatic panel"`)}
#' }
#'
#' Use [gvr_panel_genes()] to retrieve the gene vector for a specific
#' disease, or pass the panel name (or vector of names) to
#' [read.gvr()] via its `panel` argument to filter a gvr table.
#'
#' @examples
#' gvr_list_panels()
#'
#' ## The returned panel names can be passed to read.gvr() to subset
#' ## the gvr table at read time. Equivalent post-hoc filter on the
#' ## pre-parsed example table (instantaneous; no VCF re-parse):
#' gvr <- readRDS(system.file("extdata", "example_gvr.rds",
#'                            package = "germlinevaR"))
#' panel_genes <- gvr_panel_genes("breast cancer")
#' nrow(gvr[gvr$Hugo_Symbol %in% panel_genes, ])
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
#'   `"breast_cancer"` is equivalent to `"breast cancer"`; selected
#'   aliases like `"gastrointestinal stromal tumor"` -> `"gist"` are
#'   recognised). Must resolve to exactly one entry in the registry.
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
#' gvr_panel_genes("gastrointestinal stromal tumor")  # alias of "gist"
#' gvr_panel_genes("hereditary prostate cancer")
#'
#' ## Combine multiple panels and post-hoc filter the pre-parsed example
#' ## table (equivalent to read.gvr(..., panel = c(...)) but instantaneous):
#' gvr <- readRDS(system.file("extdata", "example_gvr.rds",
#'                            package = "germlinevaR"))
#' multi_panel <- unique(c(gvr_panel_genes("breast cancer"),
#'                         gvr_panel_genes("hereditary prostate cancer")))
#' nrow(gvr[gvr$Hugo_Symbol %in% multi_panel, ])
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

