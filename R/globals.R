# germlinevaR: package-level global-variable declarations
#
# data.table's non-standard evaluation and ggplot2's aes() reference column
# names as bare symbols, which `R CMD check` flags as "no visible binding for
# global variable". Declaring them here (package-scoped, registered at build
# time) silences those NOTEs for every function in the package. This is the
# idiomatic CRAN placement: one dedicated R/globals.R rather than per-file blocks.

#' @importFrom utils globalVariables
#' @keywords internal
utils::globalVariables(c(
  # --- data.table special symbols (used inside dt[...] across functions) ----
  ".", ".N", ".SD", ".__sample__", ".sr",
  # --- data.table column-vector prefix (used as `dt[, ..cols]` in gvr_summary)
  "..cols",

  # --- read.gvr(): temp/internal columns created via := ---------------------
  ".rs", ".ref", ".alt",
  # --- read.gvr(): MAF columns referenced bare inside data.table calls ------
  "Genotype", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "ABraOM_AF",
  "dbSNP_RS", "Reference_Allele",
  # --- read.gvr(): ABraOM lookup-table columns (built then keyed/joined) ----
  "avsnp147", "Ref", "Alt", "Frequencies", "rs", "ref", "alt", "af", "x.af",

  # --- gvr_summary(): section columns + nested-renderer aes() symbols -------
  "Hugo_Symbol", "Total", "Variant_Classification", "CLIN_SIG",
  "Category", "Sample", "n",
  # --- gvr_summary(): top-rsID aggregation + drill-down detail columns ------
  ".__rs__", "Chromosome", "Start_Position", "IMPACT",

  # --- gvr_plot(): matrix/summary columns ------------------------------
  "n_var", "n_samp", "N"
))
