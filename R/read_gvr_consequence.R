# =============================================================================
# read_gvr_consequence.R
#
# Internal helpers used by read.gvr() for VEP consequence dispatch:
#   * `.gvr_effect_priority` - the Ensembl SO consequence severity table
#   * `.gvr_severe_term_uncached` / `.gvr_severe_term` - most-severe-term and
#       severity rank for a Consequence string ("missense_variant&splice_..."),
#       with optional per-call memoisation
#   * `.gvr_most_severe_term` / `.gvr_consequence_rank` - thin wrappers
#       around the memoised resolver
#   * `.gvr_vep_to_class_uncached` / `.gvr_vep_to_class` - VEP SO term ->
#       Variant_Classification (MAF style), with optional per-call memoisation
#
# Cache strategy: every memoised helper accepts a `cache = NULL` argument
# typed as an `environment()`. The caller (`read.gvr()`) creates fresh
# `new.env(parent = emptyenv())` envs at the start of each call and threads
# them through. This preserves Turn-3 semantics exactly (caches are scoped
# to one `read.gvr()` invocation, not to the session).
#
# All helpers are package-internal: not exported, no Rd page (@keywords
# internal + @noRd). Promoted from read.gvr()'s body in the Turn 5a refactor.
# =============================================================================

# Package-private constants -----------------------------------------------------

# Ensembl SO consequence severity priority (vcf2maf-style; lower = more severe).
# Used by `.gvr_severe_term_uncached()` to rank multi-term Consequence strings.
.gvr_effect_priority <- c(
    "transcript_ablation" = 1, "exon_loss_variant" = 2, "splice_donor_variant" = 3,
    "splice_acceptor_variant" = 3, "stop_gained" = 4, "frameshift_variant" = 5,
    "stop_lost" = 6, "start_lost" = 7, "initiator_codon_variant" = 8,
    "transcript_amplification" = 9, "protein_altering_variant" = 10,
    "missense_variant" = 11, "conservative_missense_variant" = 11,
    "rare_amino_acid_variant" = 11, "incomplete_terminal_codon_variant" = 14,
    "splice_region_variant" = 13, "splice_donor_5th_base_variant" = 13,
    "splice_donor_region_variant" = 13, "splice_polypyrimidine_tract_variant" = 13,
    "stop_retained_variant" = 15, "synonymous_variant" = 15, "coding_sequence_variant" = 16,
    "mature_miRNA_variant" = 17, "5_prime_UTR_variant" = 18,
    "5_prime_UTR_premature_start_codon_gain_variant" = 18, "3_prime_UTR_variant" = 19,
    "non_coding_transcript_exon_variant" = 20, "non_coding_exon_variant" = 20,
    "intron_variant" = 21, "non_coding_transcript_variant" = 22, "nc_transcript_variant" = 22,
    "NMD_transcript_variant" = 23, "upstream_gene_variant" = 24, "downstream_gene_variant" = 25,
    "TFBS_ablation" = 26, "TFBS_amplification" = 27, "TF_binding_site_variant" = 28,
    "regulatory_region_ablation" = 29, "regulatory_region_amplification" = 30,
    "regulatory_region_variant" = 31, "feature_elongation" = 32, "feature_truncation" = 33,
    "intergenic_variant" = 34
)

# Helpers (consequence severity) -----------------------------------------------

#' Compute the most-severe VEP term and its rank (uncached)
#'
#' Splits a VEP Consequence string on `&`, looks each term up in
#' `.gvr_effect_priority`, replaces unknown terms with `999L`, and returns
#' the lowest-priority term and its rank. PURE in the inputs.
#'
#' @param consequence A single VEP Consequence string
#'   (e.g. `"missense_variant&splice_region_variant"`). `NA` and empty
#'   inputs return `list(term = NA_character_, rank = 999L)`.
#'
#' @return A list with elements `term` (character) and `rank` (integer).
#'
#' @keywords internal
#' @noRd
.gvr_severe_term_uncached <- function(consequence) {
    if (is.na(consequence) || consequence == "") {
        return(list(term = NA_character_, rank = 999L))
    }
    terms <- strsplit(consequence, "&", fixed = TRUE)[[1]]
    pr <- .gvr_effect_priority[terms]
    pr[is.na(pr)] <- 999L
    wm <- which.min(pr)
    list(term = terms[wm], rank = as.integer(pr[wm]))
}

#' Compute the most-severe VEP term and its rank (memoised)
#'
#' Memoising wrapper around `.gvr_severe_term_uncached()`. The cache key is
#' the Consequence string itself (with a sentinel for `NA`). Output is
#' byte-identical to the uncached path by construction.
#'
#' @param consequence A single VEP Consequence string. See
#'   `.gvr_severe_term_uncached()`.
#' @param cache An `environment()` used as a key/value memoisation cache.
#'   Pass a fresh `new.env(parent = emptyenv())` per `read.gvr()` invocation
#'   to preserve Turn-3 cache scope (caches do NOT persist across calls).
#'
#' @return A list with elements `term` (character) and `rank` (integer).
#'
#' @keywords internal
#' @noRd
.gvr_severe_term <- function(consequence, cache) {
    key <- if (is.na(consequence)) "\001NA" else consequence
    hit <- cache[[key]]
    if (!is.null(hit)) return(hit)
    val <- .gvr_severe_term_uncached(consequence)
    assign(key, val, envir = cache)
    val
}

#' Most-severe VEP term (term only)
#'
#' Thin wrapper around `.gvr_severe_term()`. Returns only the term.
#'
#' @inheritParams .gvr_severe_term
#' @return Character of length 1: the most-severe Consequence SO term, or
#'   `NA_character_` for missing input.
#'
#' @keywords internal
#' @noRd
.gvr_most_severe_term <- function(consequence, cache) {
    .gvr_severe_term(consequence, cache = cache)$term
}

#' Most-severe VEP term rank (rank only)
#'
#' Thin wrapper around `.gvr_severe_term()`. Returns only the rank.
#'
#' @inheritParams .gvr_severe_term
#' @return Integer of length 1: the severity rank in `.gvr_effect_priority`
#'   (lower is more severe), or `999L` for missing input.
#'
#' @keywords internal
#' @noRd
.gvr_consequence_rank <- function(consequence, cache) {
    .gvr_severe_term(consequence, cache = cache)$rank
}

# Helpers (Variant_Classification dispatch) ------------------------------------

#' Map a VEP SO term + variant type to MAF Variant_Classification (uncached)
#'
#' PURE switch table from VEP SO consequence terms to the MAF
#' Variant_Classification vocabulary, with three small overrides for
#' frameshift/inframe events that depend on `var_type` (INS vs DEL). The
#' default branch returns `"Targeted_Region"` for terms not enumerated here.
#'
#' @param term The most-severe VEP SO term, e.g. `"missense_variant"`.
#' @param var_type Variant type as produced by `.gvr_coords()`
#'   (`"SNP"`, `"DEL"`, `"INS"`, `"DNP"`, `"TNP"`, `"ONP"`).
#'
#' @return A character of length 1: the MAF Variant_Classification.
#'
#' @keywords internal
#' @noRd
.gvr_vep_to_class_uncached <- function(term, var_type) {
    m <- switch(term,
        "splice_acceptor_variant" = "Splice_Site", "splice_donor_variant" = "Splice_Site",
        "transcript_ablation" = "Splice_Site", "exon_loss_variant" = "Splice_Site",
        "stop_gained" = "Nonsense_Mutation", "stop_lost" = "Nonstop_Mutation",
        "start_lost" = "Translation_Start_Site", "initiator_codon_variant" = "Translation_Start_Site",
        "missense_variant" = "Missense_Mutation", "conservative_missense_variant" = "Missense_Mutation",
        "rare_amino_acid_variant" = "Missense_Mutation", "protein_altering_variant" = "Missense_Mutation",
        "transcript_amplification" = "Intron", "splice_region_variant" = "Splice_Region",
        "splice_donor_5th_base_variant" = "Splice_Region", "splice_donor_region_variant" = "Splice_Region",
        "splice_polypyrimidine_tract_variant" = "Splice_Region", "stop_retained_variant" = "Silent",
        "synonymous_variant" = "Silent", "incomplete_terminal_codon_variant" = "Silent",
        "coding_sequence_variant" = "Missense_Mutation", "mature_miRNA_variant" = "RNA",
        "5_prime_UTR_variant" = "5'UTR", "5_prime_UTR_premature_start_codon_gain_variant" = "5'UTR",
        "3_prime_UTR_variant" = "3'UTR", "non_coding_transcript_exon_variant" = "RNA",
        "non_coding_exon_variant" = "RNA", "non_coding_transcript_variant" = "RNA",
        "nc_transcript_variant" = "RNA", "NMD_transcript_variant" = "Silent",
        "intron_variant" = "Intron", "upstream_gene_variant" = "5'Flank",
        "downstream_gene_variant" = "3'Flank", "TFBS_ablation" = "Targeted_Region",
        "TFBS_amplification" = "Targeted_Region", "TF_binding_site_variant" = "IGR",
        "regulatory_region_ablation" = "Targeted_Region", "regulatory_region_amplification" = "Targeted_Region",
        "regulatory_region_variant" = "IGR", "feature_elongation" = "Targeted_Region",
        "feature_truncation" = "Targeted_Region", "intergenic_variant" = "IGR",
        "Targeted_Region"
    )
    if (term == "frameshift_variant") m <- if (var_type == "DEL") "Frame_Shift_Del" else "Frame_Shift_Ins"
    if (term == "inframe_insertion") m <- "In_Frame_Ins"
    if (term == "inframe_deletion") m <- "In_Frame_Del"
    if (term == "disruptive_inframe_insertion") m <- "In_Frame_Ins"
    if (term == "disruptive_inframe_deletion") m <- "In_Frame_Del"
    m
}

#' Map a VEP SO term + variant type to MAF Variant_Classification (memoised)
#'
#' Memoising wrapper around `.gvr_vep_to_class_uncached()`. Key is the
#' composite `(term, var_type)` with sentinels for NA. Output is
#' byte-identical to the uncached path.
#'
#' @inheritParams .gvr_vep_to_class_uncached
#' @param cache An `environment()` used as a memoisation cache. Pass a fresh
#'   `new.env(parent = emptyenv())` per `read.gvr()` invocation.
#'
#' @return The MAF Variant_Classification string.
#'
#' @keywords internal
#' @noRd
.gvr_vep_to_class <- function(term, var_type, cache) {
    key <- paste0(
        if (is.na(term)) "\001NA" else term, "\002",
        if (is.na(var_type)) "\001NA" else var_type
    )
    hit <- cache[[key]]
    if (!is.null(hit)) return(hit)
    val <- .gvr_vep_to_class_uncached(term, var_type)
    assign(key, val, envir = cache)
    val
}
