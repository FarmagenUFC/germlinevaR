# =============================================================================
# read_gvr_vcf_utils.R
#
# Internal helpers used by read.gvr() to parse generic VCF fields:
#   * INFO column accessors (single-key and parse-once-look-up-many)
#   * Per-ALT genotype code extraction
#   * Indel-aware MAF coordinate / allele conversion
#
# All helpers are package-internal: not exported, no Rd man page generated
# (@keywords internal + @noRd). Promoted from read.gvr()'s body in the
# Turn 5a refactor; behaviour is byte-identical to the previous closure form.
# =============================================================================

#' Fetch a single INFO key's value
#'
#' Scans a VCF `INFO` string for `key=value` and returns the value, or
#' `NA_character_` if the key is missing or is present as a bare flag
#' (`;key;`).
#'
#' @param info_str A single character value: the raw `INFO` column for one
#'   record.
#' @param key The INFO field name to retrieve (single character).
#'
#' @return The value as character, or `NA_character_` when the key is not
#'   present as a `key=value` pair.
#'
#' @keywords internal
#' @noRd
.gvr_info_get <- function(info_str, key) {
    pat <- paste0("(^|;)", key, "=([^;]*)")
    m <- regmatches(info_str, regexec(pat, info_str))[[1]]
    if (length(m) >= 3) m[3] else NA_character_
}

#' Parse an INFO string once into a named character vector
#'
#' Equivalent to calling `.gvr_info_get()` for every key in `info_str`, but
#' executes a single split rather than one regex per key. Bare flags
#' (`;NMD;`) are excluded from the result, matching `.gvr_info_get()`'s
#' NA-on-missing-key behaviour.
#'
#' @param info_str A single character value: the raw `INFO` column for one
#'   record.
#'
#' @return A named character vector. Names are INFO keys; values are the
#'   corresponding `value` strings. Empty input or `NA` input returns a
#'   zero-length named character vector.
#'
#' @keywords internal
#' @noRd
.gvr_info_parse <- function(info_str) {
    if (is.na(info_str) || info_str == "") {
        return(setNames(character(0), character(0)))
    }
    kv <- strsplit(info_str, ";", fixed = TRUE)[[1]]
    eq <- regexpr("=", kv, fixed = TRUE)
    has <- eq > 0L
    keys <- substr(kv, 1L, eq - 1L)
    vals <- substr(kv, eq + 1L, nchar(kv))
    setNames(vals[has], keys[has])
}

#' Genotype-derived allele codes for a single split ALT row
#'
#' For a VCF row that may have multiple ALT alleles and was split so that
#' `Allele2` corresponds to a particular ALT index `ai`, derive the two
#' integer allele codes appropriate for that ALT row. Handles homozygous-ALT
#' (returns `(ai, ai)`), heterozygous with REF (returns `(0L, ai)`), and
#' compound-het cases (returns `(other_alt_idx, ai)`). NA-coercion warnings
#' from `as.integer()` on non-numeric tokens (e.g. `.`) are intentionally
#' suppressed; those positions are dropped from the genotype.
#'
#' @param gt The genotype string for this record (e.g. `"0/1"`, `"1|2"`,
#'   or `"./."`).
#' @param ai The ALT index (1-based) that this row corresponds to.
#'
#' @return A list with elements `c1` and `c2` (integer codes for the two
#'   alleles assigned to the split row).
#'
#' @keywords internal
#' @noRd
.gvr_gt_codes_for_alt <- function(gt, ai) {
    ai <- as.integer(ai)
    # why: as.integer() on a split genotype string; '.' and other non-numeric
    #      tokens become NA and are dropped by !is.na() below.
    gidx <- suppressWarnings(as.integer(strsplit(gt, "[/|]")[[1]]))
    gidx <- gidx[!is.na(gidx)]
    if (length(gidx) == 0L) return(list(c1 = NA_integer_, c2 = ai))
    others <- gidx[gidx != ai]
    if (length(others) == 0L) return(list(c1 = ai, c2 = ai))
    partner <- if (0L %in% others) 0L else others[1]
    list(c1 = partner, c2 = ai)
}

#' MAF-style coordinates and alleles for a single REF/ALT pair
#'
#' Converts a (POS, REF, ALT) triple into the MAF-style fields
#' (`Variant_Type`, `Start_Position`, `End_Position`, `Reference_Allele`,
#' `Tumor_Seq_Allele2`). Handles SNPs, ONPs (DNP/TNP/ONP), insertions,
#' deletions, and spanning deletions (`alt == "*"`).
#'
#' @param pos VCF POS (1-based start coordinate). Coerced to integer.
#' @param ref REF allele.
#' @param alt ALT allele.
#'
#' @return A list of length 5: `var_type`, `start`, `end`, `ref_allele`,
#'   `tum_allele2`. For insertions, `start` is set to `pos + rl - 1` and
#'   `ref_allele` is `"-"`; for deletions, `start` is `pos + al` and
#'   `tum_allele2` is `"-"`, mirroring the MAF convention.
#'
#' @keywords internal
#' @noRd
.gvr_coords <- function(pos, ref, alt) {
    pos <- as.integer(pos)
    if (alt == "*") {
        return(list(
            var_type = "DEL", start = pos, end = pos,
            ref_allele = ref, tum_allele2 = "*"
        ))
    }
    rl <- nchar(ref)
    al <- nchar(alt)
    if (rl == al && rl == 1) {
        return(list(
            var_type = "SNP", start = pos, end = pos,
            ref_allele = ref, tum_allele2 = alt
        ))
    }
    if (rl == al && rl > 1) {
        vt <- switch(as.character(rl), "2" = "DNP", "3" = "TNP", "ONP")
        return(list(
            var_type = vt, start = pos, end = pos + rl - 1,
            ref_allele = ref, tum_allele2 = alt
        ))
    }
    if (al > rl) {
        ins <- substr(alt, rl + 1, al)
        return(list(
            var_type = "INS", start = pos + rl - 1, end = pos + rl,
            ref_allele = "-", tum_allele2 = ins
        ))
    } else {
        del <- substr(ref, al + 1, rl)
        return(list(
            var_type = "DEL", start = pos + al, end = pos + rl - 1,
            ref_allele = del, tum_allele2 = "-"
        ))
    }
}
