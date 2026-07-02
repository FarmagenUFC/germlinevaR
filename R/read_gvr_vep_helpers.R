# =============================================================================
# read_gvr_vep_helpers.R
#
# Internal helpers used by read.gvr() to parse VEP-specific fields in a VCF:
#   * URL decoding of VEP CSQ values
#   * HGVSp 3-letter -> 1-letter conversion (aa3to1 table + make_hgvsp_short)
#   * Candidate VEP CSQ-allele encodings for REF/ALT matching
#   * CSQ header parser and sample-name parser
#
# Promoted from read.gvr()'s body (Turn 5a refactor) to reduce the function's
# size without changing its behaviour. All helpers are package-internal:
# they are NOT exported and have NO Rd man page (@keywords internal + @noRd).
# Names use the `.gvr_` prefix to match the existing convention in
# gvr_lollipop.R (.gvr_http_get_retry, .gvr_abbrev_domain, .gvr_wrap_two_lines).
# =============================================================================

# Package-private constants -----------------------------------------------------

# Amino-acid 3-letter -> 1-letter map (27 codes incl. Ter/Sec/Pyl/ambiguity).
# Used by `.gvr_make_hgvsp_short()` to convert e.g. "p.Arg175His" to "p.R175H".
.gvr_aa3to1 <- c(
    Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C", Gln = "Q", Glu = "E",
    Gly = "G", His = "H", Ile = "I", Leu = "L", Lys = "K", Met = "M", Phe = "F",
    Pro = "P", Ser = "S", Thr = "T", Trp = "W", Tyr = "Y", Val = "V", Ter = "*",
    Sec = "U", Pyl = "O", Asx = "B", Glx = "Z", Xaa = "X", Xle = "J"
)

# Precompiled regex alternation for a single-pass 3->1 substitution. Built once
# at package load. All codes are literal alpha (regex-safe), and the 27 codes
# are mutually non-overlapping, so a single regex pass is order-independent and
# byte-identical to the original 27 sequential gsub() calls (verified by the
# upstream author over 28,686 unique real HGVSp strings, 0 diffs).
.gvr_aa3to1_pat <- paste(names(.gvr_aa3to1), collapse = "|")

# Helpers -----------------------------------------------------------------------

#' URL-decode a single VEP CSQ field value
#'
#' VEP's CSQ block percent-encodes a small set of reserved characters in
#' values (`=`, `;`, `,`, `:`). This helper reverses that encoding for a
#' single scalar character value. NA and empty-string inputs are passed
#' through unchanged so the helper can be applied uniformly to CSQ columns.
#'
#' @param x A single character value, possibly `NA`.
#'
#' @return A character of length 1 with `%3D`, `%3B`, `%2C`, `%3A`
#'   decoded back to `=`, `;`, `,`, `:`. `NA_character_` is returned for
#'   `NA` input; the empty string is returned for `""` input.
#'
#' @keywords internal
#' @noRd
.gvr_url_decode <- function(x) {
    if (is.na(x) || x == "") return(x)
    x <- gsub("%3D", "=", x, fixed = TRUE)
    x <- gsub("%3B", ";", x, fixed = TRUE)
    x <- gsub("%2C", ",", x, fixed = TRUE)
    x <- gsub("%3A", ":", x, fixed = TRUE)
    x
}

#' Build HGVSp_Short from a VEP HGVSp string
#'
#' Converts a VEP "HGVSp" string into the short `p.<oneletter>` form: drops
#' any transcript prefix (`ENSP00012345.6:p.Arg175His` -> `p.Arg175His`),
#' strips parentheses, then replaces every 3-letter amino-acid code with its
#' 1-letter equivalent in a single regex pass.
#'
#' @param hgvsp A single character value, possibly `NA`. Empty string is
#'   accepted and returned unchanged as the empty string.
#'
#' @return A character of length 1: `""` for empty/NA input, otherwise the
#'   transcript-stripped, parenthesis-stripped, 1-letter HGVSp. Byte-identical
#'   to the previous closure-form implementation in `read.gvr()`.
#'
#' @keywords internal
#' @noRd
.gvr_make_hgvsp_short <- function(hgvsp) {
    if (is.na(hgvsp) || hgvsp == "") return("")
    s <- .gvr_url_decode(hgvsp)
    s <- sub("^[^:]*:", "", s)
    s <- gsub("[()]", "", s)
    m <- gregexpr(.gvr_aa3to1_pat, s, perl = TRUE)
    regmatches(s, m) <- lapply(
        regmatches(s, m),
        function(hit) unname(.gvr_aa3to1[hit])
    )
    s
}

#' Candidate VEP CSQ-allele encodings for a REF/ALT pair
#'
#' Given a REF and a single ALT, returns the set of candidate strings that
#' VEP could have placed in the CSQ Allele column. Three strategies are
#' tried in priority order: (a) anchor-trim (drop one leading base);
#' (b) longest-common-prefix trim; (c) repeat-aware bidirectional trim
#' (common suffix, then common prefix); plus (d) the raw ALT as a fallback.
#' Duplicates are removed but order is preserved. `"-"` is used for empty
#' insertions/deletions and `"*"` is preserved for spanning deletions.
#'
#' @param ref REF allele (single character of length 1, may be multi-base).
#' @param alt ALT allele (single character of length 1, may be multi-base
#'   or `"*"`).
#'
#' @return A character vector of unique candidate allele strings, in the
#'   priority order described above.
#'
#' @keywords internal
#' @noRd
.gvr_vep_allele_candidates <- function(ref, alt) {
    if (alt == "*") return("*")
    if (nchar(ref) == 1 && nchar(alt) == 1) return(alt)
    rl <- nchar(ref)
    al <- nchar(alt)
    cand <- character(0)
    r1 <- substr(ref, 2, rl)
    a1 <- substr(alt, 2, al)
    cand <- c(cand, if (nchar(a1) == 0) "-" else if (nchar(r1) == 0) a1 else a1)
    k <- 0L
    mx <- min(rl, al)
    while (k < mx && substr(ref, k + 1, k + 1) == substr(alt, k + 1, k + 1)) k <- k + 1L
    rp <- substr(ref, k + 1, rl)
    ap <- substr(alt, k + 1, al)
    cand <- c(cand, if (nchar(ap) == 0) "-" else if (nchar(rp) == 0) ap else ap)
    rr <- ref
    aa <- alt
    while (nchar(rr) > 0 && nchar(aa) > 0 &&
        substr(rr, nchar(rr), nchar(rr)) == substr(aa, nchar(aa), nchar(aa))) {
        rr <- substr(rr, 1, nchar(rr) - 1)
        aa <- substr(aa, 1, nchar(aa) - 1)
    }
    j <- 0L
    mx2 <- min(nchar(rr), nchar(aa))
    while (j < mx2 && substr(rr, j + 1, j + 1) == substr(aa, j + 1, j + 1)) j <- j + 1L
    rb <- substr(rr, j + 1, nchar(rr))
    ab <- substr(aa, j + 1, nchar(aa))
    cand <- c(cand, if (nchar(ab) == 0) "-" else if (nchar(rb) == 0) ab else ab)
    cand <- c(cand, alt)
    unique(cand[cand != ""])
}

#' Parse the VEP CSQ "Format:" declaration from a VCF header
#'
#' Scans the header of a gzipped VCF for the `##INFO=<ID=CSQ,...,
#' Description="Consequence ... Format: A|B|C|...">` line and returns the
#' ordered field names as a character vector. Used by `read.gvr()` to map
#' each pipe-delimited CSQ block to named slots.
#'
#' @param path Path to a `.vcf.gz` file.
#'
#' @return Character vector of CSQ field names, in declaration order.
#'   Throws an error if the CSQ header is not found before `#CHROM`.
#'
#' @keywords internal
#' @noRd
.gvr_get_csq_fields <- function(path) {
    con <- gzfile(path, "r")
    on.exit(close(con))
    repeat {
        line <- readLines(con, n = 1L)
        if (length(line) == 0L) stop("CSQ header not found")
        if (startsWith(line, "##INFO=<ID=CSQ")) {
            fmt <- sub(".*Format: ", "", line)
            fmt <- sub('">$', "", fmt)
            return(strsplit(fmt, "|", fixed = TRUE)[[1]])
        }
        if (startsWith(line, "#CHROM")) stop("Reached #CHROM before CSQ header")
    }
}

#' Extract the sample name from the `#CHROM` header line of a VCF
#'
#' Reads the gzipped VCF header line by line, finds the `#CHROM ...` line,
#' tab-splits it, and returns the last column (the sample-name column for
#' single-sample VCFs). Throws an error if the `#CHROM` line is not found.
#'
#' @param path Path to a `.vcf.gz` file.
#'
#' @return A single character value: the sample name as written in the
#'   `#CHROM` line.
#'
#' @keywords internal
#' @noRd
.gvr_get_sample_name <- function(path) {
    con <- gzfile(path, "r")
    on.exit(close(con))
    repeat {
        line <- readLines(con, n = 1L)
        if (length(line) == 0L) stop("#CHROM not found")
        if (startsWith(line, "#CHROM")) {
            cols <- strsplit(line, "\t", fixed = TRUE)[[1]]
            return(cols[length(cols)])
        }
    }
}
