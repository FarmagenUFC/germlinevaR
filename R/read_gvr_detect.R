# =============================================================================
# read_gvr_detect.R
#
# Internal helpers used by read.gvr() to detect:
#   * the location of sibling .R files in standalone (non-package) mode
#       (.gvr_locate_sibling)
#   * which annotator wrote a VCF (VEP, SnpEff, or both)
#       (.gvr_detect_annotator)
#   * which reference build a VCF is from (GRCh37/GRCh38/T2T-CHM13v2.0)
#       (.gvr_detect_build)
#
# `.gvr_locate_sibling()` is used by read.gvr() in standalone mode (when the
# package is source()'d rather than loaded as a library) to auto-source its
# sibling .R files. In a package context this branch is skipped because R
# already loads every R/*.R file at namespace creation.
#
# `.gvr_detect_annotator()` and `.gvr_detect_build()` are PURE; the sibling
# read.gvr.snpeff.R file still nests its own byte-identical copy of
# `.detect_annotator()` (so its standalone-mode dispatch works without
# depending on this file). This is acceptable code duplication scoped to the
# legacy standalone path; the package code path uses only `.gvr_detect_*`.
#
# All helpers are package-internal: not exported, no Rd page (@keywords
# internal + @noRd). Promoted from read.gvr()'s body in the Turn 5a refactor.
# =============================================================================

#' Locate a sibling .R file next to this script (standalone-mode helper)
#'
#' Used by `read.gvr()` in standalone mode to discover where its sibling
#' files (`read.gvr.snpeff.R`, `read.gvr.dual.R`) live so they can be
#' `source()`d. Tries the directory containing this script first (via
#' `utils::getSrcFilename()`), then falls back to the current working
#' directory.
#'
#' @param basename_ The bare filename to look for (e.g.
#'   `"read.gvr.snpeff.R"`).
#'
#' @return Full path to the sibling file if found, or `NULL` if not present
#'   in either candidate directory.
#'
#' @keywords internal
#' @noRd
.gvr_locate_sibling <- function(basename_) {
    this_dir <- tryCatch(
        {
            fn <- function() {}
            fp <- utils::getSrcFilename(fn, full.names = TRUE)
            if (length(fp) && nzchar(fp)) {
                normalizePath(dirname(fp), mustWork = FALSE)
            } else {
                NULL
            }
        },
        error = function(e) NULL
    )
    if (is.null(this_dir) || !nzchar(this_dir)) this_dir <- getwd()
    candidate <- file.path(this_dir, basename_)
    if (file.exists(candidate)) return(candidate)
    NULL
}

#' Detect the variant annotator(s) used to produce a VCF
#'
#' Reads the header of a gzipped VCF looking for `##INFO=<ID=...>` lines,
#' and reports which of VEP (`CSQ`), SnpEff (`ANN`), or both are present.
#' Returns the annotator label plus the full list of INFO tags as an
#' attribute (used by callers for diagnostics).
#'
#' @param path Path to a `.vcf.gz` file.
#'
#' @return One of `"vep"`, `"snpeff"`, `"dual"`, or `NA_character_` when
#'   neither annotator is detected. The full list of `##INFO` tags is
#'   attached as `attr(res, "info_tags")` for downstream diagnostics.
#'
#' @keywords internal
#' @noRd
.gvr_detect_annotator <- function(path) {
    con <- gzfile(path, "r")
    on.exit(close(con))
    csq_found <- FALSE
    ann_found <- FALSE
    info_tags <- character(0L)
    repeat {
        line <- readLines(con, n = 1L, warn = FALSE)
        if (length(line) == 0L) break
        if (!startsWith(line, "##")) break
        if (startsWith(line, "##INFO=<ID=")) {
            tag <- sub("^##INFO=<ID=([^,>]+).*$", "\\1", line)
            info_tags <- c(info_tags, tag)
            if (identical(tag, "CSQ")) csq_found <- TRUE
            if (identical(tag, "ANN")) ann_found <- TRUE
        }
    }
    if (csq_found && ann_found) {
        res <- "dual"
    } else if (csq_found) {
        res <- "vep"
    } else if (ann_found) {
        res <- "snpeff"
    } else {
        res <- NA_character_
    }
    attr(res, "info_tags") <- info_tags
    res
}

#' Detect the reference build of a VCF from its header
#'
#' Combines three independent signals from a VCF header to infer the
#' reference build:
#'
#' * VEP signal: `##VEP="v..." ... assembly="GRCh3X.p..."` attribute.
#' * SnpEff signal: `##SnpEffCmd="SnpEff ... GRCh3X.NN input.vcf"`
#'   database token (also recognises legacy `hg19`/`hg38`).
#' * Contig signal: the first `##contig=<ID=(chr)?(1|M),length=N>` whose
#'   length matches a known canonical contig (GRCh37/GRCh38/T2T-CHM13v2.0).
#'
#' Confidence is `"high"` when at least two signals agree on the same
#' canonical label, `"low"` when only one signal is present, and `"none"`
#' otherwise (including the case where 2-3 signals are present but no two
#' agree).
#'
#' @param path Path to a `.vcf.gz` file.
#'
#' @return A list with elements:
#'   * `label` - one of `"GRCh38"`, `"GRCh37"`, `"T2T-CHM13v2.0"`, or
#'       `NA_character_`.
#'   * `confidence` - `"high"`, `"low"`, or `"none"`.
#'   * `signals` - a list of length 3 (`vep`, `snpeff`, `contig`), each
#'       either a canonical label or `NA_character_`.
#'
#' @keywords internal
#' @noRd
.gvr_detect_build <- function(path) {
    contig_lookup <- list(
        `1` = list(
            `248956422` = "GRCh38",
            `249250621` = "GRCh37",
            `248387328` = "T2T-CHM13v2.0"
        ),
        `M` = list(`16569` = "GRCh38", `16571` = "GRCh37")
    )

    con <- gzfile(path, "r")
    on.exit(close(con))
    vep_sig <- NA_character_
    snpeff_sig <- NA_character_
    contig_sig <- NA_character_
    first_std_contig_seen <- FALSE
    repeat {
        line <- readLines(con, n = 1L, warn = FALSE)
        if (length(line) == 0L) break
        if (!startsWith(line, "##")) break

        if (is.na(vep_sig) && startsWith(line, "##VEP=")) {
            m <- regmatches(line, regexec("assembly=\"(GRCh3[78])", line))[[1]]
            if (length(m) == 2L) vep_sig <- m[2]
        }

        if (is.na(snpeff_sig) && startsWith(line, "##SnpEffCmd=")) {
            m <- regmatches(line, regexec("\\b(GRCh3[78])\\.[0-9]+\\b", line))[[1]]
            if (length(m) == 2L) {
                snpeff_sig <- m[2]
            } else {
                m2 <- regmatches(line, regexec("\\b(hg19|hg38)\\b", line))[[1]]
                if (length(m2) == 2L) {
                    snpeff_sig <- c(hg19 = "GRCh37", hg38 = "GRCh38")[m2[2]]
                }
            }
        }

        if (!first_std_contig_seen && startsWith(line, "##contig=<ID=")) {
            m <- regmatches(line, regexec(
                "^##contig=<ID=(?:chr)?([0-9XYM]+),length=([0-9]+)", line
            ))[[1]]
            if (length(m) == 3L) {
                first_std_contig_seen <- TRUE
                chrom_key <- m[2]
                len_key <- m[3]
                if (!is.null(contig_lookup[[chrom_key]]) &&
                    !is.null(contig_lookup[[chrom_key]][[len_key]])) {
                    contig_sig <- contig_lookup[[chrom_key]][[len_key]]
                }
            }
        }
    }

    sigs <- c(vep = vep_sig, snpeff = snpeff_sig, contig = contig_sig)
    sigs_present <- sigs[!is.na(sigs)]
    if (length(sigs_present) == 0L) {
        return(list(
            label = NA_character_, confidence = "none",
            signals = as.list(sigs)
        ))
    }
    tbl <- sort(table(sigs_present), decreasing = TRUE)
    top <- names(tbl)[1L]
    top_count <- as.integer(tbl[1L])
    if (top_count >= 2L) {
        return(list(label = top, confidence = "high", signals = as.list(sigs)))
    }
    if (length(sigs_present) == 1L) {
        return(list(label = top, confidence = "low", signals = as.list(sigs)))
    }
    list(label = NA_character_, confidence = "none", signals = as.list(sigs))
}
