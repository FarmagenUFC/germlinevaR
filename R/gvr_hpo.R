# gvr_hpo.R - HPO phenotype-to-gene resolution helper for germlinevaR
#
# Resolves one or more Human Phenotype Ontology (HPO) identifiers to a
# character vector of HGNC gene symbols by looking them up in the HPO
# phenotype-to-gene association table
# (https://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt).
#
# Public helper: gvr_hpo_genes()
# Consumed by:   read.gvr(), read.gvr.snpeff(), read.gvr.dual(), gvr_filter()
#                (each via an `hpo` argument, unioned with `genes` and `panel`).
#
# The downloaded table is cached under tools::R_user_dir("germlinevaR", "cache")
# and auto-refreshed after `max_age_days` days (default 30). Users on
# air-gapped systems can pass `hpo_path = <local file>` to bypass the network.
# ----------------------------------------------------------------------------


#' Resolve HPO phenotype terms to associated genes
#'
#' @description
#' Retrieves genes associated with one or more Human Phenotype Ontology (HPO)
#' terms using the HPO phenotype-to-gene association table. This is mainly used
#' by [read.gvr()] and [gvr_filter()] through their `hpo` argument, but can
#' also be called directly to inspect the gene set before filtering.
#'
#' The HPO table is downloaded once and cached under
#' `tools::R_user_dir("germlinevaR", "cache")`. Cached files older than
#' `max_age_days` days trigger an automatic re-download; setting
#' `refresh_cache = TRUE` forces refresh regardless of age. Use
#' `hpo_path = <local file>` for offline / air-gapped operation.
#'
#' Input HPO identifiers are lenient-parsed: `"HP:0003002"`, `"hp:0003002"`,
#' `"3002"`, and `"0003002"` are all accepted and normalised to canonical
#' `"HP:0003002"` before lookup. When normalisation happens a `message()` is
#' emitted so typos remain visible.
#'
#' No ontology-descendant expansion is performed in this release: only genes
#' associated with the exact HPO terms supplied are returned. To include
#' descendants, resolve them yourself first (e.g. via `ontologyIndex::get_descendants()`).
#'
#' @param hpo Character vector of HPO identifiers, e.g. `"HP:0003002"` or
#'   `c("HP:0003002", "HP:0001939")`. Lenient input accepted (see Description).
#' @param hpo_path Optional path to a local `phenotype_to_genes.txt` file.
#'   If supplied, no download is attempted.
#' @param hpo_url URL used when `hpo_path = NULL`. Default is the canonical
#'   HPO/OBO phenotype-to-gene association file.
#' @param cache_dir Directory used to cache the downloaded HPO file. `NULL`
#'   uses `tools::R_user_dir("germlinevaR", "cache")`.
#' @param refresh_cache Logical. If `TRUE`, force re-download of the HPO file
#'   ignoring `max_age_days`. Default `FALSE`.
#' @param max_age_days Numeric. Automatic refresh threshold in days for the
#'   cached HPO file. Default 30. Set to `Inf` to disable age-based refresh.
#' @param verbose Logical. If `TRUE`, print progress messages (download,
#'   cache age, per-code resolution counts, normalisation events). Default TRUE.
#'
#' @return A sorted, upper-cased, deduplicated character vector of HGNC gene
#'   symbols associated with the supplied HPO term(s). If nothing resolves,
#'   returns `character(0)` and emits a warning listing the unresolved terms.
#'
#' @examples
#' # Runnable example: use the tiny HPO fixture shipped in inst/extdata so
#' # the example needs no network access.
#' hpo_fx <- system.file("extdata", "hpo_phenotype_to_genes_mini.tsv",
#'                       package = "germlinevaR")
#' genes <- gvr_hpo_genes("HP:0003002", hpo_path = hpo_fx)
#' head(genes)
#'
#' # Lenient input forms all normalise to the canonical HP:0003002:
#' identical(genes,
#'           gvr_hpo_genes("3002", hpo_path = hpo_fx))
#'
#' # Multiple terms are unioned:
#' gvr_hpo_genes(c("HP:0003002", "HP:0025022"), hpo_path = hpo_fx)
#'
#' \donttest{
#' # Network use (downloads and caches the full HPO table on first call):
#' # genes <- gvr_hpo_genes("HP:0003002")
#'
#' # Force refresh of the cached HPO file:
#' # genes <- gvr_hpo_genes("HP:0003002", refresh_cache = TRUE)
#'
#' # Use with read.gvr() to restrict to HPO-implicated genes:
#' # gvr <- read.gvr("path/to/vcfs", hpo = "HP:0003002")
#'
#' # Or on an already-loaded table:
#' # gvr_flt <- gvr_filter(gvr, hpo = "HP:0003002")
#' }
#'
#' @family germlinevaR
#' @export
gvr_hpo_genes <- function(hpo,
                          hpo_path      = NULL,
                          hpo_url       = "https://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt",
                          cache_dir     = NULL,
                          refresh_cache = FALSE,
                          max_age_days  = 30,
                          verbose       = TRUE) {
    if (is.null(hpo) || length(hpo) == 0L) {
        return(character(0))
    }

    # --- Lenient input coercion -----------------------------------------------
    #   Accept "HP:0003002", "hp:0003002", "3002", "0003002" -> canonical
    #   "HP:0003002". Emit an info message when any coercion happens so typos
    #   remain visible to the user.
    raw <- unique(trimws(as.character(hpo)))
    raw <- raw[!is.na(raw) & nzchar(raw)]

    hpo_norm <- .gvr_hpo_normalise(raw, verbose = verbose)

    bad_hpo <- !grepl("^HP:[0-9]{7}$", hpo_norm)
    if (any(bad_hpo)) {
        stop(
            "gvr_hpo_genes: unrecognisable HPO identifier(s) after normalisation: ",
            paste(hpo_norm[bad_hpo], collapse = ", "),
            ". Expected an integer id (e.g. '3002'), zero-padded seven-digit id ",
            "('0003002'), or full 'HP:0003002' form.",
            call. = FALSE
        )
    }

    # --- Locate the HPO table (local / cached / downloaded) -------------------
    hpo_file <- .gvr_hpo_resolve_file(
        hpo_path      = hpo_path,
        hpo_url       = hpo_url,
        cache_dir     = cache_dir,
        refresh_cache = refresh_cache,
        max_age_days  = max_age_days,
        verbose       = verbose
    )

    hpo_dt <- .gvr_hpo_read_table(hpo_file)

    # --- Look up each requested term ------------------------------------------
    hit <- hpo_dt[hpo_dt[["HPO_ID"]] %in% hpo_norm, ]

    if (isTRUE(verbose)) {
        # Per-code resolution summary: how many genes did each input term hit?
        per_code_n <- vapply(hpo_norm, function(id) {
            sum(hpo_dt[["HPO_ID"]] == id)
        }, integer(1))
        message(sprintf(
            "gvr_hpo_genes: resolved %d HPO term(s) -> %d gene rows (per-code: %s)",
            length(hpo_norm),
            nrow(hit),
            paste(sprintf("%s=%d", hpo_norm, per_code_n), collapse = ", ")
        ))
    }

    unresolved <- hpo_norm[!hpo_norm %in% hpo_dt[["HPO_ID"]]]
    if (length(unresolved) > 0L) {
        warning(
            "gvr_hpo_genes: ", length(unresolved),
            " HPO term(s) not found in the table: ",
            paste(unresolved, collapse = ", "),
            call. = FALSE
        )
    }

    if (nrow(hit) == 0L) {
        return(character(0))
    }

    genes <- unique(trimws(hit$Gene_Symbol))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    sort(unique(toupper(genes)))
}


# =============================================================================
# Internal helpers (unexported)
# =============================================================================


# Normalise HPO input strings to canonical "HP:NNNNNNN" form.
# Accepts: "HP:0003002", "hp:0003002", "hp:3002", "3002", "0003002".
# Emits a message() at verbose=TRUE listing any inputs that were coerced.
.gvr_hpo_normalise <- function(x, verbose = TRUE) {
    x <- trimws(as.character(x))
    normd <- vapply(x, function(s) {
        if (grepl("^HP:[0-9]{7}$", s)) return(s)         # already canonical
        if (grepl("^hp:[0-9]+$", s, ignore.case = TRUE)) {
            num <- sub("^hp:", "", s, ignore.case = TRUE)
            if (nchar(num) <= 7L) {
                return(sprintf("HP:%07d", as.integer(num)))
            }
        }
        if (grepl("^[0-9]+$", s) && nchar(s) <= 7L) {
            return(sprintf("HP:%07d", as.integer(s)))
        }
        s   # unrecognised; downstream validator will error
    }, character(1))

    coerced <- x[x != normd]
    if (length(coerced) > 0L && isTRUE(verbose)) {
        message(sprintf(
            "gvr_hpo_genes: normalised %d HPO input(s) to canonical form: %s",
            length(coerced),
            paste(sprintf("'%s' -> '%s'", coerced, normd[x != normd]),
                  collapse = ", ")
        ))
    }
    unname(normd)
}


# Resolve the on-disk path to the HPO phenotype_to_genes.txt file:
# - if `hpo_path` supplied and exists: return it verbatim
# - else: use / create a cache under `cache_dir` (default R_user_dir)
# - if cache file exists and is younger than `max_age_days` and !refresh_cache:
#     reuse
# - else: download from `hpo_url` (method = "libcurl" for reliable redirect
#   handling on Windows) and validate size > 1 KB.
.gvr_hpo_resolve_file <- function(hpo_path      = NULL,
                                  hpo_url,
                                  cache_dir     = NULL,
                                  refresh_cache = FALSE,
                                  max_age_days  = 30,
                                  verbose       = TRUE) {
    if (!is.null(hpo_path)) {
        if (!file.exists(hpo_path)) {
            stop("gvr_hpo_genes: local HPO file not found: ", hpo_path,
                 call. = FALSE)
        }
        return(normalizePath(hpo_path, mustWork = TRUE))
    }

    if (is.null(cache_dir)) {
        cache_dir <- tools::R_user_dir("germlinevaR", "cache")
    }

    if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    }

    dest <- file.path(cache_dir, "hpo_phenotype_to_genes.txt")

    # Age-based auto-refresh: if cached file exists but is older than
    # `max_age_days`, re-download unless the user opted out (Inf) or forced
    # refresh (TRUE, handled below).
    stale <- FALSE
    if (file.exists(dest) && is.finite(max_age_days)) {
        age_days <- as.numeric(difftime(
            Sys.time(), file.info(dest)$mtime, units = "days"
        ))
        if (isTRUE(age_days > max_age_days)) {
            stale <- TRUE
            if (isTRUE(verbose)) {
                message(sprintf(
                    "gvr_hpo_genes: cached HPO file is %.1f days old (> %.0f); refreshing",
                    age_days, max_age_days
                ))
            }
        }
    }

    if (!file.exists(dest) || isTRUE(refresh_cache) || isTRUE(stale)) {
        if (isTRUE(verbose)) {
            message("gvr_hpo_genes: downloading HPO phenotype-to-gene table from ",
                    hpo_url)
        }

        # Download to a scratch file so a partial or failed transfer never
        # clobbers a previously-good cached copy. On success we atomic-rename
        # the scratch file into place.
        scratch <- paste0(dest, ".part")
        if (file.exists(scratch)) unlink(scratch)

        # utils::download.file() can emit a warning (not an error) for DNS
        # failures, 4xx/5xx responses, and connection resets. We escalate
        # those warnings to errors for the duration of this call so tryCatch
        # can react uniformly.
        dl_error <- tryCatch(
            {
                withCallingHandlers(
                    utils::download.file(
                        url      = hpo_url,
                        destfile = scratch,
                        mode     = "wb",
                        method   = "libcurl",   # reliable redirect handling
                        quiet    = !isTRUE(verbose)
                    ),
                    warning = function(w) {
                        stop(conditionMessage(w), call. = FALSE)
                    }
                )
                NULL
            },
            error = function(e) conditionMessage(e)
        )

        # Post-download sanity: real HPO file is >1 MB. Any download that
        # produced <1 KB is treated as failure regardless of what the tool
        # reported.
        bad_download <-
            !is.null(dl_error) ||
            !file.exists(scratch) ||
            isTRUE(is.na(file.size(scratch))) ||
            isTRUE(file.size(scratch) < 1024L)

        if (isTRUE(bad_download)) {
            if (file.exists(scratch)) unlink(scratch)
            err_chunk <- if (!is.null(dl_error))
                sprintf(" (download error: %s)", dl_error) else ""
            if (file.exists(dest)) {
                warning(sprintf(
                    "could not refresh HPO table from %s; reusing existing cached copy at %s%s",
                    hpo_url, dest, err_chunk
                ), call. = FALSE)
            } else {
                stop(sprintf(
                    "could not download HPO table from %s. Use 'hpo_path=' with a local phenotype_to_genes.txt file%s.",
                    hpo_url, err_chunk
                ), call. = FALSE)
            }
        } else {
            # Atomic promote: replace `dest` with the fully-downloaded scratch
            # file. file.rename is atomic on POSIX and best-effort on Windows.
            if (file.exists(dest)) unlink(dest)
            ok <- file.rename(scratch, dest)
            if (!isTRUE(ok)) {
                # Rare cross-device rename failure; fall back to copy + unlink
                file.copy(scratch, dest, overwrite = TRUE)
                unlink(scratch)
            }
        }
    }

    dest
}


# Read the HPO phenotype_to_genes.txt file and return a two-column
# data.table (HPO_ID, Gene_Symbol). Column-name autodetect first, positional
# fallback second (HPO always has hpo_id in col 1 and gene_symbol in col 4).
.gvr_hpo_read_table <- function(path) {
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("gvr_hpo_genes requires the 'data.table' package.", call. = FALSE)
    }

    # HPO's phenotype_to_genes.txt starts with a "#format:..." comment line
    # followed by the header row. data.table::fread has no `comment.char`, so
    # pre-scan the file to find the first line that does NOT start with '#'
    # and pass that 1-based line number to fread as `skip=`. This is robust
    # to header renames upstream (whether the header is `hpo_id` or `HPO-ID`
    # or anything else).
    skip_lines <- 0L
    con <- file(path, open = "r", encoding = "UTF-8")
    on.exit(close(con), add = TRUE)
    while (TRUE) {
        ln <- readLines(con, n = 1L, warn = FALSE)
        if (length(ln) == 0L) break
        if (!grepl("^\\s*#", ln)) break
        skip_lines <- skip_lines + 1L
    }
    close(con)
    on.exit()

    dt <- data.table::fread(
        path,
        sep          = "\t",
        header       = TRUE,
        skip         = skip_lines,
        fill         = TRUE,
        quote        = "",
        na.strings   = "",   # only empty cells are NA; literal "NA" stays a string
        data.table   = TRUE
    )

    nm <- names(dt)

    hpo_col  <- intersect(nm, c("HPO-ID", "HPO_ID", "hpo_id", "HPOId", "HPO_id"))
    gene_col <- intersect(nm, c("Gene-Name", "Gene_Name", "gene_symbol",
                                "GeneSymbol", "gene_name", "Gene_Symbol"))

    if (length(hpo_col) == 0L || length(gene_col) == 0L) {
        # Positional fallback: current HPO schema is
        #   hpo_id  hpo_name  ncbi_gene_id  gene_symbol  disease_id
        # so col 1 = HPO id, col 4 = gene symbol.
        if (ncol(dt) >= 4L) {
            hpo_col  <- nm[1L]
            gene_col <- nm[4L]
        } else {
            stop(
                "gvr_hpo_genes: could not identify HPO and gene-symbol columns ",
                "in HPO table (columns found: ",
                paste(nm, collapse = ", "), ").",
                call. = FALSE
            )
        }
    }

    out <- data.table::data.table(
        HPO_ID      = as.character(dt[[hpo_col[1L]]]),
        Gene_Symbol = as.character(dt[[gene_col[1L]]])
    )

    out <- out[!is.na(out[["HPO_ID"]])      & nzchar(out[["HPO_ID"]]), ]
    out <- out[!is.na(out[["Gene_Symbol"]]) & nzchar(out[["Gene_Symbol"]]), ]

    out
}
