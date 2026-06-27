#' Clear the auto-fetched protein-domain cache used by `gvr_lollipop`
#'
#' @description
#' Deletes cached InterPro domain `.rds` files written by
#' `gvr_lollipop(..., domains = "auto")`. Useful when InterPro releases a new
#' version (roughly quarterly) and you want fresh annotations on the next
#' call. Also clears the matching in-memory session cache so the next call
#' actually re-fetches instead of returning the in-memory copy.
#'
#' @details
#' The cache directory is resolved with the same precedence chain used by
#' `gvr_lollipop()`, but here the chain only **finds** an existing directory;
#' it never creates one. The precedence (first hit wins) is:
#' \enumerate{
#'   \item `cache_dir` argument (explicit override).
#'   \item Environment variable `GVR_CACHE_DIR`.
#'   \item R option `getOption("germlinevaR.cache_dir")`.
#'   \item `tools::R_user_dir("germlinevaR", "cache")`.
#'   \item `file.path(tempdir(), "germlinevaR_cache")`.
#' }
#' If the resolved directory does not exist, the helper returns invisibly
#' with a `verbose` message and no files deleted (this is not an error;
#' a missing cache directory simply means there was nothing to clear).
#'
#' @param gene Character(1) or `NULL`. Gene symbol to delete from the cache.
#'   `NULL` (default) deletes every cached domain file in the cache directory
#'   (i.e., every file matching `^domains_interpro_.*\\.rds$`).
#' @param organism Integer or character or `NULL`. NCBI taxonomy id used as
#'   part of the cache key. When `gene` is given but `organism` is `NULL`,
#'   every taxonomy variant of that gene is deleted (i.e., the pattern is
#'   `^domains_interpro_<GENE>_.*\\.rds$`). Default `NULL`.
#' @param cache_dir Character(1) or `NULL`. Override the cache directory.
#'   `NULL` (default) triggers the precedence chain above.
#' @param verbose Logical(1). Print one line per deleted file. Default
#'   `TRUE`.
#'
#' @return Invisibly, a character vector of the file paths that were
#'   deleted (length 0 if none).
#'
#' @seealso [gvr_lollipop()]
#'
#' @examples
#' ## Clear the cache for a specific gene (safe to run; no-op if cache is empty)
#' gvr_domain_cache_clear(gene = "TP53")
#'
#' ## Clear everything (all genes, all organisms)
#' gvr_domain_cache_clear()
#'
#' ## Clear only TP53 across all organisms
#' gvr_domain_cache_clear(gene = "TP53")
#'
#' ## Clear only human TP53
#' gvr_domain_cache_clear(gene = "TP53", organism = 9606)
#' @export
gvr_domain_cache_clear <- function(gene      = NULL,
                                   organism  = NULL,
                                   cache_dir = NULL,
                                   verbose   = TRUE) {
    # ---- Resolve cache directory (READ-ONLY: never create) ----
    candidates <- character(0)
    if (!is.null(cache_dir) && is.character(cache_dir) && length(cache_dir) == 1L &&
        nzchar(cache_dir)) {
        candidates <- c(candidates, cache_dir)
    }
    env_dir <- Sys.getenv("GVR_CACHE_DIR", unset = "")
    if (nzchar(env_dir)) candidates <- c(candidates, env_dir)
    opt_dir <- getOption("germlinevaR.cache_dir", default = NULL)
    if (!is.null(opt_dir) && is.character(opt_dir) && length(opt_dir) == 1L &&
        nzchar(opt_dir)) {
        candidates <- c(candidates, opt_dir)
    }
    candidates <- c(candidates, tools::R_user_dir("germlinevaR", which = "cache"))
    candidates <- c(candidates, file.path(tempdir(), "germlinevaR_cache"))

    cdir <- NA_character_
    for (cand in candidates) {
        if (dir.exists(cand)) {
            cdir <- cand
            break
        }
    }
    if (is.na(cdir)) {
        if (isTRUE(verbose))
            message("gvr_domain_cache_clear: no cache directory found (nothing to clear).")
        return(invisible(character(0)))
    }

    # ---- Build glob pattern based on gene/organism args ----
    gene_part <- if (is.null(gene)) ".*" else gsub("([.+*?\\[\\^\\]$()=!<>|:\\-#])",
        "\\\\\\1", gene)
    org_part  <- if (is.null(organism)) ".*" else as.character(organism)
    pattern   <- sprintf("^domains_interpro_%s_%s\\.rds$", gene_part, org_part)

    files <- list.files(cdir, pattern = pattern, full.names = TRUE)
    if (length(files) == 0L) {
        if (isTRUE(verbose))
            message(sprintf("gvr_domain_cache_clear: no matching cache files in %s", cdir))
    }

    # ---- Delete on-disk files ----
    deleted <- character(0)
    for (f in files) {
        ok <- tryCatch(file.remove(f), error = function(e) FALSE,
            warning = function(w) FALSE)
        if (isTRUE(ok)) {
            deleted <- c(deleted, f)
            if (isTRUE(verbose)) message(sprintf("gvr_domain_cache_clear: removed %s", f))
        } else {
            warning(sprintf("gvr_domain_cache_clear: could not remove %s", f), call. = FALSE)
        }
    }

    # ---- Clear matching in-memory entries ----
    # The in-memory cache lives at file scope in gvr_lollipop.R; rely on the
    # symbol being available in the calling environment chain (the same way
    # gvr_lollipop accesses it).
    if (exists(".gvr_domain_mem_cache", inherits = TRUE)) {
        mem <- get(".gvr_domain_mem_cache", inherits = TRUE)
        if (is.environment(mem)) {
            keys <- ls(mem, all.names = TRUE)
            drop <- if (is.null(gene)) {
                keys
            } else {
                # keys look like "GENE|ORG"
                parts <- strsplit(keys, "|", fixed = TRUE)
                keep_g <- vapply(parts, function(p) length(p) >= 1L && p[1] == gene, logical(1))
                keep_o <- if (is.null(organism)) {
                    rep(TRUE, length(parts))
                } else {
                    vapply(parts, function(p) length(p) >= 2L && p[2] == as.character(organism),
                        logical(1))
                }
                keys[keep_g & keep_o]
            }
            if (length(drop) > 0L) {
                rm(list = drop, envir = mem)
                if (isTRUE(verbose))
                    message(sprintf("gvr_domain_cache_clear: cleared %d in-memory entr%s",
                        length(drop), if (length(drop) == 1L) "y" else "ies"))
            }
        }
    }

    invisible(deleted)
}
