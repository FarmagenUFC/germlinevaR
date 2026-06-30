#' Resolve HPO phenotype terms to associated genes
#'
#' @description
#' Retrieves genes associated with one or more Human Phenotype Ontology (HPO)
#' terms using the HPO phenotype-to-gene association table. This is mainly used
#' by [gvr_filter()] through its `hpo` argument, but can also be called directly
#' to inspect the gene set before filtering.
#'
#' @param hpo Character vector of HPO identifiers, e.g. `"HP:0003002"`.
#' @param hpo_path Optional path to a local `phenotype_to_genes.txt` file.
#'   If supplied, no download is attempted.
#' @param hpo_url URL used when `hpo_path = NULL`. Default uses the HPO/OBO
#'   phenotype-to-gene association file.
#' @param cache_dir Directory used to cache the downloaded HPO file. `NULL`
#'   uses `tools::R_user_dir("germlinevaR", "cache")`.
#' @param refresh_cache Logical. If `TRUE`, force re-download of the HPO file.
#'   Default `FALSE`.
#' @param verbose Logical. If `TRUE`, print download/cache messages.
#'
#' @return A character vector of unique HGNC gene symbols associated with the
#'   supplied HPO term(s).
#'
#' @examples
#' \dontrun{
#' genes <- gvr_hpo_genes("HP:0003002")
#' genes
#' }
#'
#' @family germlinevaR
#' @export
gvr_hpo_genes <- function(hpo,
                          hpo_path = NULL,
                          hpo_url = "https://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt",
                          cache_dir = NULL,
                          refresh_cache = FALSE,
                          verbose = TRUE) {
    if (is.null(hpo) || length(hpo) == 0L) {
        return(character(0))
    }

    hpo <- unique(trimws(as.character(hpo)))
    hpo <- hpo[nzchar(hpo)]

    bad_hpo <- !grepl("^HP:[0-9]{7}$", hpo)
    if (any(bad_hpo)) {
        stop(
            "gvr_hpo_genes: invalid HPO identifier(s): ",
            paste(hpo[bad_hpo], collapse = ", "),
            ". Expected format like 'HP:0003002'.",
            call. = FALSE
        )
    }

    hpo_file <- .gvr_hpo_resolve_file(
        hpo_path = hpo_path,
        hpo_url = hpo_url,
        cache_dir = cache_dir,
        refresh_cache = refresh_cache,
        verbose = verbose
    )

    hpo_dt <- .gvr_hpo_read_table(hpo_file)

    hit <- hpo_dt[HPO_ID %in% hpo]

    if (nrow(hit) == 0L) {
        warning(
            "gvr_hpo_genes: no genes found for HPO term(s): ",
            paste(hpo, collapse = ", "),
            call. = FALSE
        )
        return(character(0))
    }

    genes <- unique(trimws(hit$Gene_Symbol))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    sort(unique(toupper(genes)))
}


.gvr_hpo_resolve_file <- function(hpo_path = NULL,
                                  hpo_url,
                                  cache_dir = NULL,
                                  refresh_cache = FALSE,
                                  verbose = TRUE) {
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

    if (!file.exists(dest) || isTRUE(refresh_cache)) {
        if (isTRUE(verbose)) {
            message("gvr_hpo_genes: downloading HPO phenotype-to-gene table")
        }

        ok <- tryCatch({
            utils::download.file(
                url = hpo_url,
                destfile = dest,
                mode = "wb",
                quiet = !isTRUE(verbose)
            )
            TRUE
        }, error = function(e) {
            if (file.exists(dest) && file.size(dest) == 0L) {
                unlink(dest)
            }
            stop(
                "gvr_hpo_genes: could not download HPO table from ",
                hpo_url,
                ". Use 'hpo_path=' with a local phenotype_to_genes.txt file. ",
                "Original error: ", conditionMessage(e),
                call. = FALSE
            )
        })

        if (!isTRUE(ok) || !file.exists(dest)) {
            stop("gvr_hpo_genes: HPO download failed.", call. = FALSE)
        }
    }

    dest
}


.gvr_hpo_read_table <- function(path) {
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("gvr_hpo_genes requires the 'data.table' package.", call. = FALSE)
    }

    dt <- data.table::fread(
        path,
        sep = "\t",
        header = TRUE,
        comment.char = "#",
        fill = TRUE,
        quote = "",
        data.table = TRUE
    )

    nm <- names(dt)

    hpo_col <- intersect(nm, c("HPO-ID", "HPO_ID", "hpo_id", "HPOId", "HPO_id"))
    gene_col <- intersect(nm, c("Gene-Name", "Gene_Name", "gene_symbol",
                                "GeneSymbol", "gene_name"))

    if (length(hpo_col) == 0L || length(gene_col) == 0L) {
        if (ncol(dt) >= 4L) {
            hpo_col <- nm[1L]
            gene_col <- nm[4L]
        } else {
            stop(
                "gvr_hpo_genes: could not identify HPO and gene-symbol columns ",
                "in HPO table.",
                call. = FALSE
            )
        }
    }

    out <- dt[, .(
        HPO_ID = as.character(get(hpo_col[1L])),
        Gene_Symbol = as.character(get(gene_col[1L]))
    )]

    out <- out[!is.na(HPO_ID) & nzchar(HPO_ID)]
    out <- out[!is.na(Gene_Symbol) & nzchar(Gene_Symbol)]

    out
}
