#' Enrich a germlinevaR MAF with clinical-significance predictions from MyVariant.info / dbNSFP
#'
#' @description
#' Augments the MAF returned by [read.gvr()] / [read.gvr.snpeff()] with one or
#' more clinical-significance prediction columns (default: AlphaMissense
#' categorical verdict and REVEL ensemble score), placed immediately after the
#' existing `CLIN_SIG` column. Predictions are looked up from the free public
#' MyVariant.info REST API (no key, no registration) backed by dbNSFP. An
#' on-disk SQLite cache makes re-runs offline and fast.
#'
#' `read.gvr()` and `read.gvr.snpeff()` are NOT modified by this enrichment
#' step. Their Phase N+2 byte-identity guarantee remains intact.
#'
#' @param maf A `data.table` / `data.frame` produced by [read.gvr()] or
#'   [read.gvr.snpeff()]. Must contain columns `NCBI_Build`, `Chromosome`,
#'   `Start_Position`, `Reference_Allele`, `Tumor_Seq_Allele2`,
#'   `Variant_Classification`, and (for default column placement) `CLIN_SIG`.
#' @param scope One of `"missense_splice"` (default), `"all"`, `"has_clinsig"`,
#'   or `"panel"`. Controls which rows are looked up against the API; rows
#'   outside scope receive `NA` for the new columns.
#' @param panel Character vector of gene symbols. Only used when
#'   `scope = "panel"`.
#' @param dbnsfp_columns Character vector of dbNSFP field tokens. Default:
#'   `c("AlphaMissense_pred", "REVEL_score")`. See Details for the full token
#'   list and the output-column naming convention.
#' @param cache_path Path to the SQLite cache file. If `NULL` (default),
#'   `tools::R_user_dir("germlinevaR", which = "cache") / "clinsig.sqlite"` is
#'   used. The parent directory is created if needed.
#' @param batch_size Variants per MyVariant.info POST (max 1000 per their docs;
#'   default 500 = safety margin).
#' @param max_retries Maximum HTTP retries per batch on transient failures.
#'   Default 3. Backoff: 1 s, 2 s, 4 s, 8 s (capped).
#' @param request_timeout HTTP timeout per request, in seconds. Default 60.
#' @param offline_only If `TRUE`, only consult the cache; never make network
#'   calls. Uncached rows receive `NA`. Default `FALSE`.
#' @param on_collision Behavior when the input MAF already contains columns
#'   that would be created by the enrichment step. `"overwrite"` (default)
#'   drops the existing columns with a warning, `"error"` fails loud.
#' @param verbose If `TRUE`, print progress to the console. Default `FALSE`.
#'
#' @details
#' \subsection{Supported `dbnsfp_columns` tokens}{
#'   \tabular{lll}{
#'     \strong{Token}                \tab \strong{Output column}                  \tab \strong{Type} \cr
#'     `"AlphaMissense_pred"`        \tab `CLIN_SIG_AlphaMissense`               \tab character (likely_pathogenic / ambiguous / likely_benign) \cr
#'     `"AlphaMissense_score"`       \tab `CLIN_SIG_AlphaMissense_score`         \tab numeric 0-1 \cr
#'     `"REVEL_score"`               \tab `CLIN_SIG_REVEL`                       \tab numeric 0-1 \cr
#'     `"REVEL_rankscore"`           \tab `CLIN_SIG_REVEL_rankscore`             \tab numeric 0-1 \cr
#'     `"CADD_phred"`                \tab `CLIN_SIG_CADD_phred`                  \tab numeric \cr
#'     `"ClinPred_pred"`             \tab `CLIN_SIG_ClinPred`                    \tab character (D/T) \cr
#'     `"MetaRNN_pred"`              \tab `CLIN_SIG_MetaRNN`                     \tab character (D/T) \cr
#'     `"BayesDel_addAF_pred"`       \tab `CLIN_SIG_BayesDel`                    \tab character (D/T) \cr
#'     `"PrimateAI_pred"`            \tab `CLIN_SIG_PrimateAI`                   \tab character (D/T)
#'   }
#' }
#'
#' \subsection{Build dispatch}{
#'   `unique(maf$NCBI_Build)` is mapped to MyVariant.info's `assembly`
#'   parameter:
#'   \itemize{
#'     \item `"GRCh37"` or `"hg19"` -> `assembly=hg19`
#'     \item `"GRCh38"` or `"hg38"` -> `assembly=hg38`
#'     \item Other labels (e.g. `"T2T-CHM13v2.0"`) -> error before any network call
#'   }
#' }
#'
#' \subsection{Scope semantics}{
#'   \itemize{
#'     \item `"missense_splice"` (default): `Variant_Classification` in
#'       `c("Missense_Mutation","Splice_Site","Splice_Region")`. Matches the
#'       coverage profile of dbNSFP (missense + splice SNVs).
#'     \item `"all"`: every row queried; non-missense rows are mostly `NA`
#'       because dbNSFP does not cover them.
#'     \item `"has_clinsig"`: missense_splice intersected with rows whose
#'       existing `CLIN_SIG` field is non-empty.
#'     \item `"panel"`: missense_splice intersected with
#'       `Hugo_Symbol %in% panel`.
#'   }
#' }
#'
#' @return A `data.table` with the same number of rows and the same row order
#' as `maf`, with new columns inserted immediately after `CLIN_SIG`. The
#' original `CLIN_SIG` column is preserved unchanged; downstream columns shift
#' right.
#'
#' @seealso [read.gvr()], [read.gvr.snpeff()]
#'
#' @export
enrich.gvr.clinsig <- function(maf,
                               scope          = c("missense_splice", "all", "has_clinsig", "panel"),
                               panel          = NULL,
                               dbnsfp_columns = c("AlphaMissense_pred", "REVEL_score"),
                               cache_path     = NULL,
                               batch_size     = 500L,
                               max_retries    = 3L,
                               request_timeout = 60L,
                               offline_only   = FALSE,
                               on_collision   = c("overwrite", "error"),
                               verbose        = FALSE) {

  ## ---------------------------------------------------------------------- ##
  ## 0. Input validation                                                    ##
  ## ---------------------------------------------------------------------- ##
  if (!is.data.frame(maf)) {
    stop("`maf` must be a data.frame or data.table (got ",
         paste(class(maf), collapse = "/"), ")", call. = FALSE)
  }
  if (nrow(maf) == 0L) {
    warning("`maf` has zero rows; returning input unchanged.", call. = FALSE)
    return(data.table::as.data.table(maf))
  }

  scope        <- match.arg(scope)
  on_collision <- match.arg(on_collision)

  required_cols <- c("NCBI_Build", "Chromosome", "Start_Position",
                     "Reference_Allele", "Tumor_Seq_Allele2",
                     "Variant_Classification")
  missing_cols  <- setdiff(required_cols, names(maf))
  if (length(missing_cols) > 0L) {
    stop("`maf` is missing required column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  ## ---------------------------------------------------------------------- ##
  ## 1. Validate dbnsfp_columns and build the field map                     ##
  ## ---------------------------------------------------------------------- ##
  fm <- .clinsig_field_map()
  bad_tokens <- setdiff(dbnsfp_columns, fm$token)
  if (length(bad_tokens) > 0L) {
    stop("Unknown dbnsfp_columns token(s): ",
         paste(bad_tokens, collapse = ", "), "\n",
         "Valid tokens are: ", paste(fm$token, collapse = ", "),
         call. = FALSE)
  }
  fm <- fm[match(dbnsfp_columns, fm$token), , drop = FALSE]

  ## ---------------------------------------------------------------------- ##
  ## 2. Build dispatch                                                      ##
  ## ---------------------------------------------------------------------- ##
  build_label <- unique(as.character(maf$NCBI_Build))
  build_label <- build_label[!is.na(build_label) & nzchar(build_label)]
  if (length(build_label) == 0L) {
    stop("MAF has no usable values in `NCBI_Build`.", call. = FALSE)
  }
  if (length(build_label) > 1L) {
    stop("MAF has multiple values in `NCBI_Build` (",
         paste(build_label, collapse = ", "),
         "); enrichment requires a single build.", call. = FALSE)
  }
  assembly_arg <- switch(build_label,
                         "GRCh37" = "hg19",
                         "hg19"   = "hg19",
                         "GRCh38" = "hg38",
                         "hg38"   = "hg38",
                         NA_character_)
  if (is.na(assembly_arg)) {
    stop("Unsupported genome build label `", build_label,
         "`. MyVariant.info supports only hg19/GRCh37 and hg38/GRCh38.",
         " Refusing to make any network call.", call. = FALSE)
  }
  build_canonical <- switch(assembly_arg, "hg19" = "GRCh37", "hg38" = "GRCh38")

  ## ---------------------------------------------------------------------- ##
  ## 3. Resolve scope to row indices                                        ##
  ## ---------------------------------------------------------------------- ##
  missense_splice_set <- c("Missense_Mutation", "Splice_Site", "Splice_Region")
  vc <- as.character(maf$Variant_Classification)
  in_scope <- switch(scope,
    "missense_splice" = vc %in% missense_splice_set,
    "all"             = rep(TRUE, nrow(maf)),
    "has_clinsig"     = {
      cs <- if ("CLIN_SIG" %in% names(maf)) as.character(maf$CLIN_SIG) else rep(NA_character_, nrow(maf))
      cs <- trimws(cs)
      has_cs <- !is.na(cs) & nzchar(cs) & !cs %in% c("-", ".")
      (vc %in% missense_splice_set) & has_cs
    },
    "panel" = {
      if (is.null(panel) || length(panel) == 0L) {
        stop("`scope = \"panel\"` requires a non-empty `panel` argument.", call. = FALSE)
      }
      if (!"Hugo_Symbol" %in% names(maf)) {
        stop("`scope = \"panel\"` requires the MAF to contain `Hugo_Symbol`.", call. = FALSE)
      }
      (vc %in% missense_splice_set) & (as.character(maf$Hugo_Symbol) %in% panel)
    }
  )
  in_scope[is.na(in_scope)] <- FALSE

  if (verbose) {
    message(sprintf("[enrich.gvr.clinsig] scope=%s  in_scope=%d / %d  build=%s assembly=%s",
                    scope, sum(in_scope), nrow(maf), build_canonical, assembly_arg))
  }

  ## ---------------------------------------------------------------------- ##
  ## 4. Handle column collisions before any work                            ##
  ## ---------------------------------------------------------------------- ##
  new_cols    <- fm$out_col
  collisions  <- intersect(new_cols, names(maf))
  if (length(collisions) > 0L) {
    if (on_collision == "error") {
      stop("MAF already contains output column(s): ",
           paste(collisions, collapse = ", "),
           ". Pass `on_collision = \"overwrite\"` to replace them.",
           call. = FALSE)
    }
    warning("Overwriting existing column(s): ",
            paste(collisions, collapse = ", "), call. = FALSE)
  }

  ## ---------------------------------------------------------------------- ##
  ## 5. Open the SQLite cache                                               ##
  ## ---------------------------------------------------------------------- ##
  cache_db <- .clinsig_cache_open(cache_path)
  on.exit(try(DBI::dbDisconnect(cache_db), silent = TRUE), add = TRUE)

  api_version <- .clinsig_api_version()

  ## ---------------------------------------------------------------------- ##
  ## 6. Build canonical variant keys for in-scope rows                      ##
  ## ---------------------------------------------------------------------- ##
  idx_scope <- which(in_scope)
  if (length(idx_scope) > 0L) {
    keys <- .clinsig_make_keys(
      chrom = as.character(maf$Chromosome[idx_scope]),
      pos   = as.integer(maf$Start_Position[idx_scope]),
      ref   = as.character(maf$Reference_Allele[idx_scope]),
      alt   = as.character(maf$Tumor_Seq_Allele2[idx_scope]),
      build = build_canonical
    )
  } else {
    keys <- data.frame(chrom = character(), pos = integer(),
                       ref = character(), alt = character(),
                       hgvs = character(), build = character(),
                       stringsAsFactors = FALSE)
  }

  ## ---------------------------------------------------------------------- ##
  ## 7. Cache lookup                                                        ##
  ## ---------------------------------------------------------------------- ##
  cached_raw <- if (nrow(keys) > 0L) {
    .clinsig_cache_get(cache_db, keys, api_version)
  } else {
    data.frame(chrom = character(), pos = integer(), ref = character(),
               alt = character(), payload = character(), status = character(),
               stringsAsFactors = FALSE)
  }

  ## Partial-hit detection via fields_fetched bookkeeping.
  ## A cached row is a real hit iff every currently-requested dbnsfp path is
  ## a member of the comma-separated `fields_fetched` set recorded at fetch
  ## time. Rows with status="not_in_dbnsfp" or "fetch_error" are valid
  ## hits regardless of field set (no dbNSFP entry exists for the variant,
  ## so requesting different fields cannot help).
  cached <- cached_raw
  partial_hits <- integer(0)
  if (nrow(cached_raw) > 0L) {
    requested <- fm$mvi_path
    partial_hits <- which(vapply(seq_len(nrow(cached_raw)), function(k) {
      if (cached_raw$status[k] != "ok") return(FALSE)
      have <- strsplit(cached_raw$fields_fetched[k], ",", fixed = TRUE)[[1]]
      !all(requested %in% have)
    }, logical(1L)))
    if (length(partial_hits) > 0L) {
      cached <- cached_raw[-partial_hits, , drop = FALSE]
    }
  }

  if (verbose) {
    message(sprintf("[enrich.gvr.clinsig] cache hits=%d / %d  (partial-hits requeued=%d)",
                    nrow(cached), nrow(keys), length(partial_hits)))
  }

  ## ---------------------------------------------------------------------- ##
  ## 8. Network fetch of misses                                             ##
  ## ---------------------------------------------------------------------- ##
  miss_keys <- if (nrow(keys) == 0L) {
    keys[0L, , drop = FALSE]
  } else {
    cached_idx <- with(cached,
                       paste(chrom, pos, ref, alt, sep = "|"))
    all_idx    <- with(keys,
                       paste(chrom, pos, ref, alt, sep = "|"))
    keys[!(all_idx %in% cached_idx), , drop = FALSE]
  }

  fetched <- data.frame(chrom = character(), pos = integer(), ref = character(),
                        alt = character(), payload = character(),
                        status = character(), stringsAsFactors = FALSE)

  if (nrow(miss_keys) > 0L && !offline_only) {
    # Fetch the union of all supported fields so cached payloads survive
    # changes in dbnsfp_columns across calls (plan's "union re-fetch" rule).
    fm_all <- .clinsig_field_map()
    fetched <- .clinsig_fetch_batches(
      keys           = miss_keys,
      assembly_arg   = assembly_arg,
      dbnsfp_paths   = fm_all$mvi_path,
      batch_size     = batch_size,
      max_retries    = max_retries,
      request_timeout = request_timeout,
      verbose        = verbose
    )
    if (nrow(fetched) > 0L) {
      .clinsig_cache_put(cache_db, fetched, build_canonical, api_version,
                         fields_fetched = paste(fm_all$mvi_path, collapse = ","))
    }
  } else if (nrow(miss_keys) > 0L && offline_only && verbose) {
    message(sprintf("[enrich.gvr.clinsig] offline_only=TRUE; %d miss(es) left as NA",
                    nrow(miss_keys)))
  }

  ## ---------------------------------------------------------------------- ##
  ## 9. Merge cached + fetched                                              ##
  ## ---------------------------------------------------------------------- ##
  ## `cached` carries an extra fields_fetched column used only for the
  ## partial-hit decision above; drop it so its schema matches `fetched`.
  if ("fields_fetched" %in% names(cached)) {
    cached <- cached[, setdiff(names(cached), "fields_fetched"), drop = FALSE]
  }
  combined <- rbind(cached, fetched, make.row.names = FALSE)

  ## ---------------------------------------------------------------------- ##
  ## 10. Build output columns and insert them after CLIN_SIG                ##
  ## ---------------------------------------------------------------------- ##
  parsed <- .clinsig_parse_payloads(combined$payload, fm, status = combined$status)
  new_col_values <- lapply(seq_len(nrow(fm)), function(i) {
    val <- if (fm$out_type[i] == "numeric") rep(NA_real_, nrow(maf))
           else rep(NA_character_, nrow(maf))
    if (nrow(combined) > 0L && length(idx_scope) > 0L) {
      key_str_maf <- paste(keys$chrom, keys$pos, keys$ref, keys$alt, sep = "|")
      key_str_cmb <- paste(combined$chrom, combined$pos, combined$ref, combined$alt, sep = "|")
      ord         <- match(key_str_maf, key_str_cmb)
      val[idx_scope] <- parsed[[i]][ord]
    }
    val
  })
  names(new_col_values) <- fm$out_col

  out <- data.table::as.data.table(maf)
  if (length(collisions) > 0L) {
    out[, (collisions) := NULL]
  }
  for (cn in names(new_col_values)) {
    out[, (cn) := new_col_values[[cn]]]
  }
  current_names <- names(out)
  new_added     <- names(new_col_values)
  others        <- setdiff(current_names, new_added)
  if ("CLIN_SIG" %in% others) {
    pos <- which(others == "CLIN_SIG")
    final_order <- c(others[seq_len(pos)], new_added,
                     if (pos < length(others)) others[(pos + 1L):length(others)] else character())
  } else {
    warning("MAF has no `CLIN_SIG` column; new columns appended at the end.",
            call. = FALSE)
    final_order <- c(others, new_added)
  }
  data.table::setcolorder(out, final_order)

  out
}


## =========================================================================
## Internal helpers
## =========================================================================

#' @keywords internal
#' @noRd
.clinsig_field_map <- function() {
  data.frame(
    token    = c("AlphaMissense_pred", "AlphaMissense_score",
                 "REVEL_score",         "REVEL_rankscore",
                 "CADD_phred",
                 "ClinPred_pred",       "MetaRNN_pred",
                 "BayesDel_addAF_pred", "PrimateAI_pred"),
    mvi_path = c("dbnsfp.alphamissense.pred",  "dbnsfp.alphamissense.score",
                 "dbnsfp.revel.score",         "dbnsfp.revel.rankscore",
                 "dbnsfp.cadd.phred",
                 "dbnsfp.clinpred.pred",       "dbnsfp.metarnn.pred",
                 "dbnsfp.bayesdel.add_af.pred","dbnsfp.primateai.pred"),
    out_col  = c("CLIN_SIG_AlphaMissense",       "CLIN_SIG_AlphaMissense_score",
                 "CLIN_SIG_REVEL",               "CLIN_SIG_REVEL_rankscore",
                 "CLIN_SIG_CADD_phred",
                 "CLIN_SIG_ClinPred",            "CLIN_SIG_MetaRNN",
                 "CLIN_SIG_BayesDel",            "CLIN_SIG_PrimateAI"),
    out_type = c("character", "numeric",
                 "numeric",   "numeric",
                 "numeric",
                 "character", "character",
                 "character", "character"),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
.clinsig_api_version <- function() {
  paste0("myvariant_v1_", format(Sys.Date(), "%Y-%m"))
}

#' @keywords internal
#' @noRd
.clinsig_cache_open <- function(cache_path = NULL) {
  if (is.null(cache_path) || !nzchar(cache_path)) {
    cdir <- tools::R_user_dir("germlinevaR", which = "cache")
    if (!dir.exists(cdir)) dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
    cache_path <- file.path(cdir, "clinsig.sqlite")
  } else {
    parent <- dirname(cache_path)
    if (nzchar(parent) && !dir.exists(parent)) {
      dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    }
  }
  con <- DBI::dbConnect(RSQLite::SQLite(), cache_path)
  try(DBI::dbExecute(con, "PRAGMA journal_mode = WAL;"),    silent = TRUE)
  try(DBI::dbExecute(con, "PRAGMA synchronous = NORMAL;"),  silent = TRUE)
  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS myvariant_dbnsfp_cache (
      build          TEXT NOT NULL,
      chrom          TEXT NOT NULL,
      pos            INTEGER NOT NULL,
      ref            TEXT NOT NULL,
      alt            TEXT NOT NULL,
      api_version    TEXT NOT NULL,
      payload        TEXT NOT NULL,
      status         TEXT NOT NULL,
      fields_fetched TEXT NOT NULL DEFAULT '',
      cached_at      INTEGER NOT NULL,
      PRIMARY KEY (build, chrom, pos, ref, alt, api_version)
    );")
  ## In case the table pre-existed without the fields_fetched column (legacy),
  ## add it. Safe to swallow the error if the column already exists.
  try(DBI::dbExecute(con,
    "ALTER TABLE myvariant_dbnsfp_cache ADD COLUMN fields_fetched TEXT NOT NULL DEFAULT '';"),
    silent = TRUE)
  DBI::dbExecute(con, "
    CREATE INDEX IF NOT EXISTS idx_lookup
      ON myvariant_dbnsfp_cache (build, chrom, pos);")
  con
}

#' @keywords internal
#' @noRd
.clinsig_make_keys <- function(chrom, pos, ref, alt, build) {
  chrom_norm <- ifelse(grepl("^chr", chrom, ignore.case = TRUE),
                       chrom, paste0("chr", chrom))
  hgvs <- ifelse(
    nchar(ref) == 1L & nchar(alt) == 1L & ref != "-" & alt != "-",
    sprintf("%s:g.%d%s>%s", chrom_norm, pos, ref, alt),
    NA_character_
  )
  data.frame(
    chrom = chrom_norm,
    pos   = as.integer(pos),
    ref   = ref,
    alt   = alt,
    hgvs  = hgvs,
    build = build,
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
.clinsig_cache_get <- function(con, keys, api_version) {
  if (nrow(keys) == 0L) {
    return(data.frame(chrom = character(), pos = integer(), ref = character(),
                      alt = character(), payload = character(),
                      status = character(), fields_fetched = character(),
                      stringsAsFactors = FALSE))
  }
  build <- keys$build[1L]
  # Stage the lookup keys in a temporary table and inner-join. This is
  # dramatically faster than chunked WHERE OR queries when many keys are
  # in play (148K keys: ~1 s vs ~800 s for chunked OR).
  tmp_keys <- data.frame(
    chrom = as.character(keys$chrom),
    pos   = as.integer(keys$pos),
    ref   = as.character(keys$ref),
    alt   = as.character(keys$alt),
    stringsAsFactors = FALSE
  )
  DBI::dbWriteTable(con, "__clinsig_lookup", tmp_keys,
                    temporary = TRUE, overwrite = TRUE)
  on.exit(try(DBI::dbExecute(con,
              "DROP TABLE IF EXISTS __clinsig_lookup"), silent = TRUE),
          add = TRUE)
  sql <- "
    SELECT c.chrom AS chrom, c.pos AS pos, c.ref AS ref, c.alt AS alt,
           c.payload AS payload, c.status AS status,
           c.fields_fetched AS fields_fetched
    FROM myvariant_dbnsfp_cache c
    INNER JOIN __clinsig_lookup k
      ON c.chrom = k.chrom AND c.pos = k.pos
         AND c.ref = k.ref AND c.alt = k.alt
    WHERE c.build = ? AND c.api_version = ?"
  DBI::dbGetQuery(con, sql, params = list(build, api_version))
}

#' @keywords internal
#' @noRd
.clinsig_cache_put <- function(con, fetched, build, api_version, fields_fetched) {
  if (nrow(fetched) == 0L) return(invisible(NULL))
  now <- as.integer(Sys.time())
  df <- data.frame(
    build           = build,
    chrom           = fetched$chrom,
    pos             = as.integer(fetched$pos),
    ref             = fetched$ref,
    alt             = fetched$alt,
    api_version     = api_version,
    payload         = fetched$payload,
    status          = fetched$status,
    fields_fetched  = fields_fetched,
    cached_at       = now,
    stringsAsFactors = FALSE
  )
  DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, "
      INSERT OR REPLACE INTO myvariant_dbnsfp_cache
        (build, chrom, pos, ref, alt, api_version, payload, status, fields_fetched, cached_at)
      VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
      params = unname(as.list(df)))
  })
  invisible(NULL)
}

#' Walk a dot-notation path into a nested list
#'
#' @keywords internal
#' @noRd
.clinsig_nav <- function(obj, path) {
  parts <- strsplit(path, ".", fixed = TRUE)[[1]]
  for (p in parts) {
    if (is.null(obj)) return(NULL)
    if (!is.list(obj)) return(NULL)
    obj <- obj[[p]]
  }
  obj
}

#' @keywords internal
#' @noRd
.clinsig_parse_payloads <- function(payloads, fm, status = NULL) {
  out <- vector("list", nrow(fm))
  names(out) <- fm$out_col
  for (i in seq_len(nrow(fm))) {
    out[[i]] <- if (fm$out_type[i] == "numeric") rep(NA_real_, length(payloads))
                else rep(NA_character_, length(payloads))
  }
  if (length(payloads) == 0L) return(out)
  # If status is provided, only parse rows with status=="ok"; all others
  # carry no dbNSFP data and must stay NA.
  parse_idx <- if (is.null(status)) seq_along(payloads) else which(status == "ok")
  for (k in parse_idx) {
    p <- payloads[k]
    if (is.na(p) || !nzchar(p) || p == "{}") next
    obj <- tryCatch(jsonlite::fromJSON(p, simplifyVector = FALSE),
                    error = function(e) NULL)
    if (is.null(obj)) next
    for (i in seq_len(nrow(fm))) {
      v <- .clinsig_nav(obj, fm$mvi_path[i])
      if (is.null(v)) next
      if (is.list(v)) v <- unlist(v, use.names = FALSE)
      if (length(v) == 0L) next
      if (length(v) > 1L) v <- v[1L]
      if (fm$out_type[i] == "numeric") {
        out[[i]][k] <- suppressWarnings(as.numeric(v))
      } else {
        out[[i]][k] <- .clinsig_normalize_pred(fm$token[i], as.character(v))
      }
    }
  }
  out
}

#' @keywords internal
#' @noRd
.clinsig_normalize_pred <- function(token, value) {
  if (is.null(value) || is.na(value) || !nzchar(value)) return(NA_character_)
  if (token == "AlphaMissense_pred") {
    return(switch(toupper(value),
                  "LIKELY_PATHOGENIC" = "likely_pathogenic",
                  "LIKELY_BENIGN"     = "likely_benign",
                  "AMBIGUOUS"         = "ambiguous",
                  "LPA"               = "likely_pathogenic",
                  "LBE"               = "likely_benign",
                  "AMB"               = "ambiguous",
                  "P"                 = "likely_pathogenic",
                  "B"                 = "likely_benign",
                  "A"                 = "ambiguous",
                  tolower(value)))
  }
  value
}

#' @keywords internal
#' @noRd
.clinsig_fetch_batches <- function(keys, assembly_arg, dbnsfp_paths,
                                   batch_size, max_retries, request_timeout,
                                   verbose) {
  if (nrow(keys) == 0L) {
    return(data.frame(chrom = character(), pos = integer(), ref = character(),
                      alt = character(), payload = character(),
                      status = character(), stringsAsFactors = FALSE))
  }
  # Rows whose HGVS we cannot build (NA hgvs, e.g. indels) -> status not_in_dbnsfp
  no_hgvs <- is.na(keys$hgvs) | !nzchar(keys$hgvs)
  out_list <- list()
  if (any(no_hgvs)) {
    sub <- keys[no_hgvs, , drop = FALSE]
    out_list[[length(out_list) + 1L]] <- data.frame(
      chrom   = sub$chrom,
      pos     = sub$pos,
      ref     = sub$ref,
      alt     = sub$alt,
      payload = "{}",
      status  = "not_in_dbnsfp",
      stringsAsFactors = FALSE
    )
  }
  queryable <- keys[!no_hgvs, , drop = FALSE]
  if (nrow(queryable) == 0L) {
    return(do.call(rbind, c(out_list, list(make.row.names = FALSE))))
  }

  fields_str <- paste(dbnsfp_paths, collapse = ",")
  n          <- nrow(queryable)
  n_batches  <- ceiling(n / batch_size)
  if (verbose) {
    message(sprintf("[enrich.gvr.clinsig] POST batches: %d (batch_size=%d, total=%d)",
                    n_batches, batch_size, n))
  }

  for (b in seq_len(n_batches)) {
    i0 <- (b - 1L) * batch_size + 1L
    i1 <- min(b * batch_size, n)
    chunk <- queryable[i0:i1, , drop = FALSE]
    ids_str <- paste(chunk$hgvs, collapse = ",")

    body <- .clinsig_post_with_retry(
      ids_str         = ids_str,
      fields_str      = fields_str,
      assembly_arg    = assembly_arg,
      max_retries     = max_retries,
      request_timeout = request_timeout,
      verbose         = verbose,
      batch_num       = b,
      batch_total     = n_batches
    )

    if (is.null(body)) {
      # Persistent failure -> mark this batch's rows with status="fetch_error"
      out_list[[length(out_list) + 1L]] <- data.frame(
        chrom   = chunk$chrom,
        pos     = chunk$pos,
        ref     = chunk$ref,
        alt     = chunk$alt,
        payload = "{}",
        status  = "fetch_error",
        stringsAsFactors = FALSE
      )
      next
    }

    # Parse the JSON array; each element keyed by `query` == HGVS
    arr <- tryCatch(jsonlite::fromJSON(body, simplifyVector = FALSE),
                    error = function(e) NULL)
    if (is.null(arr) || !is.list(arr)) {
      out_list[[length(out_list) + 1L]] <- data.frame(
        chrom   = chunk$chrom,
        pos     = chunk$pos,
        ref     = chunk$ref,
        alt     = chunk$alt,
        payload = "{}",
        status  = "fetch_error",
        stringsAsFactors = FALSE
      )
      next
    }
    # MyVariant.info MAY return objects in a different order than requested,
    # and MAY return multiple objects per ID if the variant is multi-allelic.
    # Index by `query` field with the first match wins.
    by_query <- list()
    for (el in arr) {
      q <- el$query
      if (is.null(q) || is.na(q) || !nzchar(q)) next
      if (is.null(by_query[[q]])) by_query[[q]] <- el
    }

    payloads <- vapply(chunk$hgvs, function(h) {
      el <- by_query[[h]]
      if (is.null(el)) return("{}")
      jsonlite::toJSON(el, auto_unbox = TRUE, null = "null", na = "null")
    }, character(1L))

    statuses <- vapply(chunk$hgvs, function(h) {
      el <- by_query[[h]]
      if (is.null(el))            return("fetch_error")
      if (isTRUE(el$notfound))    return("not_in_dbnsfp")
      if (is.null(el$dbnsfp))     return("not_in_dbnsfp")
      "ok"
    }, character(1L))

    out_list[[length(out_list) + 1L]] <- data.frame(
      chrom   = chunk$chrom,
      pos     = chunk$pos,
      ref     = chunk$ref,
      alt     = chunk$alt,
      payload = unname(payloads),
      status  = unname(statuses),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, c(out_list, list(make.row.names = FALSE)))
}

#' @keywords internal
#' @noRd
.clinsig_post_with_retry <- function(ids_str, fields_str, assembly_arg,
                                     max_retries, request_timeout, verbose,
                                     batch_num, batch_total) {
  url <- "https://myvariant.info/v1/variant"
  for (attempt in seq_len(max_retries + 1L)) {
    req <- httr2::request(url) |>
           httr2::req_method("POST") |>
           httr2::req_body_form(ids      = ids_str,
                                fields   = fields_str,
                                assembly = assembly_arg) |>
           httr2::req_headers(Accept = "application/json") |>
           httr2::req_user_agent("germlinevaR/0.1 (Phase N+3)") |>
           httr2::req_timeout(request_timeout) |>
           httr2::req_error(is_error = function(resp) FALSE)
    resp <- tryCatch(httr2::req_perform(req), error = function(e) e)
    if (inherits(resp, "error")) {
      if (verbose) message(sprintf("[enrich.gvr.clinsig] batch %d/%d attempt %d -> error: %s",
                                   batch_num, batch_total, attempt,
                                   conditionMessage(resp)))
    } else {
      sc <- httr2::resp_status(resp)
      if (sc >= 200L && sc < 300L) {
        return(httr2::resp_body_string(resp))
      }
      if (verbose) message(sprintf("[enrich.gvr.clinsig] batch %d/%d attempt %d -> HTTP %d",
                                   batch_num, batch_total, attempt, sc))
      # 4xx (other than 429) are non-retriable
      if (sc >= 400L && sc < 500L && sc != 429L) {
        warning(sprintf("MyVariant.info batch %d/%d failed with HTTP %d (non-retriable); marking rows as fetch_error.",
                        batch_num, batch_total, sc), call. = FALSE)
        return(NULL)
      }
    }
    if (attempt > max_retries) break
    wait_sec <- min(2^(attempt - 1L), 8L)
    Sys.sleep(wait_sec)
  }
  warning(sprintf("MyVariant.info batch %d/%d failed after %d attempts; marking rows as fetch_error.",
                  batch_num, batch_total, max_retries + 1L), call. = FALSE)
  NULL
}
