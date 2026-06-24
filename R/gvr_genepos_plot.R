# ----------------------------------------------------------------------------
# gvr_genepos.plot  --  gene-track lollipop on a cDNA axis
#
# A per-gene "lollipop on the gene structure" plot:
#   * exon / intron / UTR track for one transcript (Ensembl REST + cache,
#     optional GTF override)
#   * one stem-and-dot per variant, placed on a cDNA axis from HGVSc,
#     colored by Variant_Classification
#   * same visual idiom as gvr_lollipop() (bar_half = 0.30, dom_half = 0.60,
#     same theme, same hotspot wash)
#
# Companion to gvr_lollipop(). Where gvr_lollipop draws protein domains on
# a protein-position axis, gvr_genepos.plot draws gene structure on a cDNA
# axis.
# ----------------------------------------------------------------------------

# ---- Internal helpers (file-scope, not exported) ---------------------------

# FUSE-safe writer; mirrors .fuse_save in gvr_lollipop. Necessary because
# direct writes to /mnt/results/* via grDevices/ggsave can produce zero-byte
# files on the S3-backed FUSE mount; we render to tempdir then `cp` over.
.gvr_genepos_fuse_save <- function(final_path, draw_fun) {
  tmp <- file.path(tempdir(), basename(final_path))
  ok <- tryCatch(
    { draw_fun(tmp); file.exists(tmp) && file.info(tmp)$size > 0 },
    error = function(e) {
      warning(sprintf("gvr_genepos.plot: render failed (%s): %s",
                      basename(final_path), conditionMessage(e)))
      FALSE
    })
  if (!ok) return(NA_character_)
  system2("cp", c("-f", shQuote(tmp), shQuote(final_path)))
  if (!file.exists(final_path) || file.info(final_path)$size == 0) {
    warning(sprintf("gvr_genepos.plot: copy to '%s' failed; left at '%s'.",
                    final_path, tmp))
    return(tmp)
  }
  final_path
}

# Cache directory resolver: same precedence chain as gvr_lollipop /
# gvr_domain_cache_clear. Returns a writable directory path, or
# NA_character_ if no disk caching should be performed (cache_dir == FALSE).
.gvr_genepos_resolve_cache_dir <- function(cache_dir_arg) {
  if (isFALSE(cache_dir_arg)) return(NA_character_)

  candidates <- character(0)
  if (!is.null(cache_dir_arg) && is.character(cache_dir_arg) &&
      length(cache_dir_arg) == 1L && nzchar(cache_dir_arg)) {
    candidates <- c(candidates, cache_dir_arg)
  }
  env_dir <- Sys.getenv("GVR_CACHE_DIR", unset = "")
  if (nzchar(env_dir)) candidates <- c(candidates, env_dir)
  opt_dir <- getOption("germlinevaR.cache_dir", default = NULL)
  if (!is.null(opt_dir) && is.character(opt_dir) && length(opt_dir) == 1L &&
      nzchar(opt_dir)) {
    candidates <- c(candidates, opt_dir)
  }
  candidates <- c(candidates,
                  tools::R_user_dir("germlinevaR", which = "cache"),
                  file.path(tempdir(), "germlinevaR_cache"))

  for (cand in candidates) {
    ok <- tryCatch({
      if (!dir.exists(cand))
        dir.create(cand, recursive = TRUE, showWarnings = FALSE)
      dir.exists(cand) && file.access(cand, mode = 2L) == 0L
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (isTRUE(ok)) return(cand)
  }
  NA_character_
}

# Cache file path for a transcript structure.
.gvr_genestruct_cache_path <- function(enst, assembly, cache_dir) {
  if (is.na(cache_dir)) return(NA_character_)
  file.path(cache_dir,
            sprintf("genestruct_ensembl_%s_%s.rds", enst, assembly))
}

# HTTP GET with retry (matches gvr_lollipop's .gvr_http_get_retry).
# Defined here at file scope so we don't depend on .gvr_http_get_retry
# being in scope from another file's nested closure.
.gvr_genepos_http_get_retry <- function(url, timeout_s = 30, tries = 3) {
  delays <- c(0.5, 1.5, 3.0)
  last_err <- NULL
  resp <- NULL
  for (i in seq_len(tries)) {
    resp <- tryCatch(
      httr::GET(url, httr::timeout(timeout_s)),
      error = function(e) { last_err <<- e; NULL }
    )
    if (!is.null(resp)) {
      sc <- httr::status_code(resp)
      if (!(sc %in% c(429L, 500L, 502L, 503L, 504L))) return(resp)
    }
    if (i < tries) Sys.sleep(delays[i])
  }
  if (!is.null(resp)) return(resp)
  stop(last_err)
}

# Parse the Ensembl REST /lookup/id response into the lean internal
# gene-structure list documented in PLAN section 3d.
.gvr_genestruct_parse_rest <- function(payload, enst, assembly) {
  exons_in <- payload$Exon
  if (is.null(exons_in) || length(exons_in) == 0L)
    stop(sprintf("gvr_genepos.plot: Ensembl REST returned no exons for %s.", enst))

  # Convert exon list to data.table; sort by genomic start.
  ex <- data.table::rbindlist(lapply(exons_in, function(e) {
    data.table::data.table(
      g_start   = as.integer(e$start),
      g_end     = as.integer(e$end),
      length_bp = as.integer(e$end) - as.integer(e$start) + 1L,
      exon_id   = if (is.null(e$id)) NA_character_ else as.character(e$id)
    )
  }))
  data.table::setorder(ex, g_start)

  # On minus-strand transcripts, the 5' exon is at the highest genomic
  # coordinate. We number exons in 5'->3' transcript order regardless of
  # strand, so the user always sees exon[1] at the visual left of the plot.
  strand <- as.integer(payload$strand)
  if (is.na(strand) || !(strand %in% c(1L, -1L))) {
    stop(sprintf("gvr_genepos.plot: Ensembl REST returned non-canonical strand for %s.",
                 enst))
  }
  if (strand == -1L) ex <- ex[order(-g_start)]
  ex[, exon_idx := seq_len(.N)]
  data.table::setcolorder(ex, c("exon_idx", "g_start", "g_end", "length_bp", "exon_id"))

  # CDS range from Translation (NULL for non-coding transcripts).
  cds_range <- NULL
  if (!is.null(payload$Translation)) {
    cds_range <- c(as.integer(payload$Translation$start),
                   as.integer(payload$Translation$end))
  }

  # UTRs: list of objects with type=five_prime_UTR/three_prime_UTR + start+end.
  utr5 <- data.table::data.table(g_start = integer(), g_end = integer(),
                                 length_bp = integer())
  utr3 <- data.table::data.table(g_start = integer(), g_end = integer(),
                                 length_bp = integer())
  if (!is.null(payload$UTR) && length(payload$UTR) > 0L) {
    for (u in payload$UTR) {
      typ <- as.character(u$object_type %||% u$type)
      typ <- tolower(typ)
      gs  <- as.integer(u$start)
      ge  <- as.integer(u$end)
      row <- data.table::data.table(g_start = gs, g_end = ge,
                                    length_bp = ge - gs + 1L)
      if (grepl("five", typ, fixed = TRUE) || grepl("5", typ, fixed = TRUE))
        utr5 <- data.table::rbindlist(list(utr5, row))
      else if (grepl("three", typ, fixed = TRUE) || grepl("3", typ, fixed = TRUE))
        utr3 <- data.table::rbindlist(list(utr3, row))
    }
  }
  # Order UTR segments in transcript 5'->3' order (consistent with exons).
  if (nrow(utr5) > 0L) {
    if (strand == -1L) data.table::setorder(utr5, -g_start) else data.table::setorder(utr5, g_start)
  }
  if (nrow(utr3) > 0L) {
    if (strand == -1L) data.table::setorder(utr3, -g_start) else data.table::setorder(utr3, g_start)
  }

  list(
    transcript_id = as.character(enst),
    gene_symbol   = as.character(payload$display_name %||% NA_character_),
    assembly      = as.character(assembly),
    strand        = strand,
    seq_region    = as.character(payload$seq_region_name %||% NA_character_),
    is_canonical  = isTRUE(as.logical(payload$is_canonical)),
    exons         = ex,
    cds_range     = cds_range,
    utr5          = utr5,
    utr3          = utr3,
    cdna_len      = as.integer(sum(ex$length_bp))
  )
}

# Local null-coalesce (file-scope helper used in parser above).
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1L && is.na(a))) b else a

# Fetch transcript structure: cache hit -> reuse; else REST + cache.
.gvr_genestruct_fetch <- function(enst, assembly, cache_dir, verbose = TRUE) {
  cache_path <- .gvr_genestruct_cache_path(enst, assembly, cache_dir)
  if (!is.na(cache_path) && file.exists(cache_path)) {
    if (isTRUE(verbose))
      message(sprintf("gvr_genepos.plot: cache hit %s", cache_path))
    return(readRDS(cache_path))
  }

  if (!requireNamespace("httr", quietly = TRUE))
    stop("gvr_genepos.plot: package 'httr' is required for Ensembl REST. ",
         "Install it or supply 'gtf_path' for offline use.")
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("gvr_genepos.plot: package 'jsonlite' is required to parse Ensembl REST.")

  host <- if (identical(assembly, "GRCh37")) "https://grch37.rest.ensembl.org" else "https://rest.ensembl.org"
  if (!identical(assembly, "GRCh37") && !identical(assembly, "GRCh38")) {
    stop(sprintf(
      "gvr_genepos.plot: unsupported NCBI_Build '%s'. Supply 'gtf_path' for offline lookup.",
      assembly))
  }
  url <- sprintf("%s/lookup/id/%s?expand=1;utr=1", host, enst)
  if (isTRUE(verbose))
    message(sprintf("gvr_genepos.plot: REST GET %s", url))

  resp <- .gvr_genepos_http_get_retry(url, timeout_s = 30, tries = 3)
  sc <- httr::status_code(resp)
  if (sc != 200L) {
    body <- tryCatch(httr::content(resp, "text", encoding = "UTF-8"),
                     error = function(e) "")
    stop(sprintf("gvr_genepos.plot: Ensembl REST returned HTTP %d for %s.\n%s",
                 sc, enst, substr(body, 1, 200)))
  }
  payload <- httr::content(resp, "parsed", type = "application/json",
                           simplifyVector = FALSE)
  out <- .gvr_genestruct_parse_rest(payload, enst, assembly)

  if (!is.na(cache_path)) {
    tryCatch({
      saveRDS(out, cache_path)
      if (isTRUE(verbose))
        message(sprintf("gvr_genepos.plot: cached -> %s", cache_path))
    }, error = function(e) {
      warning(sprintf("gvr_genepos.plot: cache write failed: %s",
                      conditionMessage(e)))
    })
  }
  out
}

# Offline GTF parser: extract exon/UTR/CDS records for one transcript.
# Prefers rtracklayer::import (Suggests); falls back to a streaming line
# parser on key=value attribute strings.
.gvr_genestruct_from_gtf <- function(gtf_path, enst, assembly = "GRCh38") {
  if (!file.exists(gtf_path))
    stop(sprintf("gvr_genepos.plot: gtf_path does not exist: %s", gtf_path))

  if (requireNamespace("rtracklayer", quietly = TRUE)) {
    gr <- rtracklayer::import(gtf_path)
    df <- as.data.frame(gr, stringsAsFactors = FALSE)
    df <- df[!is.na(df$transcript_id) & df$transcript_id == enst, , drop = FALSE]
    if (nrow(df) == 0L)
      stop(sprintf("gvr_genepos.plot: transcript_id '%s' not found in GTF %s.",
                   enst, gtf_path))
    feat_col <- "type"
  } else {
    # Streaming line parse on the GTF.
    lns <- readLines(gtf_path, warn = FALSE)
    lns <- lns[!grepl("^#", lns)]
    lns <- lns[grepl(sprintf("transcript_id \"%s\"", enst), lns, fixed = TRUE)]
    if (length(lns) == 0L)
      stop(sprintf("gvr_genepos.plot: transcript_id '%s' not found in GTF %s.",
                   enst, gtf_path))
    parts <- strsplit(lns, "\t", fixed = TRUE)
    seqid    <- vapply(parts, `[`, character(1), 1)
    type     <- vapply(parts, `[`, character(1), 3)
    start_v  <- as.integer(vapply(parts, `[`, character(1), 4))
    end_v    <- as.integer(vapply(parts, `[`, character(1), 5))
    strand_v <- vapply(parts, `[`, character(1), 7)
    df <- data.frame(seqid = seqid, type = type,
                     start = start_v, end = end_v, strand = strand_v,
                     stringsAsFactors = FALSE)
    feat_col <- "type"
  }

  exons_df <- df[df[[feat_col]] == "exon", , drop = FALSE]
  if (nrow(exons_df) == 0L)
    stop(sprintf("gvr_genepos.plot: no exon rows in GTF for %s.", enst))
  strand_chr <- as.character(exons_df$strand[1])
  strand_i   <- if (strand_chr == "-") -1L else 1L

  ex <- data.table::data.table(
    g_start   = as.integer(exons_df$start),
    g_end     = as.integer(exons_df$end),
    length_bp = as.integer(exons_df$end) - as.integer(exons_df$start) + 1L
  )
  data.table::setorder(ex, g_start)
  if (strand_i == -1L) ex <- ex[order(-g_start)]
  ex[, exon_idx := seq_len(.N)]
  data.table::setcolorder(ex, c("exon_idx", "g_start", "g_end", "length_bp"))

  cds_df <- df[df[[feat_col]] == "CDS", , drop = FALSE]
  cds_range <- if (nrow(cds_df) > 0L) c(min(cds_df$start), max(cds_df$end)) else NULL

  mk <- function(types) {
    mdt <- df[df[[feat_col]] %in% types, , drop = FALSE]
    if (nrow(mdt) == 0L)
      return(data.table::data.table(g_start = integer(), g_end = integer(),
                                    length_bp = integer()))
    out <- data.table::data.table(g_start = as.integer(mdt$start),
                                  g_end   = as.integer(mdt$end))
    out[, length_bp := g_end - g_start + 1L]
    if (strand_i == -1L) data.table::setorder(out, -g_start) else data.table::setorder(out, g_start)
    out
  }
  utr5 <- mk(c("five_prime_utr", "5UTR", "5'UTR"))
  utr3 <- mk(c("three_prime_utr", "3UTR", "3'UTR"))

  list(
    transcript_id = as.character(enst),
    gene_symbol   = NA_character_,
    assembly      = as.character(assembly),
    strand        = strand_i,
    seq_region    = as.character(exons_df$seqid[1] %||% NA_character_),
    is_canonical  = NA,
    exons         = ex,
    cds_range     = cds_range,
    utr5          = utr5,
    utr3          = utr3,
    cdna_len      = as.integer(sum(ex$length_bp))
  )
}

# Parse a HGVS coding-DNA string into (kind, pos, valid). Anchors only the
# leading integer; ignores insertion/deletion suffix lengths. PLAN 3f.
.gvr_parse_hgvsc <- function(hgvsc) {
  n <- length(hgvsc)
  out <- data.table::data.table(
    kind  = rep(NA_character_, n),
    pos   = rep(NA_integer_,   n),
    valid = rep(FALSE,         n)
  )
  if (n == 0L) return(out)

  # Strip transcript prefix "ENST...c." / "NM_...:c." / leading "c."
  s <- as.character(hgvsc)
  s <- sub("^[^:]+:", "", s)        # drop transcript:
  s <- sub("^c\\.", "", s)          # drop c.
  s <- trimws(s)

  # Splice / deep-intronic: contains '+K' or '-K' AFTER the leading anchor
  # but BEFORE any letter -- e.g. "123+2A>G", "456-1G>T". The anchor pos is
  # still the leading integer, but we mark these as splice and skip.
  is_splice <- grepl("^-?\\*?\\d+[+\\-]\\d+", s)

  # 5'UTR: starts with '-' followed by digits.
  is_utr5 <- grepl("^-\\d+", s) & !is_splice

  # 3'UTR: starts with '*' followed by digits.
  is_utr3 <- grepl("^\\*\\d+", s) & !is_splice

  # CDS: starts with a positive digit.
  is_cds <- grepl("^\\d+", s) & !is_splice & !is_utr5 & !is_utr3

  # Extract leading integer (preserving sign for utr5, dropping '*' for utr3).
  utr5_pos <- suppressWarnings(as.integer(sub("^(-\\d+).*$", "\\1", s)))
  utr3_pos <- suppressWarnings(as.integer(sub("^\\*(\\d+).*$", "\\1", s)))
  cds_pos  <- suppressWarnings(as.integer(sub("^(\\d+).*$", "\\1", s)))

  out[is_splice, `:=`(kind = "splice", valid = FALSE)]
  out[is_utr5, `:=`(kind = "utr5", pos = utr5_pos[is_utr5], valid = TRUE)]
  out[is_utr3, `:=`(kind = "utr3", pos = utr3_pos[is_utr3], valid = TRUE)]
  out[is_cds,  `:=`(kind = "cds",  pos = cds_pos[is_cds],  valid = TRUE)]

  # Any row that didn't match any of the above stays kind=NA, valid=FALSE
  # (unparsed). Also flag valid=FALSE for any row that parsed but produced NA.
  out[is.na(pos) & valid == TRUE, valid := FALSE]
  out
}

# Build the per-region layout table on the cDNA axis. Returns a list with
# elements:
#   regions  -- data.table of (exon_idx, region_kind, g_start, g_end,
#               x_start, x_end, fill) where region_kind in
#               {"utr5","cds","utr3","intron"}.
#   anchors  -- numeric vector indexed by exon_idx giving the cDNA-axis
#               offset where each exon starts (used to map HGVSc CDS pos
#               to x). Specifically anchors[i] = visual x of the first base
#               of exon i (post-intron-scaling).
#   cds_len, utr5_len, utr3_len, x_atg (x of c.1), x_max.
.gvr_genepos_layout <- function(struct,
                                intron_scale     = "fixed",
                                intron_visual_bp = 200L,
                                utr_visual_bp    = NULL) {
  ex <- struct$exons
  n_ex <- nrow(ex)
  if (n_ex == 0L) stop("gvr_genepos.plot: structure has 0 exons.")

  # Determine per-exon CDS/UTR overlap. We work on genomic coordinates first
  # then map back to per-exon-segment cDNA widths.
  cds_g_lo <- if (!is.null(struct$cds_range)) min(struct$cds_range) else NA_integer_
  cds_g_hi <- if (!is.null(struct$cds_range)) max(struct$cds_range) else NA_integer_
  has_cds  <- !is.na(cds_g_lo) && !is.na(cds_g_hi)

  # For each exon, split its real cDNA length into utr5 / cds / utr3 pieces.
  # Note: at exon-level granularity 5'UTR/CDS/3'UTR can co-exist in a single
  # exon (terminal exons); we split on genomic overlap.
  ex_split <- vector("list", n_ex)
  for (i in seq_len(n_ex)) {
    g_lo <- ex$g_start[i]; g_hi <- ex$g_end[i]
    bp <- g_hi - g_lo + 1L
    if (!has_cds) {
      # Non-coding transcript: all exon length is "cds" by convention
      # (we won't plot UTR but the axis still works).
      ex_split[[i]] <- list(utr5_bp = 0L, cds_bp = bp, utr3_bp = 0L)
      next
    }
    # 3 cases for the exon relative to the genomic CDS interval.
    if (g_hi < cds_g_lo || g_lo > cds_g_hi) {
      # entirely outside CDS -> all UTR. Side depends on transcript strand.
      side_5p <- if (struct$strand == 1L) (g_hi < cds_g_lo) else (g_lo > cds_g_hi)
      if (side_5p)
        ex_split[[i]] <- list(utr5_bp = bp, cds_bp = 0L, utr3_bp = 0L)
      else
        ex_split[[i]] <- list(utr5_bp = 0L, cds_bp = 0L, utr3_bp = bp)
      next
    }
    # exon overlaps CDS: split into up to 3 pieces.
    cds_lo_in <- max(g_lo, cds_g_lo)
    cds_hi_in <- min(g_hi, cds_g_hi)
    cds_bp <- cds_hi_in - cds_lo_in + 1L
    left_bp  <- max(0L, cds_lo_in - g_lo)
    right_bp <- max(0L, g_hi - cds_hi_in)
    if (struct$strand == 1L) {
      ex_split[[i]] <- list(utr5_bp = left_bp,  cds_bp = cds_bp, utr3_bp = right_bp)
    } else {
      ex_split[[i]] <- list(utr5_bp = right_bp, cds_bp = cds_bp, utr3_bp = left_bp)
    }
  }
  utr5_lens <- vapply(ex_split, `[[`, integer(1), "utr5_bp")
  cds_lens  <- vapply(ex_split, `[[`, integer(1), "cds_bp")
  utr3_lens <- vapply(ex_split, `[[`, integer(1), "utr3_bp")

  utr5_total <- sum(utr5_lens)
  cds_total  <- sum(cds_lens)
  utr3_total <- sum(utr3_lens)

  # Optional visual-cap on UTR: rescale each exon's UTR contribution
  # proportionally so total UTR visual length = utr_visual_bp.
  if (!is.null(utr_visual_bp) && is.finite(utr_visual_bp) && utr_visual_bp > 0) {
    if (utr5_total > 0L) {
      scl <- as.numeric(utr_visual_bp) / utr5_total
      utr5_lens_vis <- utr5_lens * scl
    } else utr5_lens_vis <- utr5_lens
    if (utr3_total > 0L) {
      scl <- as.numeric(utr_visual_bp) / utr3_total
      utr3_lens_vis <- utr3_lens * scl
    } else utr3_lens_vis <- utr3_lens
  } else {
    utr5_lens_vis <- utr5_lens
    utr3_lens_vis <- utr3_lens
  }

  # Intron visual widths (one fewer than exons).
  intron_real_bp <- integer(max(0L, n_ex - 1L))
  if (n_ex > 1L) {
    for (i in seq_len(n_ex - 1L)) {
      # Real intron bp = gap between genomic neighbours, regardless of strand.
      a <- ex$g_start[i]; b <- ex$g_end[i]
      c2 <- ex$g_start[i + 1L]; d2 <- ex$g_end[i + 1L]
      gap_lo <- min(b, d2); gap_hi <- max(a, c2)
      intron_real_bp[i] <- max(0L, gap_hi - gap_lo - 1L)
    }
  }
  intron_scale <- match.arg(intron_scale, c("fixed", "proportional", "log"))
  intron_vis_bp <- switch(
    intron_scale,
    fixed        = rep(as.numeric(intron_visual_bp), length(intron_real_bp)),
    proportional = as.numeric(intron_real_bp),
    log          = {
      if (length(intron_real_bp) == 0L) numeric(0) else {
        ln <- log10(pmax(intron_real_bp, 1) + 1)
        scl <- if (mean(ln) > 0) mean(cds_lens[cds_lens > 0]) / mean(ln) else 1
        ln * scl
      }
    }
  )

  # Walk exons and lay out segments left-to-right in transcript 5'->3' order.
  regions <- list()
  cursor  <- 0.0
  anchors <- numeric(n_ex)
  for (i in seq_len(n_ex)) {
    anchors[i] <- cursor
    # Per-exon segments in 5'->3' order: utr5 -> cds -> utr3
    if (utr5_lens_vis[i] > 0) {
      regions[[length(regions) + 1L]] <- data.table::data.table(
        exon_idx = ex$exon_idx[i], region_kind = "utr5",
        g_start = ex$g_start[i],  g_end = ex$g_end[i],
        x_start = cursor, x_end = cursor + utr5_lens_vis[i])
      cursor <- cursor + utr5_lens_vis[i]
    }
    if (cds_lens[i] > 0) {
      regions[[length(regions) + 1L]] <- data.table::data.table(
        exon_idx = ex$exon_idx[i], region_kind = "cds",
        g_start = ex$g_start[i],  g_end = ex$g_end[i],
        x_start = cursor, x_end = cursor + cds_lens[i])
      cursor <- cursor + cds_lens[i]
    }
    if (utr3_lens_vis[i] > 0) {
      regions[[length(regions) + 1L]] <- data.table::data.table(
        exon_idx = ex$exon_idx[i], region_kind = "utr3",
        g_start = ex$g_start[i],  g_end = ex$g_end[i],
        x_start = cursor, x_end = cursor + utr3_lens_vis[i])
      cursor <- cursor + utr3_lens_vis[i]
    }
    # Intron (after this exon, except after the last one).
    if (i < n_ex) {
      regions[[length(regions) + 1L]] <- data.table::data.table(
        exon_idx = ex$exon_idx[i], region_kind = "intron",
        g_start = NA_integer_,    g_end = NA_integer_,
        x_start = cursor, x_end = cursor + intron_vis_bp[i])
      cursor <- cursor + intron_vis_bp[i]
    }
  }
  regions_dt <- data.table::rbindlist(regions)

  # X anchor for c.1 (first base of CDS): scan the first cds segment.
  cds_segs <- regions_dt[region_kind == "cds"]
  x_atg <- if (nrow(cds_segs) > 0L) cds_segs$x_start[1] else 0
  x_max <- cursor

  list(
    regions  = regions_dt,
    anchors  = anchors,
    cds_len  = as.integer(cds_total),
    utr5_len = as.integer(utr5_total),
    utr3_len = as.integer(utr3_total),
    x_atg    = x_atg,
    x_max    = x_max,
    intron_scale = intron_scale,
    intron_vis_bp = intron_vis_bp
  )
}

# Resolve variant_palette for this function: mirrors gvr_lollipop's resolver
# but stand-alone (so we don't depend on the closure-scoped one).
# ----- variant-class color table (package-internal) -----------------------
# Copy of the same table used inside gvr_lollipop(); kept here at package
# scope so .gvr_genepos_resolve_variant_palette can resolve symbolic names
# without invoking gvr_lollipop's local closure.
GVR_CLASS_COLORS <- c(
  "Translation_Start_Site" = "#000000", "Nonsense_Mutation" = "#D55E00",
  "Nonstop_Mutation"       = "#882255", "Splice_Site"       = "#CC79A7",
  "Frame_Shift_Del"        = "#E69F00", "Frame_Shift_Ins"   = "#F0E442",
  "In_Frame_Del"           = "#56B4E9", "In_Frame_Ins"      = "#0072B2",
  "Missense_Mutation"      = "#009E73", "Splice_Region"     = "#44AA99",
  "Protein_altering_variant" = "#117733", "Silent"          = "#999933",
  "5'UTR"                 = "#AA4499", "3'UTR"            = "#DDCC77",
  "5'Flank"               = "#88CCEE", "3'Flank"          = "#332288",
  "RNA"                    = "#BBBBBB", "Intron"            = "#DDDDDD",
  "IGR"                    = "#777777", "Targeted_Region"   = "#666666",
  "Other"                  = "#CCCCCC"
)

.gvr_genepos_resolve_variant_palette <- function(vp, classes_present) {
  classes_present <- as.character(classes_present)
  if (length(classes_present) == 0L) return(character(0))
  if (is.character(vp) && length(vp) == 1L && is.null(names(vp))) {
    if (identical(vp, "gvr")) {
      out <- GVR_CLASS_COLORS[classes_present]
      out[is.na(out)] <- GVR_CLASS_COLORS[["Other"]]
      names(out) <- classes_present
      return(out)
    }
    cols <- gvr_color_palette(vp, length(classes_present))
    return(stats::setNames(cols, classes_present))
  }
  if (!is.character(vp))
    stop("`variant_palette` must be a character string or character vector")
  nm <- names(vp); if (is.null(nm)) nm <- rep("", length(vp))
  overrides <- vp[nzchar(nm)]
  fallbacks <- vp[!nzchar(nm)]
  if (length(fallbacks) > 1L)
    stop("`variant_palette` accepts at most one unnamed element")
  fallback_name <- if (length(fallbacks) == 1L) unname(fallbacks) else "gvr"

  out <- if (identical(fallback_name, "gvr")) {
    tmp <- GVR_CLASS_COLORS[classes_present]
    tmp[is.na(tmp)] <- GVR_CLASS_COLORS[["Other"]]
    names(tmp) <- classes_present; tmp
  } else {
    stats::setNames(gvr_color_palette(fallback_name, length(classes_present)),
                    classes_present)
  }
  keep <- intersect(names(overrides), classes_present)
  if (length(keep)) out[keep] <- overrides[keep]
  out
}

# ---- Public: gvr_genepos.plot ----------------------------------------------

#' Gene-track lollipop plot on a cDNA axis
#'
#' Draws a per-gene track plot with exon, intron, and UTR segments for one
#' transcript, overlaid with lollipops placed on a cDNA-position x-axis using
#' the `HGVSc` field of each table row. Colours follow `Variant_Classification`
#' using the same palette set as [gvr_lollipop()].
#'
#' Companion to [gvr_lollipop()] which draws protein-domain rectangles on a
#' protein-position axis. This function instead draws gene structure on a
#' cDNA axis.
#'
#' @section Transcript resolution:
#' If `transcript_id` is `NULL`, the chosen transcript is, in order: the
#' first non-empty `MANE_SELECT` among table rows for `gene`; otherwise the
#' first `CANONICAL == "YES"` row; otherwise the transcript with the most
#' rows for that gene. The genome build is read from `NCBI_Build`.
#'
#' @section HGVSc parsing:
#' The leading `c.` integer of `HGVSc` is the variant anchor: positive
#' integers map to CDS positions, `-N` maps to the 5' UTR (negative axis),
#' `*N` maps to the 3' UTR (axis position `cds_len + N`). Strings with
#' `+K` or `-K` immediately after the leading anchor (splice / deep
#' intronic) are not plotted and counted under "splice"; anything else that
#' fails to parse is counted under "unparsed". The caption lists both
#' counts.
#'
#' @section Ensembl source:
#' By default the function calls `<host>/lookup/id/<ENST>?expand=1;utr=1`
#' on `rest.ensembl.org` (GRCh38) or `grch37.rest.ensembl.org` (GRCh37) and
#' caches the lean parsed structure as
#' `<cache_dir>/genestruct_ensembl_<ENST>_<assembly>.rds`. The cache
#' directory follows the same resolution chain as [gvr_lollipop()] (explicit
#' arg -> env `GVR_CACHE_DIR` -> option `germlinevaR.cache_dir` ->
#' `tools::R_user_dir("germlinevaR","cache")` -> `tempdir()`). Set
#' `gtf_path` to a Gencode/Ensembl GTF to skip REST entirely; this requires
#' suggesting `rtracklayer` or falls back to a tiny streaming parser.
#'
#' @param gvr A MAF-like data.table produced by [read.gvr()]. Required columns:
#'   `Hugo_Symbol`, `Transcript_ID`, `HGVSc`, `Variant_Classification`,
#'   `Tumor_Sample_Barcode`, `NCBI_Build`. `MANE_SELECT` and `CANONICAL`
#'   are used when `transcript_id` is auto-resolved.
#' @param gene Character(1). HGNC symbol (matches `Hugo_Symbol`).
#' @param transcript_id Character(1) or `NULL`. Ensembl stable id
#'   (`ENST...`). `NULL` triggers auto-resolution described above.
#' @param vc_keep Character vector or `NULL`. If non-`NULL`, only variants
#'   whose `Variant_Classification` is in `vc_keep` are kept. Otherwise the
#'   same default filter as [gvr_lollipop()] is applied (non-synonymous
#'   plus splice plus UTR).
#' @param color_by Character(1). One of `"vc"` (colour by
#'   `Variant_Classification`; default) or `"region"` (colour by track
#'   region: cds / utr5 / utr3 / splice / intron). Currently only `"vc"`
#'   has a stable visual contract.
#' @param intron_scale Character(1). One of `"fixed"` (default; constant
#'   visual width per intron, `intron_visual_bp` pixels), `"proportional"`
#'   (true bp), `"log"` (log10-scaled).
#' @param intron_visual_bp Integer(1). Visual width of every intron when
#'   `intron_scale == "fixed"`. Default `200`.
#' @param utr_visual_bp Integer(1) or `NULL`. When non-`NULL`, total UTR
#'   visual length is rescaled to this many bp regardless of real UTR
#'   length (useful for very long UTRs). `NULL` (default) means real
#'   length.
#' @param label_top Integer(1). Number of top-counted cDNA positions to
#'   label (matches [gvr_lollipop()]).
#' @param hotspot_window Numeric(1). Sliding-window width (bp on cDNA axis)
#'   for hotspot detection. Default `20`.
#' @param hotspot_min_n Numeric(1). Minimum distinct cDNA positions inside
#'   a window for it to be drawn as a hotspot band. Default `4`. Pass
#'   `Inf` to disable.
#' @param stem_alpha Numeric(1). Lollipop stem opacity. Default `0.6`.
#' @param point_size Numeric(1). Dot size. Default `3`.
#' @param bar_color,bar_border Character(1). Fill / border colour for the
#'   protein-body bar drawn over coding exons. Default `"grey85"` /
#'   `"grey40"`.
#' @param utr_color,intron_color,exon_color Character(1). Track region
#'   fill colours. Default `"grey70"` / `"grey60"` / `"steelblue"`.
#' @param variant_palette See [gvr_lollipop()] for the full grammar
#'   (`"gvr"`, palette name, or named override vector).
#' @param base_size,axis_text_size,axis_title_size Numeric. Text sizes.
#' @param ensembl_release Integer or `NULL`. Future use; currently unused
#'   (REST `latest` is always queried).
#' @param gtf_path Character(1) or `NULL`. Offline GTF override. If
#'   non-`NULL`, REST is not called.
#' @param cache_dir Character(1), `NULL`, or `FALSE`. Cache directory.
#'   `FALSE` disables on-disk caching.
#' @param out_dir Character(1) or `NULL`. Output directory root.
#' @param out_subdir Character(1). Subfolder under `out_dir`. Default
#'   `"gvr_genepos"`.
#' @param out_prefix Character(1) or `NULL`. File prefix; default
#'   `paste(gene, transcript_id, sep = "_")`.
#' @param format Character(1). One of `"png"`, `"svg"`, `"pdf"`, `"tiff"`.
#'   Default `"png"`.
#' @param width,height,dpi Numeric. Plot size and resolution.
#' @param verbose Logical(1). If `TRUE`, emit progress messages.
#'
#' @return A ggplot object, returned invisibly. As a side effect the plot
#'   is written to `<out_dir>/<out_subdir>/<out_prefix>.<format>` when
#'   `out_dir` is non-`NULL`.
#'
#' @examples
#' \dontrun{
#' # Auto-resolve MANE/CANONICAL transcript for BRCA1
#' p <- gvr_genepos.plot(gvr, "BRCA1")
#'
#' # Pin transcript and use proportional intron scaling
#' gvr_genepos.plot(gvr, "BRCA1",
#'                  transcript_id = "ENST00000357654",
#'                  intron_scale  = "proportional")
#'
#' # Fully offline using a local GTF
#' gvr_genepos.plot(gvr, "BRCA1",
#'                  gtf_path = "gencode.v44.annotation.gtf.gz")
#' }
#'
#' @export
gvr_genepos.plot <- function(gvr,
                             gene,
                             transcript_id    = NULL,
                             vc_keep          = NULL,
                             color_by         = c("vc", "region"),
                             intron_scale     = c("fixed", "proportional", "log"),
                             intron_visual_bp = 200L,
                             utr_visual_bp    = NULL,
                             label_top        = 5L,
                             hotspot_window   = 20L,
                             hotspot_min_n    = 4,
                             stem_alpha       = 0.6,
                             point_size       = 3,
                             bar_color        = "grey85",
                             bar_border       = "grey40",
                             utr_color        = "grey70",
                             intron_color     = "grey60",
                             exon_color       = "steelblue",
                             variant_palette  = "gvr",
                             base_size        = 12,
                             axis_text_size   = 11,
                             axis_title_size  = 12,
                             ensembl_release  = NULL,
                             gtf_path         = NULL,
                             cache_dir        = NULL,
                             out_dir          = ".",
                             out_subdir       = "gvr_genepos",
                             out_prefix       = NULL,
                             format           = c("png", "svg", "pdf", "tiff"),
                             width            = 10,
                             height           = 4,
                             dpi              = 300,
                             verbose          = TRUE) {

  color_by    <- match.arg(color_by)
  intron_scale <- match.arg(intron_scale)
  # `ensembl_release` is reserved for future use; the explicit reference
  # below silences an R CMD check 'parameter may not be used' NOTE without
  # changing behaviour.
  invisible(ensembl_release)
  format      <- match.arg(format)


  # ---- Argument validation -------------------------------------------------
  if (!data.table::is.data.table(gvr)) {
    if (is.data.frame(gvr)) gvr <- data.table::as.data.table(gvr)
    else stop("gvr_genepos.plot: 'gvr' must be a data.frame or data.table.")
  }
  required_cols <- c("Hugo_Symbol", "Transcript_ID", "HGVSc",
                     "Variant_Classification", "Tumor_Sample_Barcode",
                     "NCBI_Build")
  miss <- setdiff(required_cols, names(gvr))
  if (length(miss) > 0L)
    stop(sprintf("gvr_genepos.plot: 'gvr' missing required column(s): %s",
                 paste(miss, collapse = ", ")))
  if (!is.character(gene) || length(gene) != 1L || !nzchar(gene))
    stop("gvr_genepos.plot: 'gene' must be a single non-empty character.")

  # ---- Subset to this gene -------------------------------------------------
  mdt <- gvr[Hugo_Symbol == gene]
  if (nrow(mdt) == 0L)
    stop(sprintf("gvr_genepos.plot: gene '%s' not found in gvr$Hugo_Symbol.",
                 gene))

  # ---- Resolve transcript --------------------------------------------------
  if (is.null(transcript_id)) {
    # 1. MANE_SELECT
    enst_pick <- NA_character_
    if ("MANE_SELECT" %in% names(mdt)) {
      mane <- mdt[!is.na(MANE_SELECT) & nzchar(MANE_SELECT)]
      if (nrow(mane) > 0L) {
        # MANE_SELECT carries an ENST id; if it carries a NM_ id instead,
        # fall back to Transcript_ID for the same row.
        m1 <- mane$MANE_SELECT[1]
        if (grepl("^ENST", m1)) enst_pick <- sub("\\..*$", "", m1)
        else                    enst_pick <- mane$Transcript_ID[1]
      }
    }
    # 2. CANONICAL == "YES"
    if (is.na(enst_pick) && "CANONICAL" %in% names(mdt)) {
      can <- mdt[CANONICAL == "YES" & !is.na(Transcript_ID) & nzchar(Transcript_ID)]
      if (nrow(can) > 0L) enst_pick <- can$Transcript_ID[1]
    }
    # 3. Most-variants
    if (is.na(enst_pick)) {
      tt <- mdt[!is.na(Transcript_ID) & nzchar(Transcript_ID),
                .N, by = Transcript_ID]
      if (nrow(tt) > 0L) {
        data.table::setorder(tt, -N)
        enst_pick <- tt$Transcript_ID[1]
      }
    }
    if (is.na(enst_pick) || !nzchar(enst_pick))
      stop(sprintf("gvr_genepos.plot: could not resolve a transcript for '%s'.",
                   gene))
    # Strip version suffix if present (ENST00000357654.4 -> ENST00000357654)
    transcript_id <- sub("\\..*$", "", enst_pick)
    if (isTRUE(verbose))
      message(sprintf("gvr_genepos.plot: resolved transcript -> %s", transcript_id))
  }
  if (!grepl("^ENST", transcript_id))
    stop(sprintf("gvr_genepos.plot: transcript_id must be Ensembl ENST id; got '%s'.",
                 transcript_id))

  # ---- Resolve assembly ----------------------------------------------------
  builds <- unique(mdt$NCBI_Build[!is.na(mdt$NCBI_Build) & nzchar(mdt$NCBI_Build)])
  if (length(builds) == 0L) {
    assembly <- "GRCh38"
    if (isTRUE(verbose))
      message("gvr_genepos.plot: NCBI_Build empty for this gene; assuming GRCh38.")
  } else if (length(builds) > 1L) {
    stop(sprintf("gvr_genepos.plot: mixed NCBI_Build values for '%s': %s",
                 gene, paste(builds, collapse = ", ")))
  } else {
    assembly <- builds[1]
  }
  # Normalize common build aliases.
  assembly <- switch(toupper(assembly),
                     HG38 = "GRCh38", GRCH38 = "GRCh38",
                     HG19 = "GRCh37", GRCH37 = "GRCh37",
                     assembly)

  # ---- Fetch gene structure (REST or GTF) ----------------------------------
  cache_dir_resolved <- .gvr_genepos_resolve_cache_dir(cache_dir)
  if (isTRUE(verbose) && !is.na(cache_dir_resolved))
    message(sprintf("gvr_genepos.plot: cache_dir = %s", cache_dir_resolved))

  if (!is.null(gtf_path)) {
    struct <- .gvr_genestruct_from_gtf(gtf_path, transcript_id, assembly)
  } else {
    struct <- .gvr_genestruct_fetch(transcript_id, assembly, cache_dir_resolved, verbose)
  }

  # ---- Filter table rows for this transcript ---------------------------------
  # Keep rows whose Transcript_ID matches struct$transcript_id (post-versioning).
  mdt[, .__enst__ := sub("\\..*$", "", as.character(Transcript_ID))]
  mdt <- mdt[.__enst__ == transcript_id]
  if (nrow(mdt) == 0L)
    warning(sprintf("gvr_genepos.plot: no rows match transcript %s after filtering.",
                    transcript_id))

  # ---- Variant-class filter (default mirrors gvr_lollipop) -----------------
  if (is.null(vc_keep)) {
    vc_keep <- c("Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
                 "Translation_Start_Site", "Splice_Site", "Splice_Region",
                 "In_Frame_Del", "In_Frame_Ins",
                 "Frame_Shift_Del", "Frame_Shift_Ins",
                 "5'UTR", "3'UTR", "5_prime_UTR_variant", "3_prime_UTR_variant",
                 "5'Flank", "3'Flank")
  }
  mdt <- mdt[Variant_Classification %in% vc_keep]


  # ---- Parse HGVSc to (kind, pos) ------------------------------------------
  parsed <- .gvr_parse_hgvsc(mdt$HGVSc)
  mdt[, .__kind__ := parsed$kind]
  mdt[, .__pos__  := parsed$pos]
  mdt[, .__valid__ := parsed$valid]

  # Footnote counts
  # n_total intentionally unused (footnote only counts not-placed rows)
  n_splice   <- sum(mdt$.__kind__ == "splice", na.rm = TRUE)
  n_unparsed <- sum(is.na(mdt$.__kind__) | (!mdt$.__valid__ & mdt$.__kind__ != "splice"),
                    na.rm = TRUE)

  # Keep only placeable rows.
  dt <- mdt[.__valid__ == TRUE & .__kind__ %in% c("cds", "utr5", "utr3")]

  # ---- Build layout --------------------------------------------------------
  lay <- .gvr_genepos_layout(struct,
                             intron_scale     = intron_scale,
                             intron_visual_bp = intron_visual_bp,
                             utr_visual_bp    = utr_visual_bp)
  regions    <- lay$regions
  cds_len    <- lay$cds_len
  utr5_len   <- lay$utr5_len
  utr3_len   <- lay$utr3_len
  x_atg      <- lay$x_atg

  # ---- Map each variant to its x coordinate --------------------------------
  # For CDS variants, x is the visual location of CDS position pos. We need
  # the running CDS offset across exons. Walk CDS segments in order; for each
  # variant, find the segment whose CDS-position interval contains pos.
  cds_segs <- regions[region_kind == "cds"]
  data.table::setorder(cds_segs, x_start)
  cds_seg_len <- cds_segs$x_end - cds_segs$x_start
  cds_seg_cum_end <- cumsum(cds_seg_len)   # CDS-positions covered by segment i
  cds_seg_cum_start <- c(0, cds_seg_cum_end[-length(cds_seg_cum_end)])

  utr5_segs <- regions[region_kind == "utr5"]
  data.table::setorder(utr5_segs, x_start)
  utr5_seg_len <- utr5_segs$x_end - utr5_segs$x_start
  utr5_seg_cum_end <- cumsum(utr5_seg_len)
  utr5_seg_cum_start <- c(0, utr5_seg_cum_end[-length(utr5_seg_cum_end)])

  utr3_segs <- regions[region_kind == "utr3"]
  data.table::setorder(utr3_segs, x_start)
  utr3_seg_len <- utr3_segs$x_end - utr3_segs$x_start
  utr3_seg_cum_end <- cumsum(utr3_seg_len)
  utr3_seg_cum_start <- c(0, utr3_seg_cum_end[-length(utr3_seg_cum_end)])

  map_cds_pos_to_x <- function(p) {
    if (is.na(p) || p < 1L || p > cds_len) return(NA_real_)
    # Find segment i such that cds_seg_cum_start[i] < p <= cds_seg_cum_end[i]
    i <- which(p <= cds_seg_cum_end)[1]
    if (is.na(i)) return(NA_real_)
    offset_in_seg <- p - cds_seg_cum_start[i]
    cds_segs$x_start[i] + offset_in_seg
  }
  map_utr5_pos_to_x <- function(p) {
    # p is negative; |p| is distance upstream of c.1 along 5UTR.
    ap <- abs(p)
    if (is.na(ap) || ap < 1L || ap > utr5_len) return(NA_real_)
    # 5UTR segments are laid 5'->3' on the axis; segment 1 holds the most
    # upstream UTR positions. Position p = -ap corresponds to offset
    # (utr5_len - ap + 1) from start of segment 1.
    target <- utr5_len - ap + 1L
    i <- which(target <= utr5_seg_cum_end)[1]
    if (is.na(i)) return(NA_real_)
    offset_in_seg <- target - utr5_seg_cum_start[i]
    utr5_segs$x_start[i] + offset_in_seg
  }
  map_utr3_pos_to_x <- function(p) {
    if (is.na(p) || p < 1L || p > utr3_len) return(NA_real_)
    i <- which(p <= utr3_seg_cum_end)[1]
    if (is.na(i)) return(NA_real_)
    offset_in_seg <- p - utr3_seg_cum_start[i]
    utr3_segs$x_start[i] + offset_in_seg
  }

  dt[, .__x__ := NA_real_]
  for (kk in c("cds", "utr5", "utr3")) {
    idx <- which(dt$.__kind__ == kk)
    if (length(idx) == 0L) next
    mapper <- switch(kk,
                     cds  = map_cds_pos_to_x,
                     utr5 = map_utr5_pos_to_x,
                     utr3 = map_utr3_pos_to_x)
    dt$.__x__[idx] <- vapply(dt$.__pos__[idx], mapper, numeric(1))
  }
  # Drop any row whose mapping failed (out-of-bounds).
  n_outofbounds <- sum(is.na(dt$.__x__))
  if (n_outofbounds > 0L) {
    dt <- dt[!is.na(.__x__)]
    n_unparsed <- n_unparsed + n_outofbounds
  }

  # ---- Stack dots at the same x position -----------------------------------
  data.table::setorder(dt, .__x__, Variant_Classification, Tumor_Sample_Barcode)
  dt[, .__y__ := seq_len(.N), by = .__x__]
  pos_height <- dt[, list(.__top__ = max(.__y__)), by = .__x__]

  # ---- Top-N position labels -----------------------------------------------
  label_n <- if (is.infinite(label_top)) Inf else as.integer(label_top)
  label_df <- NULL
  if (label_n > 0L && nrow(dt) > 0L) {
    pos_counts <- dt[, list(.n = .N, kind = .__kind__[1], pos = .__pos__[1]),
                     by = .__x__]
    pos_counts <- pos_counts[order(-.n)]
    keep_n <- if (is.infinite(label_n)) nrow(pos_counts) else min(label_n, nrow(pos_counts))
    top_x  <- pos_counts$.__x__[seq_len(keep_n)]
    label_rows <- vector("list", length(top_x))
    for (i in seq_along(top_x)) {
      xv <- top_x[i]
      sub_dt <- dt[.__x__ == xv]
      # Always label with HGVSc. Pick the most-frequent HGVSc string at this
      # position; fall back to reconstructed c.<pos> notation only if HGVSc
      # is absent or blank for every row at this position.
      pretty <- if ("HGVSc" %in% names(sub_dt) &&
                    any(!is.na(sub_dt$HGVSc) & nzchar(sub_dt$HGVSc))) {
        # Strip transcript prefix e.g. "ENST00000269305.9:c.215C>G" -> "c.215C>G"
        raw <- sub_dt$HGVSc[!is.na(sub_dt$HGVSc) & nzchar(sub_dt$HGVSc)]
        raw <- sub("^[^:]+:", "", raw)
        tab <- table(raw)
        names(sort(tab, decreasing = TRUE))[1]
      } else {
        k <- sub_dt$.__kind__[1]; p <- sub_dt$.__pos__[1]
        switch(k, cds  = sprintf("c.%d",   p),
                  utr5 = sprintf("c.-%d",  abs(p)),
                  utr3 = sprintf("c.*%d",  p),
                  NA_character_)
      }
      raw_all <- if ("HGVSc" %in% names(sub_dt))
        sub("^[^:]+:", "", sub_dt$HGVSc[!is.na(sub_dt$HGVSc) & nzchar(sub_dt$HGVSc)])
      else character(0)
      n_distinct <- length(unique(raw_all))
      lbl <- if (n_distinct > 1L) sprintf("%s (+%d more)", pretty, n_distinct - 1L) else pretty
      label_rows[[i]] <- data.frame(
        x     = xv,
        top   = pos_height$.__top__[pos_height$.__x__ == xv],
        label = lbl,
        stringsAsFactors = FALSE
      )
    }
    label_df <- do.call(rbind, label_rows)
  }

  # ---- Variant-class palette -----------------------------------------------
  classes_present <- sort(unique(as.character(dt$Variant_Classification)))
  col_map <- .gvr_genepos_resolve_variant_palette(variant_palette, classes_present)


  # ---- Geometry constants (match gvr_lollipop) -----------------------------
  bar_half <- 0.30
  dom_half <- 0.60          # exon/UTR box half-height
  utr_half <- 0.30          # UTR box half-height (half of CDS exon height)
  intron_y <- 0
  hot_ymin <- -dom_half
  y_max_stack <- if (nrow(pos_height) > 0L) max(pos_height$.__top__) else 1
  y_upper <- y_max_stack + max(1, y_max_stack * 0.10)
  y_lower <- -dom_half - 0.20

  # ---- Hotspot detection (mirrors gvr_lollipop logic; on cDNA axis) --------
  .hw <- suppressWarnings(as.numeric(hotspot_window))
  .hm <- suppressWarnings(as.numeric(hotspot_min_n))
  .hotspots_enabled <- isTRUE(is.finite(.hw) && .hw > 0) &&
                       isTRUE(!is.na(.hm) && .hm > 0)
  hotspot_df <- data.frame(xmin = numeric(0), xmax = numeric(0),
                           xmid = numeric(0), n_in_window = integer(0),
                           stringsAsFactors = FALSE)
  if (.hotspots_enabled && is.finite(.hm) && nrow(dt) > 0L) {
    # Use only CDS-domain positions for hotspot detection (consistent with
    # the gvr_lollipop convention -- AA positions only).
    cds_x <- sort(unique(dt$.__x__[dt$.__kind__ == "cds"]))
    if (length(cds_x) >= .hm) {
      .half_w <- .hw / 2
      .in_win <- vapply(cds_x, function(p) {
        sum(cds_x >= (p - .half_w) & cds_x <= (p + .half_w))
      }, integer(1))
      .seed_pos <- cds_x[.in_win >= .hm]
      if (length(.seed_pos) > 0L) {
        .cand <- data.frame(
          xmin        = pmax(0, .seed_pos - .half_w),
          xmax        = .seed_pos + .half_w,
          n_in_window = .in_win[.in_win >= .hm],
          stringsAsFactors = FALSE
        )
        .cand <- .cand[order(.cand$xmin), , drop = FALSE]
        .merged <- vector("list", nrow(.cand))
        .cur <- .cand[1, , drop = FALSE]; .k <- 1L
        if (nrow(.cand) > 1L) {
          for (.i in 2:nrow(.cand)) {
            if (.cand$xmin[.i] <= .cur$xmax) {
              .cur$xmax        <- max(.cur$xmax, .cand$xmax[.i])
              .cur$n_in_window <- max(.cur$n_in_window, .cand$n_in_window[.i])
            } else {
              .merged[[.k]] <- .cur; .k <- .k + 1L
              .cur <- .cand[.i, , drop = FALSE]
            }
          }
        }
        .merged[[.k]] <- .cur
        hotspot_df <- do.call(rbind, .merged[seq_len(.k)])
        hotspot_df$xmid <- (hotspot_df$xmin + hotspot_df$xmax) / 2
        hotspot_df <- hotspot_df[, c("xmin","xmax","xmid","n_in_window"), drop = FALSE]
      }
    }
  }

  # ---- Re-anchor x-axis so c.1 visually sits at zero -----------------------
  # We translate every visual x so that x_atg -> 0; 5UTR ends up negative,
  # CDS sits in [0, cds_len], 3UTR sits in [cds_len, cds_len + utr3_len_vis].
  regions[, x_start := x_start - x_atg]
  regions[, x_end   := x_end   - x_atg]
  dt[, .__x__ := .__x__ - x_atg]
  pos_height[, .__x__ := .__x__ - x_atg]
  if (!is.null(label_df)) label_df$x <- label_df$x - x_atg
  if (nrow(hotspot_df) > 0L) {
    hotspot_df$xmin <- hotspot_df$xmin - x_atg
    hotspot_df$xmax <- hotspot_df$xmax - x_atg
    hotspot_df$xmid <- hotspot_df$xmid - x_atg
  }

  # ---- Track rectangles + intron lines + per-region fills ------------------
  rect_df <- regions[region_kind %in% c("cds", "utr5", "utr3")]
  rect_df[, ymin := data.table::fcase(
    region_kind == "cds",                    -dom_half,
    region_kind %in% c("utr5", "utr3"),  -utr_half
  )]
  rect_df[, ymax := data.table::fcase(
    region_kind == "cds",                     dom_half,
    region_kind %in% c("utr5", "utr3"),   utr_half
  )]
  rect_df[, fill := data.table::fcase(
    region_kind == "cds",  exon_color,
    region_kind == "utr5", utr_color,
    region_kind == "utr3", utr_color
  )]

  intron_df <- regions[region_kind == "intron"]
  intron_df[, y := intron_y]

  # ---- ggplot composition (drawing order matches gvr_lollipop) -------------
  p <- ggplot2::ggplot() +
    # 1. hotspot wash bands (behind everything except panel grid)
    {
      if (nrow(hotspot_df) > 0L) {
        ggplot2::geom_rect(data = hotspot_df,
                           ggplot2::aes(xmin = xmin, xmax = xmax,
                                        ymin = hot_ymin, ymax = y_upper),
                           inherit.aes = FALSE,
                           fill = "#ffeeba", alpha = 0.35)
      } else NULL
    } +
    # 2a. CDS exon rectangles (with border)
    ggplot2::geom_rect(
      data = rect_df[rect_df$region_kind == "cds", ],
      ggplot2::aes(xmin = x_start, xmax = x_end, ymin = ymin, ymax = ymax),
      fill = exon_color, color = bar_border, linewidth = 0.3,
      inherit.aes = FALSE) +
    # 2b. UTR rectangles (no border)
    ggplot2::geom_rect(
      data = rect_df[rect_df$region_kind %in% c("utr5", "utr3"), ],
      ggplot2::aes(xmin = x_start, xmax = x_end, ymin = ymin, ymax = ymax),
      fill = utr_color, color = NA,
      inherit.aes = FALSE) +
    # 3. intron lines connecting consecutive exons
    {
      if (nrow(intron_df) > 0L) {
        ggplot2::geom_segment(data = intron_df,
                              ggplot2::aes(x = x_start, xend = x_end,
                                           y = y, yend = y),
                              color = intron_color, linewidth = 0.4,
                              inherit.aes = FALSE)
      } else NULL
    } +
    # 5. lollipop stems
    {
      if (nrow(pos_height) > 0L) {
        stem_df <- data.frame(x = pos_height$.__x__, top = pos_height$.__top__)
        ggplot2::geom_segment(data = stem_df,
                              ggplot2::aes(x = x, xend = x, y = bar_half, yend = top),
                              alpha = stem_alpha, color = "grey40", linewidth = 0.4,
                              inherit.aes = FALSE)
      } else NULL
    } +
    # 6. dots
    {
      if (nrow(dt) > 0L) {
        dot_df <- data.frame(x = dt$.__x__, y = dt$.__y__,
                             vc = as.character(dt$Variant_Classification),
                             stringsAsFactors = FALSE)
        ggplot2::geom_point(data = dot_df,
                            ggplot2::aes(x = x, y = y, color = vc),
                            size = point_size, inherit.aes = FALSE)
      } else NULL
    } +
    ggplot2::scale_color_manual(values = col_map, name = "Variant class") +
    ggplot2::scale_x_continuous(
      breaks = function(lim) {
        # Same `nice_step` approach as gvr_lollipop: aim for ~10 ticks
        L <- diff(lim)
        if (!is.finite(L) || L <= 0) return(base::pretty(lim, n = 10))
        raw <- L / 10; pw <- 10 ^ floor(log10(raw)); b <- raw / pw
        step <- (if (b < 1.5) 1 else if (b < 3) 2 else if (b < 7) 5 else 10) * pw
        seq(floor(lim[1] / step) * step, ceiling(lim[2] / step) * step, by = step)
      },
      expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::scale_y_continuous(
      limits = c(y_lower, y_upper),
      breaks = function(lim) {
        ub <- max(1, ceiling(lim[2]))
        seq(0, ub, by = max(1, ceiling(ub / 5)))
      },
      expand = c(0, 0))


  # ---- Title, subtitle, caption, theme -------------------------------------
  n_samples <- length(unique(dt$Tumor_Sample_Barcode))
  cap_parts <- character(0)
  if (n_splice > 0L)   cap_parts <- c(cap_parts, sprintf("%d splice/intronic", n_splice))
  if (n_unparsed > 0L) cap_parts <- c(cap_parts, sprintf("%d unparsed/out-of-bounds", n_unparsed))
  caption_txt <- if (length(cap_parts) > 0L) {
    sprintf("%d variant(s) not placed: %s",
            n_splice + n_unparsed, paste(cap_parts, collapse = ", "))
  } else NULL

  p <- p +
    ggplot2::labs(
      title    = sprintf("%s -- %s (%s)", gene, transcript_id, assembly),
      subtitle = sprintf("cDNA length: %d nt | %d coding exon segment(s) | %d variant(s) in %d sample(s)",
                         struct$cdna_len,
                         sum(regions$region_kind == "cds"),
                         nrow(dt), n_samples),
      x        = "cDNA position (c.)",
      y        = "Number of sample-variants",
      color    = "Variant class",
      caption  = caption_txt) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(face = "bold",
                                                 size = base_size * 1.25,
                                                 hjust = 0,
                                                 margin = ggplot2::margin(b = 2)),
      plot.subtitle      = ggplot2::element_text(size = base_size * 0.85,
                                                 color = "grey35",
                                                 hjust = 0,
                                                 margin = ggplot2::margin(b = 8)),
      plot.caption       = ggplot2::element_text(size = base_size * 0.75,
                                                 color = "grey45",
                                                 hjust = 1,
                                                 margin = ggplot2::margin(t = 4)),
      axis.title         = ggplot2::element_text(size = axis_title_size),
      axis.text          = ggplot2::element_text(size = axis_text_size),
      legend.title       = ggplot2::element_text(face = "bold", size = base_size * 0.85),
      legend.text        = ggplot2::element_text(size = base_size * 0.80),
      legend.key.size    = ggplot2::unit(0.7, "lines"),
      legend.position    = "right",
      legend.justification = c(0, 0.5),
      legend.box.spacing = ggplot2::unit(0.3, "lines"),
      plot.margin        = ggplot2::margin(t = 8, r = 12, b = 6, l = 6)
    )

  # ---- Top-N labels --------------------------------------------------------
  if (!is.null(label_df) && nrow(label_df) > 0L) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(data = label_df,
                                        ggplot2::aes(x = x, y = top, label = label),
                                        size = base_size * 0.3,
                                        nudge_y = 0.5,
                                        segment.size = 0.2,
                                        max.overlaps = Inf,
                                        seed = 42L,
                                        inherit.aes = FALSE)
    } else {
      p <- p + ggplot2::geom_text(data = label_df,
                                  ggplot2::aes(x = x, y = top, label = label),
                                  size = base_size * 0.3,
                                  vjust = -0.6, hjust = 0.5,
                                  inherit.aes = FALSE)
    }
  }

  # ---- File output ---------------------------------------------------------
  if (!is.null(out_dir)) {
    target_dir <- if (!is.null(out_subdir) && nzchar(out_subdir))
      file.path(out_dir, out_subdir) else out_dir
    if (!dir.exists(target_dir))
      dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
    if (is.null(out_prefix) || !nzchar(out_prefix))
      out_prefix <- paste(gene, transcript_id, sep = "_")
    out_path <- file.path(target_dir, sprintf("%s.%s", out_prefix, format))

    if (file.exists(out_path) && isTRUE(verbose))
      message(sprintf("gvr_genepos.plot: overwriting existing %s", out_path))

    saver <- function(tmp) {
      dev_arg <- switch(format,
                        png  = "png",
                        svg  = "svg",
                        pdf  = "pdf",
                        tiff = "tiff")
      args <- list(filename = tmp, plot = p,
                   width = width, height = height,
                   units = "in", dpi = dpi,
                   device = dev_arg)
      if (identical(format, "svg") &&
          !requireNamespace("svglite", quietly = TRUE)) {
        # Fall back to base grDevices::svg
        grDevices::svg(tmp, width = width, height = height)
        print(p); grDevices::dev.off()
        return(invisible(NULL))
      }
      if (identical(format, "tiff")) args$compression <- "lzw"
      do.call(ggplot2::ggsave, args)
    }
    written <- .gvr_genepos_fuse_save(out_path, saver)
    if (isTRUE(verbose) && !is.na(written))
      message(sprintf("gvr_genepos.plot: wrote %s", written))
  }

  invisible(p)
}

