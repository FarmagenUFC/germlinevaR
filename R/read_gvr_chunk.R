## ---------------------------------------------------------------------------
## read_gvr_chunk.R
##
## Package-internal helpers that decompose the per-chunk parsing pipeline used
## inside `read.gvr()`'s `convert_chunk()` (see R/read.gvr.R). These were lifted
## verbatim from `convert_chunk()`'s body in Turn-5b of the Bioconductor
## submission refactor. No behavioural changes vs the prior version --- the
## PERF tricks (O1/O2/O6/O7/O8/A3/P1/P2) and `# why:` comments are preserved
## byte-identical so the numeric-equality gate and the hero-PDF SHA256 gate
## both stay green.
##
## All helpers below are PURE: they take state as explicit arguments, return
## an updated state object, and never mutate parent-frame bindings. The per-
## chunk caches (`.gvr_mstr_cache`, `.gvr_v2m_cache`) ARE passed in as explicit
## arguments to `.gvr_chunk_build_record()` so the helpers do not capture them
## from the enclosing closure.
##
## These helpers are intentionally NOT exported (no `@export`, `@noRd`,
## `@keywords internal`) -- they are an implementation detail of the chunked
## VCF reader.
## ---------------------------------------------------------------------------


## ---- Helper 1 -------------------------------------------------------------
#' Lift VCF columns and CSQ-field positions out of a chunk data.table.
#'
#' Internal helper for `read.gvr()`'s `convert_chunk()`. Hoists the 10
#' standard VCF columns to plain vectors (C3 PERF trick: kills per-record
#' `$`/is.factor overhead inside the row loop) and pre-splits the ALT
#' column once for the whole chunk (O2 PERF trick: ~30x faster on
#' multi-ALT edge cases). Also caches the `match()` positions of the 9
#' CSQ subfields the row loop needs.
#'
#' @param dt A data.table chunk from the on-the-fly VCF reader.
#' @param csq_fields Character vector of CSQ subfield names parsed from
#'   the VCF header's CSQ INFO description.
#'
#' @return A named list with the chunk-level vectors (`CHROM_v`...`SAMPLE_v`,
#'   `nrow_dt`, `alts_all`) plus `n_csq` and a list `P` of cached CSQ
#'   subfield positions (`P_Allele`, `P_Cons`, `P_SYMBOL`, `P_Feature`,
#'   `P_HGVSc`, `P_HGVSp`, `P_Existing`, `P_CANONICAL`, `P_MANESEL`).
#'
#' @keywords internal
#' @noRd
.gvr_chunk_setup_columns <- function(dt, csq_fields) {
    n_csq <- length(csq_fields)
    ci <- function(name) match(name, csq_fields)
    P_Allele <- ci("Allele")
    P_Cons <- ci("Consequence")
    P_SYMBOL <- ci("SYMBOL")
    P_Feature <- ci("Feature")
    P_HGVSc <- ci("HGVSc")
    P_HGVSp <- ci("HGVSp")
    P_Existing <- ci("Existing_variation")
    P_CANONICAL <- ci("CANONICAL")
    P_MANESEL <- ci("MANE_SELECT")

    # C3: hoist columns to plain vectors once (kills per-record `$`/is.factor overhead)
    CHROM_v <- dt$CHROM
    POS_v <- dt$POS
    ID_v <- dt$ID
    REF_v <- dt$REF
    ALT_v <- dt$ALT
    QUAL_v <- dt$QUAL
    FILTER_v <- dt$FILTER
    INFO_v <- dt$INFO
    FORMAT_v <- dt$FORMAT
    SAMPLE_v <- dt$SAMPLE
    nrow_dt <- nrow(dt)

    # O2 (PERF): split the ALT column ONCE for the whole chunk instead of calling
    # strsplit(altf, ",") per record. base strsplit() over a character vector returns
    # a list-of-vectors; alts_all[[r]] is byte-identical to the old per-record
    # strsplit(ALT_v[r], ",", fixed=TRUE)[[1]]. ~30x faster on the isolated split
    # (multi-ALT is rare, so per-record calls were almost pure call-overhead).
    alts_all <- strsplit(ALT_v, ",", fixed = TRUE)

    list(
        CHROM_v = CHROM_v, POS_v = POS_v, ID_v = ID_v, REF_v = REF_v, ALT_v = ALT_v,
        QUAL_v = QUAL_v, FILTER_v = FILTER_v, INFO_v = INFO_v,
        FORMAT_v = FORMAT_v, SAMPLE_v = SAMPLE_v,
        nrow_dt = nrow_dt, alts_all = alts_all, n_csq = n_csq,
        P = list(
            P_Allele = P_Allele, P_Cons = P_Cons, P_SYMBOL = P_SYMBOL,
            P_Feature = P_Feature, P_HGVSc = P_HGVSc, P_HGVSp = P_HGVSp,
            P_Existing = P_Existing, P_CANONICAL = P_CANONICAL,
            P_MANESEL = P_MANESEL
        )
    )
}


## ---- Helper 2 -------------------------------------------------------------
#' Resolve FORMAT/SAMPLE into chunk-level GT/AD/DP/GQ vectors when constant.
#'
#' Internal helper for `read.gvr()`'s `convert_chunk()`. Implements the O1
#' PERF trick: when the FORMAT string is constant across the chunk (the
#' usual case for VEP/GATK germline VCFs, e.g. `"GT:AD:DP:GQ:PL"`), the
#' GT/AD/DP/GQ field POSITIONS are resolved once and the values pulled with
#' a single vectorised `data.table::tstrsplit(SAMPLE, ":")` over the whole
#' column (~4x faster than per-record splits). If FORMAT varies anywhere
#' in the chunk, the caller falls back to the per-record `strsplit()` +
#' name-based lookup path (the original Turn-3 path), so output is
#' byte-identical in every case (verified by the numeric-equality gate).
#'
#' @param FORMAT_v Character vector of FORMAT cells for the chunk.
#' @param SAMPLE_v Character vector of SAMPLE cells for the chunk.
#'
#' @return Named list with `fmt_constant` (logical), the four pre-extracted
#'   GT/AD/DP/GQ columns (or NULL when missing), and the matched positions
#'   `pos_GT/pos_AD/pos_DP/pos_GQ`.
#'
#' @keywords internal
#' @noRd
.gvr_chunk_resolve_format <- function(FORMAT_v, SAMPLE_v) {
    # O1 (PERF): the per-record loop split FORMAT and SAMPLE on ":" every iteration and
    # then looked GT/AD/DP/GQ up BY NAME. In real VEP/GATK germline VCFs the FORMAT
    # string is almost always constant within a chunk (e.g. "GT:AD:DP:GQ:PL"). When it
    # IS constant we resolve the GT/AD/DP/GQ field POSITIONS once and pull them with a
    # single vectorized tstrsplit(SAMPLE, ":") over the whole column (~4x faster). If
    # FORMAT varies anywhere in the chunk we fall back to the exact per-record path, so
    # output is byte-identical in every case (verified by the validation gate).
    fmt_u <- unique(FORMAT_v)
    fmt_constant <- length(fmt_u) == 1L && !is.na(fmt_u[1L])
    GT_col <- AD_col <- DP_col <- GQ_col <- NULL
    pos_GT <- pos_AD <- pos_DP <- pos_GQ <- NA_integer_
    if (fmt_constant) {
        fmt_keys0 <- strsplit(fmt_u[1L], ":", fixed = TRUE)[[1]]
        pos_GT <- match("GT", fmt_keys0)
        pos_AD <- match("AD", fmt_keys0)
        pos_DP <- match("DP", fmt_keys0)
        pos_GQ <- match("GQ", fmt_keys0)
        smp_split <- data.table::tstrsplit(SAMPLE_v, ":", fixed = TRUE)
        n_smp_fields <- length(smp_split)
        # Pre-extract per-record GT/AD/DP/GQ vectors. A position is only valid if it
        # exists in FORMAT (pos_* not NA) AND the SAMPLE column actually carried that
        # many fields. tstrsplit pads short rows with NA, matching the by-name lookup
        # (a field absent from a given SAMPLE yielded NA / the default there too).
        GT_col <- if (!is.na(pos_GT) && pos_GT <= n_smp_fields) smp_split[[pos_GT]] else NULL
        AD_col <- if (!is.na(pos_AD) && pos_AD <= n_smp_fields) smp_split[[pos_AD]] else NULL
        DP_col <- if (!is.na(pos_DP) && pos_DP <= n_smp_fields) smp_split[[pos_DP]] else NULL
        GQ_col <- if (!is.na(pos_GQ) && pos_GQ <= n_smp_fields) smp_split[[pos_GQ]] else NULL
    }
    list(
        fmt_constant = fmt_constant,
        GT_col = GT_col, AD_col = AD_col, DP_col = DP_col, GQ_col = GQ_col,
        pos_GT = pos_GT, pos_AD = pos_AD, pos_DP = pos_DP, pos_GQ = pos_GQ
    )
}


## ---- Helper 3 -------------------------------------------------------------
#' Vectorised DP/GQ pre-filter (O6 PERF trick).
#'
#' Internal helper for `read.gvr()`'s `convert_chunk()`. Applies the O6
#' vectorised DP/GQ pre-filter BEFORE the heavy per-record CSQ loop. When
#' DP/GQ filtering is active, ~45% of records typically fail; subsetting
#' chunk-level vectors here lets the loop skip them entirely. The DP/GQ
#' values used are the SAME ones the original loop would extract (from
#' the O1 pre-split when FORMAT is constant), so the filter result is
#' byte-identical. Missing/non-numeric DP or GQ never causes a drop
#' (treated as "pass"), matching the original per-record behaviour.
#'
#' @param state Chunk state list produced by `.gvr_chunk_setup_columns()`
#'   merged with the output of `.gvr_chunk_resolve_format()` (so it carries
#'   both the chunk column vectors and the O1 fast-path fields).
#' @param filter_dp,filter_gq Logical scalars (the `read.gvr()` args).
#' @param min_DP,min_GQ Numeric thresholds (the `read.gvr()` args).
#'
#' @return List `(state = <updated state list>, n_dropped_dpgq = <int>)`.
#'   When no rows are dropped, the state is returned unchanged.
#'
#' @keywords internal
#' @noRd
.gvr_chunk_filter_dpgq <- function(state, filter_dp, filter_gq, min_DP, min_GQ) {
    # --- O6 (PERF): vectorized DP/GQ pre-filter BEFORE the heavy loop ---------------
    # When DP/GQ filtering is active, ~45% of records fail and are discarded. Doing
    # this check INSIDE the loop means all the expensive CSQ work (strsplit, block
    # assignment, consequence ranking) runs on records that are immediately dropped.
    # By vectorizing the check over the whole chunk and subsetting BEFORE the loop,
    # we skip ~45% of the heaviest work. The DP/GQ values used here are the SAME ones
    # the loop would extract (from the O1 pre-split), so the filter is byte-identical.
    # Missing/non-numeric DP or GQ never causes a drop (treated as "pass"), matching
    # the original per-record behaviour.
    n_dropped_dpgq <- 0L
    if (!(filter_dp || filter_gq)) {
        return(list(state = state, n_dropped_dpgq = 0L))
    }

    nrow_dt <- state$nrow_dt
    fmt_constant <- state$fmt_constant
    keep <- rep(TRUE, nrow_dt)

    if (fmt_constant) {
        # Fast path: DP/GQ already extracted as chunk-level vectors
        if (filter_dp && !is.null(state$DP_col)) {
            # why: as.numeric() on the DP field which may contain '' for missing; NAs are passed through the keep filter via is.na(dp_num)|dp_num>min_DP.
            dp_num <- suppressWarnings(as.numeric(state$DP_col))
            keep <- keep & (is.na(dp_num) | dp_num > min_DP)
        }
        if (filter_gq && !is.null(state$GQ_col)) {
            # why: as.numeric() on the GQ field which may contain ''; NAs propagated by the same is.na()|>thr pattern.
            gq_num <- suppressWarnings(as.numeric(state$GQ_col))
            keep <- keep & (is.na(gq_num) | gq_num > min_GQ)
        }
    } else {
        # Slow path: FORMAT varies -- extract DP/GQ per record (same logic as the loop)
        FORMAT_v <- state$FORMAT_v
        SAMPLE_v <- state$SAMPLE_v
        for (r in seq_len(nrow_dt)) {
            fmt_keys <- strsplit(FORMAT_v[r], ":", fixed = TRUE)[[1]]
            smp_vals <- strsplit(SAMPLE_v[r], ":", fixed = TRUE)[[1]]
            names(smp_vals) <- fmt_keys[seq_along(smp_vals)]
            if (filter_dp) {
                sdp <- if ("DP" %in% names(smp_vals)) smp_vals[["DP"]] else NA_character_
                # why: as.numeric() on per-sample DP when it is read row-by-row in fallback mode; non-numeric becomes NA and the row is kept.
                dp_num <- suppressWarnings(as.numeric(sdp))
                if (!is.na(dp_num) && dp_num <= min_DP) keep[r] <- FALSE
            }
            if (filter_gq) {
                sgq <- if ("GQ" %in% names(smp_vals)) smp_vals[["GQ"]] else NA_character_
                # why: as.numeric() on per-sample GQ in fallback mode; same NA-keep semantics as DP.
                gq_num <- suppressWarnings(as.numeric(sgq))
                if (!is.na(gq_num) && gq_num <= min_GQ) keep[r] <- FALSE
            }
        }
    }
    n_dropped_dpgq <- sum(!keep)
    if (n_dropped_dpgq > 0L) {
        # Subset ALL chunk-level vectors to surviving records only
        idx <- which(keep)
        state$CHROM_v <- state$CHROM_v[idx]
        state$POS_v <- state$POS_v[idx]
        state$ID_v <- state$ID_v[idx]
        state$REF_v <- state$REF_v[idx]
        state$ALT_v <- state$ALT_v[idx]
        state$QUAL_v <- state$QUAL_v[idx]
        state$FILTER_v <- state$FILTER_v[idx]
        state$INFO_v <- state$INFO_v[idx]
        state$FORMAT_v <- state$FORMAT_v[idx]
        state$SAMPLE_v <- state$SAMPLE_v[idx]
        state$alts_all <- state$alts_all[idx]
        if (fmt_constant) {
            state$GT_col <- if (!is.null(state$GT_col)) state$GT_col[idx] else NULL
            state$AD_col <- if (!is.null(state$AD_col)) state$AD_col[idx] else NULL
            state$DP_col <- if (!is.null(state$DP_col)) state$DP_col[idx] else NULL
            state$GQ_col <- if (!is.null(state$GQ_col)) state$GQ_col[idx] else NULL
        }
        state$nrow_dt <- length(idx)
    }
    list(state = state, n_dropped_dpgq = n_dropped_dpgq)
}


## ---- Helper 4 -------------------------------------------------------------
#' Rough gene/vc_nonSyn pre-filter on raw INFO string (O7 + O8 PERF tricks).
#'
#' Internal helper for `read.gvr()`'s `convert_chunk()`. Applies two cheap
#' `grepl()`-based pre-filters BEFORE the expensive per-record CSQ expansion:
#' (O7) the gene rough-filter uses a pipe-delimited `|GENE|` pattern; (O8)
#' the `vc_nonSyn` rough-filter scans for the 18 VEP consequence terms that
#' map to the 9 non-synonymous MAF classes. Both filters subset the
#' chunk-level vectors in-place. False positives are tolerated and removed
#' by the exact post-conversion filters in `read.gvr()` step 3d-bis / 3d-ter.
#' Neither filter EVER drops a record that the exact filter would keep.
#'
#' @param state Chunk state list (post-DPGQ filter).
#' @param genes Optional character vector of gene symbols (the `read.gvr()` arg).
#' @param vc_nonSyn Logical or character vector (the `read.gvr()` arg). When
#'   `FALSE`, this rough filter is skipped.
#'
#' @return List `(state = <updated state list>, n_dropped_gene = <int>,
#'   n_dropped_vc = <int>)`.
#'
#' @keywords internal
#' @noRd
.gvr_chunk_filter_rough <- function(state, genes, vc_nonSyn) {
    # --- O7 (PERF): rough gene pre-filter on raw INFO string -----------------------
    # When the `genes` argument is active, most records don't carry any of the
    # requested genes. A cheap grepl on the raw INFO string (which contains the CSQ
    # field with pipe-delimited SYMBOL sub-fields) catches the vast majority of hits
    # and drops the rest BEFORE the expensive CSQ expansion.
    #
    # CRITICAL: gene names in the CSQ field are pipe-delimited (field 4 of each
    # CSQ block), e.g. "T|missense_variant|MODIFIER|TP53|ENSG...". A bare
    # grepl("RET", INFO) would match "INTERVAL", "INTER", etc. -- keeping ~65% of
    # records instead of ~0.05%. The fix is to match |GENE| (gene name between
    # pipes), which is precise enough for a rough filter while still being
    # over-inclusive (a gene name could appear in another pipe-delimited field like
    # DOMAINS). The exact post-conversion gene filter (step 3d-bis) removes false
    # positives. The rough filter never drops a record that the exact filter would
    # keep.
    n_dropped_gene_rough <- 0L
    if (!is.null(genes)) {
        genes_chr <- as.character(genes)
        genes_chr <- genes_chr[!is.na(genes_chr) & nzchar(genes_chr)]
        if (length(genes_chr) > 0L) {
            # Pipe-delimited pattern: |GENE| for each gene.
            # This matches gene names as CSQ field values (between pipes) without
            # matching substrings like "RET" inside "INTERVAL".
            pat <- paste0(paste0("\\|", genes_chr, "\\|"), collapse = "|")
            gene_hit <- grepl(pat, state$INFO_v, ignore.case = TRUE)
            n_dropped_gene_rough <- sum(!gene_hit)
            if (n_dropped_gene_rough > 0L) {
                idx <- which(gene_hit)
                state$CHROM_v <- state$CHROM_v[idx]
                state$POS_v <- state$POS_v[idx]
                state$ID_v <- state$ID_v[idx]
                state$REF_v <- state$REF_v[idx]
                state$ALT_v <- state$ALT_v[idx]
                state$QUAL_v <- state$QUAL_v[idx]
                state$FILTER_v <- state$FILTER_v[idx]
                state$INFO_v <- state$INFO_v[idx]
                state$FORMAT_v <- state$FORMAT_v[idx]
                state$SAMPLE_v <- state$SAMPLE_v[idx]
                state$alts_all <- state$alts_all[idx]
                if (isTRUE(state$fmt_constant)) {
                    state$GT_col <- if (!is.null(state$GT_col)) state$GT_col[idx] else NULL
                    state$AD_col <- if (!is.null(state$AD_col)) state$AD_col[idx] else NULL
                    state$DP_col <- if (!is.null(state$DP_col)) state$DP_col[idx] else NULL
                    state$GQ_col <- if (!is.null(state$GQ_col)) state$GQ_col[idx] else NULL
                }
                state$nrow_dt <- length(idx)
            }
        }
    }

    # --- O8 (PERF): rough vc_nonSyn pre-filter on raw INFO string ------------------
    # When vc_nonSyn is active, ~90% of records are silent (Intron, Silent, IGR, etc.)
    # and will be dropped by the post-conversion filter. But by then the expensive CSQ
    # expansion has already run on ALL of them. The VEP consequence terms that map to
    # the 9 non-synonymous MAF classes are known and finite. If NONE of these terms
    # appear anywhere in the CSQ string, the record CANNOT produce a non-synonymous row.
    # A grepl on the raw INFO string drops ~90% of records before the heavy loop.
    # The exact post-conversion vc_nonSyn filter (step 3d-ter) then tightens any
    # edge cases. This is the same rough-grepl + exact-post-filter architecture as O7.
    n_dropped_vc_rough <- 0L
    if (!identical(vc_nonSyn, FALSE)) {
        # VEP consequence terms that map to the 9 non-synonymous MAF classes.
        # These are the ONLY terms whose most-severe-rank can produce a protein-altering
        # Variant_Classification. If none appear in the Consequence field (field 2 of
        # each CSQ block), the record is guaranteed silent.
        # Note: we use pipe-delimited matching (|TERM| or |TERM&) to avoid false
        # positives like "variant" inside "splice_region_variant". The Consequence
        # field (CSQ field #2) is pipe-delimited like all CSQ fields, and compound
        # consequences use & as separator (e.g. "missense_variant&splice_region_variant").
        nonSyn_vep_terms <- c(
            "splice_acceptor_variant", "splice_donor_variant",
            "transcript_ablation", "exon_loss_variant",
            "stop_gained", "stop_lost", "start_lost",
            "initiator_codon_variant",
            "missense_variant", "conservative_missense_variant",
            "rare_amino_acid_variant", "protein_altering_variant",
            "coding_sequence_variant",
            "frameshift_variant",
            "inframe_insertion", "inframe_deletion",
            "disruptive_inframe_insertion", "disruptive_inframe_deletion"
        )
        # Unlike gene names (O7), VEP consequence terms are long and specific enough
        # that bare substring matching is safe -- they do not appear as substrings
        # of other CSQ field values (verified empirically). Bare matching is also
        # faster than pipe-delimited matching and correctly handles compound
        # consequences (e.g. "missense_variant&splice_region_variant") where the
        # term is preceded by & rather than |.
        pat_vc <- paste0(nonSyn_vep_terms, collapse = "|")
        vc_hit <- grepl(pat_vc, state$INFO_v, ignore.case = FALSE)  # VEP terms are lowercase
        n_dropped_vc_rough <- sum(!vc_hit)
        if (n_dropped_vc_rough > 0L) {
            idx <- which(vc_hit)
            state$CHROM_v <- state$CHROM_v[idx]
            state$POS_v <- state$POS_v[idx]
            state$ID_v <- state$ID_v[idx]
            state$REF_v <- state$REF_v[idx]
            state$ALT_v <- state$ALT_v[idx]
            state$QUAL_v <- state$QUAL_v[idx]
            state$FILTER_v <- state$FILTER_v[idx]
            state$INFO_v <- state$INFO_v[idx]
            state$FORMAT_v <- state$FORMAT_v[idx]
            state$SAMPLE_v <- state$SAMPLE_v[idx]
            state$alts_all <- state$alts_all[idx]
            if (isTRUE(state$fmt_constant)) {
                state$GT_col <- if (!is.null(state$GT_col)) state$GT_col[idx] else NULL
                state$AD_col <- if (!is.null(state$AD_col)) state$AD_col[idx] else NULL
                state$DP_col <- if (!is.null(state$DP_col)) state$DP_col[idx] else NULL
                state$GQ_col <- if (!is.null(state$GQ_col)) state$GQ_col[idx] else NULL
            }
            state$nrow_dt <- length(idx)
        }
    }
    list(state = state, n_dropped_gene = n_dropped_gene_rough, n_dropped_vc = n_dropped_vc_rough)
}


## ---- Helper 5 -------------------------------------------------------------
#' Build per-ALT variant-row lists for a single VCF record.
#'
#' Internal helper for `read.gvr()`'s `convert_chunk()`. Runs the CSQ-block
#' expansion (A3 PERF trick: one vectorised `strsplit()` over all blocks,
#' then `length(f) <- n_csq` to pad), the rank-aware block-to-ALT owner
#' assignment, the canonical-only filter, and the per-ALT row construction
#' for one record. Returns a list of zero or more row lists ready for
#' `data.table::rbindlist()`.
#'
#' This is the inner body that used to live inside
#' `for (r in seq_len(nrow_dt)) { ... }` of `convert_chunk()`. Per-record
#' scalar extraction (chrom/pos/.../GT/AD/DP/GQ) stays inline in the
#' orchestrator because it consumes both the chunk state and the loop
#' index `r`; this helper receives those scalars via `record_ctx`.
#'
#' @param record_ctx Named list with the per-record scalars: `chrom, pos,
#'   vid, ref, altf, qual, filt, info, alts, gt, ad, sdp, gq, ad_vec,
#'   sample_name`.
#' @param csq_fields Character vector of CSQ subfield names.
#' @param n_csq `length(csq_fields)`.
#' @param P Named list of CSQ subfield positions from
#'   `.gvr_chunk_setup_columns()` (`P_Allele`, `P_Cons`, etc.).
#' @param canonical_only Logical scalar (the `read.gvr()` arg).
#' @param ncbi_build Character scalar (the `read.gvr()` arg).
#' @param mstr_cache,v2m_cache Per-call cache environments
#'   (`.gvr_mstr_cache`, `.gvr_v2m_cache`).
#'
#' @return List `(rows = list(<row_1>, <row_2>, ...),
#'   n_dropped_canonical = <int>)`. `rows` may be empty if every ALT was
#'   dropped by `canonical_only`.
#'
#' @keywords internal
#' @noRd
.gvr_chunk_build_record <- function(record_ctx, csq_fields, n_csq, P,
                                    canonical_only, ncbi_build,
                                    mstr_cache, v2m_cache,
                                    normalize_alleles = TRUE) {
    chrom <- record_ctx$chrom
    pos <- record_ctx$pos
    vid <- record_ctx$vid
    ref <- record_ctx$ref
    altf <- record_ctx$altf
    qual <- record_ctx$qual
    filt <- record_ctx$filt
    info <- record_ctx$info
    alts <- record_ctx$alts
    gt <- record_ctx$gt
    ad <- record_ctx$ad
    sdp <- record_ctx$sdp
    gq <- record_ctx$gq
    ad_vec <- record_ctx$ad_vec
    sample_name <- record_ctx$sample_name

    P_Allele <- P$P_Allele
    P_Cons <- P$P_Cons
    P_SYMBOL <- P$P_SYMBOL
    P_Feature <- P$P_Feature
    P_HGVSc <- P$P_HGVSc
    P_HGVSp <- P$P_HGVSp
    P_Existing <- P$P_Existing
    P_CANONICAL <- P$P_CANONICAL
    P_MANESEL <- P$P_MANESEL

    ip <- .gvr_info_parse(info)                       # C1: parse INFO once, then index
    info_DP <- unname(ip["DP"])
    info_AC <- unname(ip["AC"])
    info_AF <- unname(ip["AF"])
    info_MQ <- unname(ip["MQ"])
    info_QD <- unname(ip["QD"])
    info_CNN <- unname(ip["CNN_1D"])
    ac_vec <- if (!is.na(info_AC)) strsplit(info_AC, ",", fixed = TRUE)[[1]] else NA
    af_vec <- if (!is.na(info_AF)) strsplit(info_AF, ",", fixed = TRUE)[[1]] else NA

    csq_raw <- unname(ip["CSQ"])
    csq_blocks <- if (!is.na(csq_raw)) strsplit(csq_raw, ",", fixed = TRUE)[[1]] else character(0)
    # A3 (PERF): split ALL CSQ blocks in ONE vectorized strsplit() call rather than
    # one strsplit() per block in an lapply. base strsplit() drops trailing empty
    # fields, so each result is padded back to n_csq with `length(f) <- n_csq`
    # (NA fill). Byte-identical to the per-block path (verified on 186k real blocks
    # + edge cases); ~22% faster on the isolated split. fixed=TRUE, no new deps.
    block_fields <- lapply(strsplit(csq_blocks, "|", fixed = TRUE), function(f) {
        length(f) <- n_csq
        f
    })
    block_allele <- vapply(block_fields, function(f) {
        a <- f[P_Allele]
        if (is.na(a)) "" else a
    }, character(1))

    # rank-aware CSQ block -> ALT owner assignment
    cand_by_alt <- lapply(alts, function(a) .gvr_vep_allele_candidates(ref, a))
    block_owner <- integer(length(block_allele))
    if (length(block_allele) > 0L) {
        for (bi in seq_along(block_allele)) {
            bs <- block_allele[bi]
            best_alt <- 0L
            best_rank <- .Machine$integer.max
            for (ai2 in seq_along(alts)) {
                rk <- match(bs, cand_by_alt[[ai2]])
                if (!is.na(rk) && rk < best_rank) {
                    best_rank <- rk
                    best_alt <- ai2
                }
            }
            block_owner[bi] <- best_alt
        }
    }

    rows <- list()
    n_dropped_canonical <- 0L
    for (ai in seq_along(alts)) {
        alt <- alts[ai]
        sel <- which(block_owner == ai)
        coords <- .gvr_coords(pos, ref, alt, normalize_alleles = normalize_alleles)
        vt <- coords$var_type

        if (length(sel) > 0L) {
            # vN+11 (P1): coalesce the four per-block vapply()s into a single pre-allocated
            # walk over `sel`, indexing each block_fields[[k]] exactly once. Replaces 4x
            # closure-per-call vapply overhead with one direct loop. Output is byte-identical
            # to the original (same scalar logic, same inputs, same order(...)).
            ns <- length(sel)
            ranks <- integer(ns)
            canon <- integer(ns)
            mane <- integer(ns)
            feat <- character(ns)
            for (i in seq_len(ns)) {
                bf <- block_fields[[sel[i]]]
                ranks[i] <- .gvr_consequence_rank(bf[P_Cons], cache = mstr_cache)
                v <- bf[P_CANONICAL]
                canon[i] <- if (!is.na(v) && v == "YES") 0L else 1L
                v <- bf[P_MANESEL]
                mane[i]  <- if (!is.na(v) && v != "") 0L else 1L
                v <- bf[P_Feature]
                feat[i]  <- if (is.na(v)) "zzz" else v
            }
            ord <- order(ranks, canon, mane, feat)
            chosen <- block_fields[[sel[ord[1]]]]
        } else {
            chosen <- rep(NA_character_, n_csq)
        }

        # vN+4: canonical_only filter -- drop rows whose chosen block is not canonical.
        # Skips the per-ALT row construction entirely (saves Hugo lookup,
        # HGVSp shortening, AD parsing, list() build). Also skips ALTs with no
        # CSQ block at all (chosen all-NA), which is desired: no CANONICAL info => drop.
        if (canonical_only) {
            canon_chosen <- chosen[P_CANONICAL]
            if (is.na(canon_chosen) || canon_chosen != "YES") {
                n_dropped_canonical <- n_dropped_canonical + 1L
                next
            }
        }

        cons_str <- chosen[P_Cons]
        top_term <- .gvr_most_severe_term(cons_str, cache = mstr_cache)
        var_class <- if (is.na(top_term)) "Targeted_Region" else .gvr_vep_to_class(top_term, vt, cache = v2m_cache)
        symbol <- chosen[P_SYMBOL]
        hugo <- if (is.na(symbol) || symbol == "") "Unknown" else symbol

        gc <- .gvr_gt_codes_for_alt(gt, ai)
        # vN+11 (P2): inline what was `map_code` to avoid creating a closure per
        # (record x alt). Identical scalar branches; output byte-identical.
        .c1 <- gc$c1
        t_allele1 <- if (is.na(.c1)) "."
        else if (.c1 == 0L) coords$ref_allele
        else if (.c1 == ai) coords$tum_allele2
        else if (.c1 >= 1L && .c1 <= length(alts)) .gvr_coords(pos, ref, alts[.c1], normalize_alleles = normalize_alleles)$tum_allele2
        else "."
        .c2 <- gc$c2
        t_allele2 <- if (is.na(.c2)) "."
        else if (.c2 == 0L) coords$ref_allele
        else if (.c2 == ai) coords$tum_allele2
        else if (.c2 >= 1L && .c2 <= length(alts)) .gvr_coords(pos, ref, alts[.c2], normalize_alleles = normalize_alleles)$tum_allele2
        else "."

        t_ref_count <- if (length(ad_vec) >= 1 && !any(is.na(ad_vec))) ad_vec[1] else NA_character_
        t_alt_count <- if (length(ad_vec) >= (ai + 1) && !any(is.na(ad_vec))) ad_vec[ai + 1] else NA_character_

        exist <- chosen[P_Existing]
        dbsnp <- NA_character_
        if (!is.na(exist) && exist != "") {
            rs <- grep("^rs", strsplit(exist, "&", fixed = TRUE)[[1]], value = TRUE)
            if (length(rs)) dbsnp <- rs[1]
        }
        if (is.na(dbsnp) && !is.na(vid) && vid != ".") {
            rs <- grep("^rs", strsplit(vid, ";", fixed = TRUE)[[1]], value = TRUE)
            if (length(rs)) dbsnp <- rs[1]
        }

        hgvsp <- .gvr_url_decode(chosen[P_HGVSp])
        hgvsc <- .gvr_url_decode(chosen[P_HGVSc])
        vep_vals <- as.list(chosen)
        names(vep_vals) <- csq_fields

        rows[[length(rows) + 1L]] <- c(
            list(
                Hugo_Symbol = hugo, Entrez_Gene_Id = "0", Center = ".", NCBI_Build = ncbi_build,
                Chromosome = chrom, Start_Position = coords$start, End_Position = coords$end,
                Strand = "+", Variant_Classification = var_class, Variant_Type = vt,
                Reference_Allele = coords$ref_allele, Tumor_Seq_Allele1 = t_allele1,
                Tumor_Seq_Allele2 = t_allele2, dbSNP_RS = if (is.na(dbsnp)) "" else dbsnp,
                Tumor_Sample_Barcode = sample_name, Match_Norm_Seq_Allele1 = "",
                Match_Norm_Seq_Allele2 = "", HGVSc = if (is.na(hgvsc)) "" else hgvsc,
                HGVSp = if (is.na(hgvsp)) "" else hgvsp, HGVSp_Short = .gvr_make_hgvsp_short(hgvsp),
                Transcript_ID = {
                    v <- chosen[P_Feature]
                    if (is.na(v)) "" else v
                },
                Consequence = if (is.na(cons_str)) "" else cons_str,
                t_depth = if (is.na(sdp)) "" else sdp,
                t_ref_count = if (is.na(t_ref_count)) "" else t_ref_count,
                t_alt_count = if (is.na(t_alt_count)) "" else t_alt_count
            ),
            vep_vals,
            list(
                FILTER = filt,
                QUAL = qual,
                INFO_DP = if (is.na(info_DP)) "" else info_DP,
                INFO_AC = if (length(ac_vec) >= ai && !any(is.na(ac_vec))) ac_vec[ai] else if (!is.na(info_AC)) info_AC else "",
                INFO_AF = if (length(af_vec) >= ai && !any(is.na(af_vec))) af_vec[ai] else if (!is.na(info_AF)) info_AF else "",
                INFO_MQ = if (is.na(info_MQ)) "" else info_MQ,
                INFO_QD = if (is.na(info_QD)) "" else info_QD,
                CNN_1D = if (is.na(info_CNN)) "" else info_CNN,
                GT = gt, AD = if (is.na(ad)) "" else ad,
                sample_DP = if (is.na(sdp)) "" else sdp, GQ = if (is.na(gq)) "" else gq
            )
        )
    }
    list(rows = rows, n_dropped_canonical = n_dropped_canonical)
}
