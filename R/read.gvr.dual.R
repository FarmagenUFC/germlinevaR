# =============================================================================
# read.gvr.dual.R
# -----------------------------------------------------------------------------
# Reader for VCFs annotated with BOTH Ensembl VEP (##INFO=<ID=CSQ>) and SnpEff
# (##INFO=<ID=ANN>) on the same records. Sibling to read.gvr() (VEP-only) and
# read.gvr.snpeff() (SnpEff-only).
#
# Public function: read.gvr.dual(folder, vcf_path, ...) -- same signature and
# defaults as read.gvr(). Returns a data.table with the canonical table schema
# (one row per ALT allele, VEP-driven transcript pick), augmented with:
#   - all VEP CSQ columns including FREQS (81st field, when present in header)
#   - 4 SnpEff parallel comparison columns: snpeff_consequence, snpeff_impact,
#     snpeff_gene, snpeff_hgvsc
#   - 4 SnpEff LoF/NMD prediction columns: LOF_Gene, LOF_Pct_Transcripts,
#     NMD_Gene, NMD_Pct_Transcripts
# Output is tagged with attr(., 'annotator') = 'dual'.
#
# Implementation approach: delegates the VEP CSQ parse + table construction to
# read.gvr() via the private force-annotator flag in .gvr_internal_env
# routing loop). Then runs a second, lightweight pass over the same VCF files
# to extract ANN/LOF/NMD INFO fields and joins them onto the MAF on
# (Tumor_Sample_Barcode, Chromosome, Start_Position, Reference_Allele,
# Tumor_Seq_Allele2). The VEP-driven row count and ordering are preserved
# exactly; the SnpEff annotations are added as new columns only.
#
# Field-level policy: VEP wins. The canonical VEP columns (Hugo_Symbol,
# Consequence, IMPACT, HGVSc, HGVSp, etc.) hold VEP's pick exclusively. The
# 4 snpeff_* columns hold SnpEff's parallel pick from the matching ANN block
# (same allele, same gene as VEP; falls back to allele-only match if no
# gene-matching block exists). The 4 LOF_*/NMD_* columns are extracted from
# SnpEff's INFO LOF= and NMD= keys.
# =============================================================================


#' Convert dual-annotated germline VCF(s) (VEP + SnpEff) to an MAF-like data.table
#'
#' @description
#' Converts single-sample germline VCFs that carry **both** Ensembl VEP `CSQ`
#' and SnpEff `ANN` annotation INFO fields on the same records (typical
#' workflow: SnpEff -> VEP, or vice versa). Returns an MAF-like
#' `data.table` with one row per ALT allele, using VEP's transcript pick as
#' the spine and adding SnpEff-derived comparison columns and LoF/NMD
#' predictions.
#'
#' This function is the dual-annotator sibling of [read.gvr()] (VEP-only) and
#' [read.gvr.snpeff()] (SnpEff-only). All three are normally invoked through
#' [read.gvr()], which auto-routes inputs based on the INFO tags found in the
#' VCF header.
#'
#' @details
#' Field-level priority: **VEP wins.** The canonical table columns and all 80
#' (or 81 with `FREQS`) VEP-style CSQ columns hold VEP's pick exclusively.
#' SnpEff data is added in additional columns:
#' \itemize{
#'   \item `snpeff_consequence`, `snpeff_impact`, `snpeff_gene`, `snpeff_hgvsc`
#'     -- the SnpEff ANN block matching the same ALT allele and gene as VEP's
#'     pick (most-deleterious block when multiple match). Empty when SnpEff
#'     has no ANN block for that allele.
#'   \item `LOF_Gene`, `LOF_Pct_Transcripts` -- SnpEff's loss-of-function
#'     prediction (gene flagged, fraction of transcripts affected) from the
#'     INFO `LOF=` field. Empty when no LoF call.
#'   \item `NMD_Gene`, `NMD_Pct_Transcripts` -- SnpEff's nonsense-mediated-
#'     decay prediction from the INFO `NMD=` field. Empty when no NMD call.
#' }
#'
#' Transcript pick: VEP drives the row spine (one row per ALT allele). If VEP
#' was run with `--per_gene` (one transcript per gene per allele), the chosen
#' transcript is VEP's most-severe across genes for that allele; otherwise the
#' standard VEP most-severe-block ranking applies. SnpEff annotations are
#' attached *to* that row from the ANN block matching the same allele and
#' (preferentially) the same `Gene_Name`.
#'
#' Schema additions vs `read.gvr()` (VEP-only):
#' \itemize{
#'   \item `FREQS` -- 81st canonical CSQ field, populated when VEP was run
#'     with `--everything` (v113+). Names the gnomAD population where MAX_AF
#'     was observed (e.g. "gnomADg_ASJ"). Empty when the CSQ header does not
#'     include this field.
#'   \item 4 LOF/NMD columns (above)
#'   \item 4 snpeff_* parallel columns (above)
#' }
#' Total: ~124 columns (vs ~116 for `read.gvr()` and `read.gvr.snpeff()`).
#'
#' All other behaviour, defaults, arguments, and filtering options are
#' identical to [read.gvr()]; see that function's documentation for the full
#' option reference.
#'
#' @inheritParams read.gvr
#' @return A `data.table` with one row per ALT allele, ~124 columns. Carries
#'   `attr(., "annotator") = "dual"`.
#' @seealso [read.gvr()], [read.gvr.snpeff()]
#' @export
#' @importFrom data.table data.table as.data.table fread setDT setnames setkey rbindlist :=
#' @examples
#' ## The function signature is exported and callable:
#' is.function(read.gvr.dual)
#'
#' \dontrun{
#'   ## read.gvr.dual() expects VCFs annotated with BOTH VEP (CSQ INFO field)
#'   ## AND SnpEff (ANN INFO field) in the same record. The shipped example
#'   ## VCF is VEP-only, so a real dual-annotated example needs your own
#'   ## VCFs:
#'   gvr <- read.gvr.dual(folder = "/path/to/dual-annotated-vcfs/")
#'
#'   ## Or via the auto-router in read.gvr() when the VCF header declares
#'   ## both VEP and SnpEff INFO fields:
#'   gvr <- read.gvr("/path/to/dual-annotated-vcfs/")
#'
#'   ## Compare VEP vs SnpEff picks on high-impact variants:
#'   gvr[IMPACT == "HIGH" & snpeff_impact != "" & IMPACT != snpeff_impact,
#'       .(Hugo_Symbol, Consequence, IMPACT, snpeff_gene, snpeff_consequence,
#'         snpeff_impact)]
#' }
read.gvr.dual <- function(folder = ".",
                          vcf_path   = NULL,
                          file       = NULL,
                          pattern    = "\\.vcf\\.gz$",
                          write_tsv  = FALSE,
                          write_rds  = FALSE,
                          write_xlsx = FALSE,
                          out_dir    = NULL,
                          out_prefix = NULL,
                          chunk_size = 25000L,
                          ncbi_build = "auto",
                          add_genotype      = TRUE,
                          strip_hgvs_prefix = TRUE,
                          dedup_columns     = TRUE,
                          drop_empty_cols   = FALSE,
                          add_abraom        = TRUE,
                          abraom_path       = NULL,
                          abraom_url        = "https://abraom.ib.usp.br/download/ABRaOM_60+_SABE_609_exomes_annotated.gz",
                          cache_dir         = NULL,
                          min_DP            = 10,
                          min_GQ            = 30,
                          genes             = NULL,
                          panel             = NULL,
                          vc_nonSyn         = FALSE,
                          canonical_only    = TRUE,
                          ncores            = 1L,
                          verbose    = TRUE) {

  # ---------------------------------------------------------------------------
  # Phase 1: VEP-driven table construction.
  # ---------------------------------------------------------------------------
  # Delegate to read.gvr() with the force-annotator flag set to "vep" so the
  # VEP body runs even though .detect_annotator() would return "dual" for
  # these files. This gives us the canonical table with one row per ALT allele,
  # the canonical 80 (or 81 with FREQS) CSQ columns auto-detected from the
  # actual CSQ header, ABraOM join, build detection, DP/GQ filter, etc.
  # ---------------------------------------------------------------------------
  # Set the force-annotator flag in the package's private env, with on.exit
  # cleanup as a belt-and-braces guarantee that the flag is cleared even if
  # read.gvr() errors (it would normally be consumed-and-cleared inside read.gvr()).
  .gvr_set_force_annotator("vep")
  on.exit(.gvr_set_force_annotator(NULL), add = TRUE)

  gvr <- read.gvr(
    folder            = folder,
    vcf_path          = vcf_path,
    file              = file,
    pattern           = pattern,
    write_tsv         = FALSE,     # we'll write at the end after SnpEff join
    write_rds         = FALSE,
    write_xlsx        = FALSE,
    out_dir           = out_dir,
    out_prefix        = out_prefix,
    chunk_size        = chunk_size,
    ncbi_build        = ncbi_build,
    add_genotype      = add_genotype,
    strip_hgvs_prefix = strip_hgvs_prefix,
    dedup_columns     = dedup_columns,
    drop_empty_cols   = FALSE,     # keep all cols; we'll drop at end if requested
    add_abraom        = add_abraom,
    abraom_path       = abraom_path,
    abraom_url        = abraom_url,
    cache_dir         = cache_dir,
    min_DP            = min_DP,
    min_GQ            = min_GQ,
    genes             = genes,
    panel             = panel,
    vc_nonSyn         = vc_nonSyn,
    canonical_only    = canonical_only,
    ncores            = ncores,
    verbose           = verbose
  )

  if (!data.table::is.data.table(gvr) || nrow(gvr) == 0L) {
    # Nothing to do: empty result, tag and return.
    data.table::setattr(gvr, "annotator", "dual")
    return(gvr)
  }

  # ---------------------------------------------------------------------------
  # Phase 2: SnpEff ANN/LOF/NMD extraction from the same VCF files.
  # ---------------------------------------------------------------------------
  # We need to know which VCF files contributed to `gvr`. read.gvr() doesn't
  # expose this directly, so re-resolve the same set using the same arguments.
  # ---------------------------------------------------------------------------
  vcf_paths <- .gvr_dual_resolve_files(folder = folder, vcf_path = vcf_path,
                                       file = file, pattern = pattern,
                                       verbose = FALSE)

  if (isTRUE(verbose)) message("read.gvr.dual: Phase 2 -- extracting SnpEff ANN/LOF/NMD from ",
                               length(vcf_paths), " file(s) ...")

  # Build one long lookup table of SnpEff annotations keyed by
  # (sample, chrom, pos, ref, alt). Each VCF file is processed in chunks so
  # we never hold the whole decompressed file in memory.
  snpeff_tab <- .gvr_dual_extract_snpeff(
    vcf_paths   = vcf_paths,
    chunk_size  = chunk_size,
    verbose     = verbose
  )

  if (is.null(snpeff_tab) || nrow(snpeff_tab) == 0L) {
    if (isTRUE(verbose))
      message("  No SnpEff ANN/LOF/NMD records extracted; snpeff_*/LOF_*/NMD_* columns will be empty.")
    snpeff_tab <- data.table::data.table(
      Tumor_Sample_Barcode = character(0),
      Chromosome           = character(0),
      Start_Position       = integer(0),
      Reference_Allele     = character(0),
      Tumor_Seq_Allele2    = character(0),
      snpeff_consequence   = character(0),
      snpeff_impact        = character(0),
      snpeff_gene          = character(0),
      snpeff_hgvsc         = character(0),
      LOF_genes_list       = list(),
      LOF_pcts_list        = list(),
      NMD_genes_list       = list(),
      NMD_pcts_list        = list()
    )
  }

  # ---------------------------------------------------------------------------
  # Phase 3: Join SnpEff annotations onto the VEP MAF spine.
  # ---------------------------------------------------------------------------
  # Key: (Tumor_Sample_Barcode, Chromosome, Start_Position, Reference_Allele,
  #       Tumor_Seq_Allele2).
  #
  # NOTE: read.gvr()'s gvr_coords() converts indels to an MAF-like
  # representation (e.g. deletion "TC->-", insertion "-AT"). We must apply
  # the SAME normalization to the SnpEff side before joining.
  # .gvr_dual_extract_snpeff() already emits MAF-like ref/alt using a
  # verbatim copy of gvr_coords(), so the keys align.
  # ---------------------------------------------------------------------------
  matched <- .gvr_dual_attach_snpeff(gvr, snpeff_tab, verbose = verbose)

  if (isTRUE(verbose)) {
    n_match_gene <- sum(matched$snpeff_gene != "" & matched$Hugo_Symbol != "" &
                         matched$snpeff_gene == matched$Hugo_Symbol, na.rm = TRUE)
    n_match_any  <- sum(matched$snpeff_gene != "", na.rm = TRUE)
    pct_gene <- if (nrow(matched)) round(100 * n_match_gene / nrow(matched), 1) else 0
    pct_any  <- if (nrow(matched)) round(100 * n_match_any  / nrow(matched), 1) else 0
    message(sprintf("  SnpEff cross-annotation: gene match %d/%d rows (%.1f%%); allele match %d/%d (%.1f%%).",
                    n_match_gene, nrow(matched), pct_gene,
                    n_match_any,  nrow(matched), pct_any))
  }

  # ---------------------------------------------------------------------------
  # Phase 4: Drop empty columns + write side-effect files, if requested.
  # ---------------------------------------------------------------------------
  if (isTRUE(drop_empty_cols)) {
    is_empty_col <- function(x) {
      if (is.numeric(x)) all(is.na(x))
      else all(is.na(x) | x == "")
    }
    keep <- !vapply(matched, is_empty_col, logical(1L))
    if (any(!keep)) {
      dropped <- names(matched)[!keep]
      matched <- matched[, .SD, .SDcols = which(keep)]
      if (isTRUE(verbose))
        message("  drop_empty_cols=TRUE: removed ", sum(!keep), " all-empty column(s): ",
                paste(dropped, collapse = ", "))
    }
  }

  # Tag and final-dimension diagnostic
  data.table::setattr(matched, "annotator", "dual")
  data.table::setattr(matched, "ncbi_build", attr(gvr, "ncbi_build"))

  if (isTRUE(verbose))
    message(sprintf("Final Table Dimensions: %d rows x %d columns.",
                    nrow(matched), ncol(matched)))

  # Side-effect files (TSV / RDS / XLSX) -- mirror read.gvr() behaviour.
  if (isTRUE(write_tsv) || isTRUE(write_rds) || isTRUE(write_xlsx)) {
    out_dir_eff <- if (is.null(out_dir)) {
      if (!is.null(vcf_path)) dirname(vcf_path[1L]) else folder
    } else out_dir
    if (!dir.exists(out_dir_eff)) dir.create(out_dir_eff, recursive = TRUE, showWarnings = FALSE)
    prefix_eff <- if (is.null(out_prefix)) "gvr_dual" else out_prefix
    base <- file.path(out_dir_eff, prefix_eff)
    if (isTRUE(write_tsv)) {
      tsv_path <- paste0(base, ".tsv")
      data.table::fwrite(matched, tsv_path, sep = "\t", na = "")
      if (isTRUE(verbose)) message("  Wrote TSV: ", tsv_path)
    }
    if (isTRUE(write_rds)) {
      rds_path <- paste0(base, ".rds")
      saveRDS(matched, rds_path)
      if (isTRUE(verbose)) message("  Wrote RDS: ", rds_path)
    }
    if (isTRUE(write_xlsx)) {
      xlsx_path <- paste0(base, ".xlsx")
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "gvr_table")
      openxlsx::writeData(wb, "gvr_table", as.data.frame(matched))
      openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
      if (isTRUE(verbose)) message("  Wrote XLSX: ", xlsx_path)
    }
  }

  matched
}


# =============================================================================
# Nested helpers (not exported)
# -----------------------------------------------------------------------------
# .gvr_dual_resolve_files()       -- locate VCFs from folder/pattern/file/vcf_path
# .gvr_dual_extract_snpeff()      -- second-pass parse: ANN + LOF + NMD per record
# .gvr_dual_attach_snpeff()       -- join SnpEff lookup onto VEP MAF spine
# .gvr_dual_gvr_coords()          -- pos/ref/alt -> MAF-like (start, ref, alt)
# .gvr_dual_strip_snpeff_allele() -- normalize SnpEff non-standard allele forms
# .gvr_dual_pull_info()           -- vectorised extractor for ;-separated INFO KEY=VAL
# .gvr_dual_header_sample()       -- sample name from #CHROM header
# .gvr_dual_header_ann_fields()   -- ANN field names from ##INFO=<ID=ANN,...>
# .gvr_dual_parse_lof_nmd_record() -- parse a record's LOF= / NMD= INFO value
# =============================================================================


# Locate VCF files using the same conventions as read.gvr() / read.gvr.snpeff().
# Used for Phase 2 (we know which files contributed to the Phase-1 MAF).
.gvr_dual_resolve_files <- function(folder, vcf_path, file, pattern, verbose) {
  if (!is.null(vcf_path)) {
    vp <- vcf_path
    bad <- !file.exists(vp)
    if (any(bad)) stop("read.gvr.dual: file(s) not found: ",
                        paste(vp[bad], collapse = ", "), call. = FALSE)
    return(normalizePath(vp))
  }
  if (!is.null(file)) {
    vp <- file.path(folder, file)
    bad <- !file.exists(vp)
    if (any(bad)) stop("read.gvr.dual: file(s) not found in folder: ",
                        paste(vp[bad], collapse = ", "), call. = FALSE)
    return(normalizePath(vp))
  }
  if (!dir.exists(folder))
    stop("read.gvr.dual: folder not found: ", folder, call. = FALSE)
  vp <- list.files(folder, pattern = pattern, full.names = TRUE,
                   recursive = FALSE, ignore.case = TRUE)
  if (!length(vp))
    stop("read.gvr.dual: no files matching pattern '", pattern,
         "' in folder: ", folder, call. = FALSE)
  normalizePath(vp)
}


# MAF-like coords for one (pos, ref, alt). VERBATIM COPY of read.gvr.R's
# nested gvr_coords() helper (defined inside read.gvr() at lines ~910-930).
# Keeping this identical is REQUIRED so the SnpEff side of the dual reader
# produces the same join keys (Reference_Allele, Tumor_Seq_Allele2,
# Start_Position) that VEP path emits. If read.gvr.R's gvr_coords ever
# changes, mirror the change here.
#
# Returns a list: (var_type, start, end, ref_allele, tum_allele2)
.gvr_dual_gvr_coords <- function(pos, ref, alt) {
  pos <- as.integer(pos)
  if (is.na(alt) || alt == "*")
    return(list(var_type = "DEL", start = pos, end = pos,
                ref_allele = ref, tum_allele2 = "*"))
  rl <- nchar(ref); al <- nchar(alt)
  if (rl == al && rl == 1L)
    return(list(var_type = "SNP", start = pos, end = pos,
                ref_allele = ref, tum_allele2 = alt))
  if (rl == al && rl > 1L) {
    vt <- switch(as.character(rl), "2" = "DNP", "3" = "TNP", "ONP")
    return(list(var_type = vt, start = pos, end = pos + rl - 1L,
                ref_allele = ref, tum_allele2 = alt))
  }
  if (al > rl) {                                   # insertion (left-anchored)
    ins <- substr(alt, rl + 1L, al)
    return(list(var_type = "INS", start = pos + rl - 1L, end = pos + rl,
                ref_allele = "-", tum_allele2 = ins))
  } else {                                         # deletion (left-anchored)
    del <- substr(ref, al + 1L, rl)
    return(list(var_type = "DEL", start = pos + al, end = pos + rl - 1L,
                ref_allele = del, tum_allele2 = "-"))
  }
}


# Strip SnpEff non-standard allele forms (Cingolani spec):
#   "G-C"                -> "G"     (cancer-somatic-vs-germline)
#   "C-chr1:123456_A>T"  -> "C"     (compound)
# Standard alleles (no dash) pass through unchanged.
.gvr_dual_strip_snpeff_allele <- function(a) {
  if (is.na(a) || !nzchar(a)) return(a)
  idx <- regexpr("-", a, fixed = TRUE)
  if (idx > 0L) substr(a, 1L, idx - 1L) else a
}


# Pull the value of one INFO key (e.g. "ANN") from a vector of INFO strings.
# Returns a character vector the same length as `info`; NA where the key is
# absent. KEY=VAL pairs are ;-separated; INFO strings are not URL-decoded.
.gvr_dual_pull_info <- function(info, key) {
  # Anchor at start-or-semicolon to avoid matching "<other>ANN=" suffixes.
  pat <- paste0("(?:^|;)", key, "=([^;]*)")
  m <- regmatches(info, regexec(pat, info, perl = TRUE))
  vapply(m, function(x) if (length(x) >= 2L) x[2L] else NA_character_, character(1L))
}


# Header parsing: sample name from the #CHROM line. Matches the convention
# used in read.gvr.R (last column of the #CHROM header line).
.gvr_dual_header_sample <- function(path) {
  con <- gzfile(path, "r"); on.exit(close(con))
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) stop("read.gvr.dual: #CHROM not found in ", basename(path))
    if (startsWith(line, "#CHROM")) {
      cols <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      return(cols[length(cols)])
    }
  }
}


# Header parsing: ANN field-name vector from ##INFO=<ID=ANN ... Description=...>.
# Uses the same approach as read.gvr.snpeff.R's get_ann_fields(): pull the
# substring between the FIRST and LAST single-quote on the line, then split
# on "|" and trim. Robust to extra single-quotes elsewhere on the line.
.gvr_dual_header_ann_fields <- function(path) {
  con <- gzfile(path, "r"); on.exit(close(con))
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) return(NULL)
    if (startsWith(line, "#CHROM")) return(NULL)
    if (startsWith(line, "##INFO=<ID=ANN,")) {
      qpos <- as.integer(gregexpr("'", line, fixed = TRUE)[[1L]])
      if (length(qpos) < 2L || qpos[1L] < 0L) return(NULL)
      inner <- substr(line, qpos[1L] + 1L, qpos[length(qpos)] - 1L)
      return(trimws(strsplit(inner, "|", fixed = TRUE)[[1L]]))
    }
  }
}


# Parse one record's LOF= or NMD= INFO value into a per-gene lookup table.
# Returns list(gene=character, pct=numeric); empty vectors if the value is NA.
# SnpEff format (Cingolani): "(Gene|GeneID|N|Pct)" for single, or
# "(g1|...|0.5),(g2|...|1.0)" for multi-gene calls.
.gvr_dual_parse_lof_nmd_record <- function(s) {
  if (is.na(s) || !nzchar(s)) return(list(gene = character(0L), pct = numeric(0L)))
  s <- sub("^\\(", "", s); s <- sub("\\)$", "", s)
  tuples <- strsplit(s, "),(", fixed = TRUE)[[1L]]
  genes <- character(length(tuples))
  pcts  <- numeric(length(tuples))
  for (i in seq_along(tuples)) {
    parts <- strsplit(tuples[i], "|", fixed = TRUE)[[1L]]
    genes[i] <- if (length(parts) >= 1L) parts[1L] else ""
    pcts[i]  <- if (length(parts) >= 4L) suppressWarnings(as.numeric(parts[4L])) else NA_real_
  }
  list(gene = genes, pct = pcts)
}


# Extract SnpEff ANN + LOF + NMD from VCF files, returning a long-format
# data.table keyed by (sample, chrom, pos, ref, alt) using MAF-like alleles
# (deletions as "TC"->"-", insertions as "-"->"AT"; matches read.gvr()'s
# gvr_coords() output exactly).
#
# Per-row pre-selection of the best ANN block requires the VEP Hugo_Symbol
# for gene-aware matching. We do that at attach time, not here -- here we
# store all gene-matched candidate blocks as a list-column ("snpeff_blocks")
# along with the field-position lookup ("snpeff_field_idx"), and the full
# LOF/NMD gene-lookup maps as list-columns ("LOF_genes_list",
# "LOF_pcts_list", "NMD_genes_list", "NMD_pcts_list"). Attach then picks
# the gene-matching candidate per row.
.gvr_dual_extract_snpeff <- function(vcf_paths, chunk_size, verbose) {
  all_rows <- vector("list", length(vcf_paths))

  for (fi in seq_along(vcf_paths)) {
    path <- vcf_paths[fi]
    if (isTRUE(verbose))
      message(sprintf("  [file %d/%d] SnpEff scan %s",
                      fi, length(vcf_paths), basename(path)))

    # Read header to (a) find sample name and (b) get ANN field positions.
    sample_name <- .gvr_dual_header_sample(path)
    ann_fields  <- .gvr_dual_header_ann_fields(path)
    if (is.null(ann_fields))
      stop("read.gvr.dual: ##INFO=<ID=ANN> header not found in ",
           basename(path), call. = FALSE)

    A_Allele       <- match("Allele",            ann_fields)
    A_Annotation   <- match("Annotation",        ann_fields)
    A_Impact       <- match("Annotation_Impact", ann_fields)
    A_Gene_Name    <- match("Gene_Name",         ann_fields)
    A_HGVSc        <- match("HGVS.c",            ann_fields)
    field_idx_const <- c(A_Allele, A_Annotation, A_Impact, A_Gene_Name, A_HGVSc)

    # Stream the body in chunks. We only need columns CHROM POS REF ALT INFO.
    con <- gzfile(path, "r"); on.exit(close(con), add = TRUE)
    # Skip header
    repeat {
      line <- readLines(con, n = 1L, warn = FALSE)
      if (length(line) == 0L) break
      if (startsWith(line, "#CHROM")) break
    }

    file_rows <- list()
    rep_lines <- 0L
    repeat {
      lines <- readLines(con, n = chunk_size, warn = FALSE)
      if (!length(lines)) break
      # Pull columns 1, 2, 4, 5, 8 (CHROM POS REF ALT INFO).
      # tstrsplit(keep=...) returns one list element per *kept* column, in the
      # order given. parts[[3]] is the third kept column == source col 4 (REF).
      parts <- data.table::tstrsplit(lines, "\t", fixed = TRUE,
                                     keep = c(1L, 2L, 4L, 5L, 8L),
                                     type.convert = FALSE)
      chrom <- parts[[1L]]; pos <- as.integer(parts[[2L]])
      ref   <- parts[[3L]]; alt <- parts[[4L]]
      info  <- parts[[5L]]

      # Extract ANN / LOF / NMD substrings from INFO. INFO is ;-separated KEY=VAL.
      ann_str <- .gvr_dual_pull_info(info, "ANN")
      lof_str <- .gvr_dual_pull_info(info, "LOF")
      nmd_str <- .gvr_dual_pull_info(info, "NMD")

      # Expand multi-ALT into one row per ALT (matches MAF-row granularity).
      n_rec <- length(chrom)
      for (ri in seq_len(n_rec)) {
        if (is.na(pos[ri])) next
        alts_ri <- strsplit(alt[ri], ",", fixed = TRUE)[[1L]]
        # Parse ANN blocks once per record
        if (!is.na(ann_str[ri]) && nzchar(ann_str[ri])) {
          blocks_str <- strsplit(ann_str[ri], ",", fixed = TRUE)[[1L]]
          block_fields <- strsplit(blocks_str, "|", fixed = TRUE)
          block_allele <- if (!is.na(A_Allele))
            vapply(block_fields, function(b)
              if (length(b) >= A_Allele) .gvr_dual_strip_snpeff_allele(b[A_Allele]) else "",
              character(1L))
          else character(0L)
        } else {
          block_fields <- list()
          block_allele <- character(0L)
        }

        # Parse LOF / NMD into (gene -> pct) lookups
        lof_map <- .gvr_dual_parse_lof_nmd_record(lof_str[ri])
        nmd_map <- .gvr_dual_parse_lof_nmd_record(nmd_str[ri])

        for (ai in seq_along(alts_ri)) {
          this_alt <- alts_ri[ai]
          # MAF-like coords for this (ref, alt) pair: VERBATIM read.gvr() logic
          mc <- .gvr_dual_gvr_coords(pos[ri], ref[ri], this_alt)
          # Find ANN blocks matching this ALT (after SnpEff allele stripping)
          if (length(block_allele)) {
            sel <- which(block_allele == this_alt)
            matched_blocks <- if (length(sel)) block_fields[sel] else list()
          } else {
            matched_blocks <- list()
          }

          file_rows[[length(file_rows) + 1L]] <- list(
            Tumor_Sample_Barcode = sample_name,
            Chromosome           = chrom[ri],
            Start_Position       = mc$start,
            Reference_Allele     = mc$ref_allele,
            Tumor_Seq_Allele2    = mc$tum_allele2,
            snpeff_blocks        = list(matched_blocks),
            snpeff_field_idx     = list(field_idx_const),
            LOF_genes_list       = list(lof_map$gene),
            LOF_pcts_list        = list(lof_map$pct),
            NMD_genes_list       = list(nmd_map$gene),
            NMD_pcts_list        = list(nmd_map$pct)
          )
        }
      }
      rep_lines <- rep_lines + length(lines)
      if (isTRUE(verbose) && (rep_lines %% (chunk_size * 4L) == 0L))
        message(sprintf("     %s records scanned", format(rep_lines, big.mark = ",")))
    }
    close(con); on.exit()

    if (length(file_rows))
      all_rows[[fi]] <- data.table::rbindlist(file_rows, fill = TRUE)
  }

  out <- data.table::rbindlist(all_rows, fill = TRUE)
  if (nrow(out) == 0L) return(out)
  data.table::setkey(out, Tumor_Sample_Barcode, Chromosome, Start_Position,
                    Reference_Allele, Tumor_Seq_Allele2)
  out
}


# Attach SnpEff lookup table to VEP MAF spine. snpeff_tab has one row per
# (sample, chrom, pos, ref, alt) -- the MAF granularity. For each table row we
# look up the matching snpeff_tab row, then pick the best ANN block:
#   1) ANN block whose Gene_Name matches Hugo_Symbol (preferred)
#   2) Else the first ANN block matching the ALT allele
#   3) Else empty
# LOF/NMD: prefer the gene matching Hugo_Symbol; else first gene listed.
.gvr_dual_attach_snpeff <- function(gvr, snpeff_tab, verbose) {
  m <- data.table::copy(gvr)

  # Default-empty columns first (so every table row has them even if no snpeff_tab match).
  m[, snpeff_consequence := ""]
  m[, snpeff_impact      := ""]
  m[, snpeff_gene        := ""]
  m[, snpeff_hgvsc       := ""]
  m[, LOF_Gene           := ""]
  m[, LOF_Pct_Transcripts := NA_real_]
  m[, NMD_Gene           := ""]
  m[, NMD_Pct_Transcripts := NA_real_]

  if (!nrow(snpeff_tab)) return(m)

  # Build a left-join: m -> snpeff_tab on the 5 key columns. We want one
  # snpeff_tab row per m row in m's order. data.table syntax x[i, on=]
  # right-joins on x using i's keys, but emits rows in i's order. We use
  # snpeff_tab[m, on=...] which returns nrow(m) rows in m's order, with
  # snpeff_tab's non-key cols attached as NA where unmatched.
  join_keys <- c("Tumor_Sample_Barcode", "Chromosome", "Start_Position",
                 "Reference_Allele", "Tumor_Seq_Allele2")

  # Only carry the snpeff_tab columns we need (not the m cols, which we
  # already have on `m`). The join result `joined` has all of m's cols PLUS
  # snpeff_tab's non-key cols, in m's row order.
  joined <- snpeff_tab[m, on = join_keys, allow.cartesian = FALSE]

  # Per-row resolution: walk the snpeff_blocks list-column, pick best block.
  n <- nrow(joined)
  s_cons   <- character(n)
  s_imp    <- character(n)
  s_gene   <- character(n)
  s_hgvsc  <- character(n)
  l_gene   <- character(n)
  l_pct    <- rep(NA_real_, n)
  d_gene   <- character(n)
  d_pct    <- rep(NA_real_, n)

  blocks_col   <- joined$snpeff_blocks
  idx_col      <- joined$snpeff_field_idx
  hugo_col     <- joined$Hugo_Symbol
  lof_genes_col <- joined$LOF_genes_list
  lof_pcts_col  <- joined$LOF_pcts_list
  nmd_genes_col <- joined$NMD_genes_list
  nmd_pcts_col  <- joined$NMD_pcts_list

  for (i in seq_len(n)) {
    blocks_i <- blocks_col[[i]]
    if (length(blocks_i)) {
      idx_i    <- idx_col[[i]]
      # idx_i: c(A_Allele, A_Annotation, A_Impact, A_Gene_Name, A_HGVSc)
      A_Anno  <- idx_i[2L]
      A_Imp   <- idx_i[3L]
      A_Gene  <- idx_i[4L]
      A_Hgvsc <- idx_i[5L]

      # Try gene match against Hugo_Symbol; else default to first block (= most-severe)
      pick <- 1L
      if (!is.na(A_Gene) && !is.na(hugo_col[i]) && nzchar(hugo_col[i])) {
        gene_v <- vapply(blocks_i, function(b)
          if (length(b) >= A_Gene) b[A_Gene] else "", character(1L))
        gm <- which(gene_v == hugo_col[i])
        if (length(gm)) pick <- gm[1L]  # most-severe gene-matching block
      }
      bk <- blocks_i[[pick]]
      s_cons[i]  <- if (!is.na(A_Anno)  && length(bk) >= A_Anno)  bk[A_Anno]  else ""
      s_imp[i]   <- if (!is.na(A_Imp)   && length(bk) >= A_Imp)   bk[A_Imp]   else ""
      s_gene[i]  <- if (!is.na(A_Gene)  && length(bk) >= A_Gene)  bk[A_Gene]  else ""
      s_hgvsc[i] <- if (!is.na(A_Hgvsc) && length(bk) >= A_Hgvsc) bk[A_Hgvsc] else ""
    }

    # LOF: prefer gene match against Hugo_Symbol
    lof_genes <- lof_genes_col[[i]]
    lof_pcts  <- lof_pcts_col[[i]]
    if (length(lof_genes)) {
      lpick <- 1L
      if (!is.na(hugo_col[i]) && nzchar(hugo_col[i])) {
        gm <- which(lof_genes == hugo_col[i])
        if (length(gm)) lpick <- gm[1L]
      }
      l_gene[i] <- lof_genes[lpick]
      l_pct[i]  <- lof_pcts[lpick]
    }
    # NMD: same logic
    nmd_genes <- nmd_genes_col[[i]]
    nmd_pcts  <- nmd_pcts_col[[i]]
    if (length(nmd_genes)) {
      npick <- 1L
      if (!is.na(hugo_col[i]) && nzchar(hugo_col[i])) {
        gm <- which(nmd_genes == hugo_col[i])
        if (length(gm)) npick <- gm[1L]
      }
      d_gene[i] <- nmd_genes[npick]
      d_pct[i]  <- nmd_pcts[npick]
    }
  }

  # Write the resolved columns onto m IN m's row order. `joined` is in m's
  # row order because data.table's x[i, on=] preserves i's order when i is
  # the rhs of the join.
  m[, snpeff_consequence := s_cons]
  m[, snpeff_impact      := s_imp]
  m[, snpeff_gene        := s_gene]
  m[, snpeff_hgvsc       := s_hgvsc]
  m[, LOF_Gene           := l_gene]
  m[, LOF_Pct_Transcripts := l_pct]
  m[, NMD_Gene           := d_gene]
  m[, NMD_Pct_Transcripts := d_pct]

  m
}
