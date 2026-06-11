#' Convert VEP-annotated germline VCF(s) to a maftools-style MAF data.table
#'
#' @description
#' Converts VEP-annotated, single-sample germline VCFs (GATK HaplotypeCaller ->
#' CNN tranches -> Ensembl VEP, hg38) into a maftools-style table and returns it as
#' an in-memory `data.table` for downstream filtering ([gvr_filter()]) and
#' summarisation ([gvr_summary()]). In folder mode it finds every per-sample VCF,
#' converts each, and row-binds them into one combined MAF. The conversion uses base R
#' + \pkg{data.table} only - it does NOT depend on \pkg{maftools}.
#'
#' @details
#' Output and behaviour:
#' \itemize{
#'   \item Returns the final MAF `data.table`, one row per variant ALLELE
#'     (multi-allelic sites are split).
#'   \item A single most-severe transcript is chosen per allele (VEP severity ->
#'     `CANONICAL` -> `MANE_SELECT` -> transcript id).
#'   \item Columns include the MAF core fields, ALL VEP CSQ fields (read from the VCF
#'     header), and key GATK QC fields. `FILTER` is retained as a column and ALL
#'     variants (PASS and non-PASS) are kept.
#'   \item `Tumor_Seq_Allele1`/`Tumor_Seq_Allele2` are zygosity-aware (vcf2maf-style),
#'     and an optional `Genotype` column (`Tumor_Seq_Allele1/Tumor_Seq_Allele2`, e.g.
#'     `"T/C"`) is added next to the alleles.
#'   \item Each variant keeps its source sample in `Tumor_Sample_Barcode`.
#'   \item Absent values are written as the empty string `""` (not `NA`); downstream
#'     [gvr_filter()] / [gvr_summary()] treat `NA` and `""` identically as "missing".
#' }
#'
#' Processing options:
#' \itemize{
#'   \item MULTI-FILE: in folder mode, every file matching `pattern` (default
#'     `"*_NN.vcf.gz"`) is converted and row-bound.
#'   \item HGVS CLEANUP (`strip_hgvs_prefix`): strips the Ensembl feature prefix from
#'     `HGVSc`/`HGVSp` (e.g. `"ENST00000831140.1:n.1889G>A"` -> `"n.1889G>A"`).
#'   \item DEDUP (`dedup_columns`): removes duplicate-named columns ONLY when their
#'     values are byte-for-byte identical across all rows (otherwise keeps + warns).
#'   \item ABraOM (`add_abraom`): joins the Brazilian ABraOM SABE-609 allele frequency
#'     as the `ABraOM_AF` column (downloaded/cached from `abraom_url`).
#'   \item GENOTYPE-QUALITY FILTER (`min_DP`/`min_GQ`): keeps a record iff
#'     `DP > min_DP` AND `GQ > min_GQ`; mirrors
#'     `bcftools view -e 'FORMAT/DP<=X | FORMAT/GQ<=Y'`. Set either to `NULL` to
#'     disable that field; set both to `NULL` to disable the genotype filter entirely.
#'   \item GENE SUBSET (`genes`): restrict to a set of `Hugo_Symbol`s (exact,
#'     case-insensitive).
#' }
#'
#' @param folder Directory to scan in folder mode; every file matching `pattern` is
#'   converted and row-bound. Default `"."`. Ignored when `vcf_path` is supplied.
#' @param vcf_path Path to a single `.vep.vcf.gz` to convert (single-file mode).
#'   `NULL` (default) selects folder mode.
#' @param pattern Regular expression identifying per-sample VCFs in folder mode.
#'   Default `"_\\d+(\\.(vep|snp[eE]ff))?\\.vcf\\.gz$"` (matches `_01.vep.vcf.gz`, `_01.snpeff.vcf.gz`, `_01.vcf.gz`, etc.).
#' @param write_tsv Logical; if `TRUE`, also write the MAF as a TSV to `out_dir`.
#'   Default `FALSE`.
#' @param write_rds Logical; if `TRUE`, also write the MAF as an `.rds` to `out_dir`.
#'   Default `FALSE`.
#' @param write_xlsx Logical; if `TRUE`, also write the MAF as an `.xlsx` workbook
#'   (single `"MAF"` sheet) to `out_dir`. Requires the \pkg{openxlsx} package (a
#'   `Suggests` dependency); if it is not installed the export is skipped with a
#'   warning. Default `FALSE`. Note: germline MAFs can be very large; Excel handles
#'   them but the file is big and slow to open, so `write_rds`/`write_tsv` are better
#'   for large tables. Excel also caps a sheet at 1,048,576 rows.
#' @param out_dir Output directory for written TSV/RDS/XLSX. `NULL` (default) uses the
#'   input location/working directory. Only used when `write_tsv`/`write_rds`/`write_xlsx` is `TRUE`.
#' @param out_prefix Filename prefix for written outputs. `NULL` (default) derives one
#'   from the input.
#' @param chunk_size Integer; number of VCF records processed per chunk (controls peak
#'   memory and progress granularity). Default `25000L`.
#' @param ncbi_build Reference build label recorded in the MAF (`NCBI_Build`). Default
#'   `"GRCh38"`.
#' @param add_genotype Logical; if `TRUE` (default) add the `Genotype` column.
#' @param strip_hgvs_prefix Logical; if `TRUE` (default) strip the Ensembl feature
#'   prefix from `HGVSc`/`HGVSp`.
#' @param dedup_columns Logical; if `TRUE` (default) drop duplicate-named columns when
#'   byte-for-byte identical (otherwise keep and warn).
#' @param drop_empty_cols Logical; if `TRUE`, drop columns that are entirely `NA`/blank.
#'   Default `FALSE`.
#' @param add_abraom Logical; if `TRUE` (default) join the ABraOM SABE-609 allele
#'   frequency as `ABraOM_AF`.
#' @param abraom_path Path to a local ABraOM annotation file. `NULL` (default) uses an
#'   auto-managed cache (see `cache_dir`), downloading from `abraom_url` if needed.
#' @param abraom_url URL of the ABraOM SABE-609 annotated release used when
#'   `abraom_path` is `NULL`.
#' @param cache_dir Directory for the ABraOM reference cache (used only when
#'   `abraom_path` is `NULL`). `NULL` (default) uses
#'   `tools::R_user_dir("germlinevaR", "cache")`, which resolves to a
#'   platform-appropriate directory (e.g. `~/.cache/R/germlinevaR` on Linux).
#'   The directory is created on first download. Set explicitly if you prefer a
#'   custom location.
#' @param min_DP Numeric; keep only records with `DP > min_DP`. `NULL` (or `NA`)
#'   disables the depth filter. Default `10`.
#' @param min_GQ Numeric; keep only records with `GQ > min_GQ`. `NULL` (or `NA`)
#'   disables the genotype-quality filter. Default `30`.
#' @param genes Character vector of `Hugo_Symbol`s to keep (exact, case-insensitive),
#'   or `NULL` (default) to keep all genes.
#' @param vc_nonSyn Logical or character vector. Controls which
#'   `Variant_Classification` values are retained (mirroring maftools'
#'   `vc_nonSyn` argument). `FALSE` (default) keeps ALL variant classifications.
#'   `TRUE` keeps only protein-altering classes (High/Moderate VEP consequences):
#'   `"Frame_Shift_Del"`, `"Frame_Shift_Ins"`, `"Splice_Site"`,
#'   `"Translation_Start_Site"`, `"Nonsense_Mutation"`, `"Nonstop_Mutation"`,
#'   `"In_Frame_Del"`, `"In_Frame_Ins"`, `"Missense_Mutation"`. Alternatively,
#'   pass a custom character vector of classifications to keep. Rows with
#'   missing/blank `Variant_Classification` are always removed when this filter
#'   is active.
#' @param ncores Integer; number of worker processes for converting MULTIPLE input
#'   files in parallel via [parallel::mclapply()] (fork-based; Unix/macOS only).
#'   Default `1L` runs sequentially and is byte-identical to previous behaviour.
#'   Values `> 1` only help when more than one VCF is being read (each file is an
#'   independent task) and are clamped to `min(ncores, detectCores(), n_files)`. On
#'   non-fork platforms it falls back to sequential. A single file is unaffected.
#' @param verbose Logical; if `TRUE` (default) print per-file and per-chunk progress
#'   (file i/N, cumulative records, elapsed seconds).
#'
#' @return A `data.table` MAF: one row per variant allele, with MAF core columns, all
#'   VEP CSQ fields, key GATK QC fields, `Tumor_Sample_Barcode`, and (when enabled) the
#'   `Genotype` and `ABraOM_AF` columns. TSV/RDS files are written as a side
#'   effect when `write_tsv`/`write_rds` is `TRUE`.
#'
#' @seealso [gvr_filter()] to filter the returned MAF, [gvr_summary()] to summarise it.
#' @family germlinevaR
#' @author germlinevaR authors
#'
#' @examples
#' \dontrun{
#' ## Folder mode: merge ALL *_NN.vcf.gz into one MAF
#' maf <- read.gvr("/path/to/folder")
#'
#' ## Also write TSV + RDS outputs
#' maf <- read.gvr("/path/to/folder", write_tsv = TRUE, write_rds = TRUE,
#'                 out_dir = "/path/to/results")
#'
#' ## Single-file mode
#' maf <- read.gvr(vcf_path = "/path/to/SAMPLE_01.vep.vcf.gz")
#'
#' ## DP/GQ genotype filter (ON by default; mirrors
#' ##   bcftools view -e 'FORMAT/DP<=10 | FORMAT/GQ<=30')
#' maf <- read.gvr("/path/to/folder")                               # DP>10 & GQ>30
#' maf <- read.gvr("/path/to/folder", min_DP = NULL, min_GQ = NULL) # no DP/GQ filter
#' maf <- read.gvr("/path/to/folder", min_DP = 20,  min_GQ = 50)    # stricter
#'
#' ## Restrict to genes of interest (e.g. an MEN1/parathyroid panel)
#' maf <- read.gvr("/path/to/folder",
#'                 genes = c("MEN1", "RET", "CDKN1B", "CDC73", "CASR", "AIP"))
#'
#' ## Then filter freely, e.g.:
#' maf[FILTER == "PASS" & Variant_Classification == "Missense_Mutation"]
#' }
#'
#' @importFrom data.table data.table as.data.table rbindlist fread fwrite setnames
#'   setcolorder setDT set setattr tstrsplit setkey :=
#' @importFrom stats setNames
#' @importFrom utils download.file
#' @importFrom openxlsx createWorkbook
#' @export
read.gvr <- function(folder = ".",
                           vcf_path   = NULL,
                           pattern    = "_\\d+(\\.(vep|snp[eE]ff))?\\.vcf\\.gz$",   # v2/v9: _01.vep.vcf.gz, _01.snpeff.vcf.gz, _01.vcf.gz, ...
                           write_tsv  = FALSE,
                           write_rds  = FALSE,
                           write_xlsx = FALSE,   # v6: also write the MAF as .xlsx (one "MAF" sheet)
                           out_dir    = NULL,
                           out_prefix = NULL,
                           chunk_size = 25000L,
                           ncbi_build = "GRCh38",
                           add_genotype      = TRUE,   # v2
                           strip_hgvs_prefix = TRUE,   # v2
                           dedup_columns     = TRUE,   # v2
                           drop_empty_cols   = FALSE,  # v3: drop all-NA/blank cols
                           add_abraom        = TRUE,   # v3: join ABraOM SABE609 AF
                           abraom_path       = NULL,   # v3: local file; NULL=auto cache
                           abraom_url        = "https://abraom.ib.usp.br/download/ABRaOM_60+_SABE_609_exomes_annotated.gz",
                           cache_dir         = NULL,   # v7: ABraOM cache dir; NULL = tools::R_user_dir()
                           min_DP            = 10,     # v4: keep only DP > min_DP (NULL = no DP filter)
                           min_GQ            = 30,     # v4: keep only GQ > min_GQ (NULL = no GQ filter)
                           genes             = NULL,   # v4: keep only these Hugo_Symbols
                           vc_nonSyn         = FALSE,  # v8: keep only protein-altering Variant_Classification
                           ncores            = 1L,     # v6: parallel files (>1 forks mclapply; 1 = sequential, default)
                           verbose    = TRUE) {
  # ===========================================================================
  # Nested helper: .read_gvr_locate_sibling() + auto-source sibling
  # ---------------------------------------------------------------------------
  # In standalone mode (source()'d from a directory), read.gvr() tries to
  # auto-source its sibling file read.gvr.snpeff.R so that SnpEff dispatch
  # works without manual setup. In a package, all R/*.R files are loaded
  # into the namespace automatically, so the auto-source is skipped.
  # ===========================================================================
  .read_gvr_locate_sibling <- function() {
    # 1) Use the active source file's directory if R can tell us where THIS
    #    script lives (works for source() and Rscript). Recent R: getSrcFilename
    #    on a function defined HERE returns this file's path.
    this_dir <- tryCatch({
      fn <- function() {}
      fp <- utils::getSrcFilename(fn, full.names = TRUE)
      if (length(fp) && nzchar(fp)) normalizePath(dirname(fp), mustWork = FALSE) else NULL
    }, error = function(e) NULL)
    # 2) Fallback: look in current working directory.
    if (is.null(this_dir) || !nzchar(this_dir)) this_dir <- getwd()
    candidate <- file.path(this_dir, "read.gvr.snpeff.R")
    if (file.exists(candidate)) return(candidate)
    NULL
  }

  # Auto-source sibling if read.gvr.snpeff is not already available
  # (in a package, it always is; in standalone mode, this provides it)
  if (!exists("read.gvr.snpeff", mode = "function", inherits = TRUE)) {
    sib <- .read_gvr_locate_sibling()
    if (!is.null(sib)) {
      tryCatch(
        source(sib, local = FALSE, chdir = FALSE),
        error = function(e) warning(sprintf(
          "read.gvr: failed to source sibling '%s': %s", sib, conditionMessage(e)))
      )
    }
  }

  # ===========================================================================
  # Nested helper: .detect_annotator()
  # ---------------------------------------------------------------------------
  # Was previously top-level (defined in read.gvr.snpeff.R and source()d via
  # the sibling auto-loader above). Nested here so read.gvr() is self-
  # contained for package hygiene. The sibling read.gvr.snpeff.R file also
  # nests its own copy. Both are byte-identical pure functions.
  # ===========================================================================

  .detect_annotator <- function(path) {
    con <- gzfile(path, "r")
    on.exit(close(con))
    csq_found <- FALSE
    ann_found <- FALSE
    info_tags <- character(0L)
    repeat {
      line <- readLines(con, n = 1L, warn = FALSE)
      if (length(line) == 0L) break              # EOF before any record line
      if (!startsWith(line, "##")) break          # hit body / #CHROM header
      if (startsWith(line, "##INFO=<ID=")) {
        tag <- sub("^##INFO=<ID=([^,>]+).*$", "\\1", line)
        info_tags <- c(info_tags, tag)
        if (identical(tag, "CSQ")) csq_found <- TRUE
        if (identical(tag, "ANN")) ann_found <- TRUE
      }
    }
    if (csq_found && ann_found) {
      warning(sprintf(
        "read.gvr: '%s' has both VEP CSQ and SnpEff ANN INFO tags. Using VEP (CSQ takes priority).",
        basename(path)),
        call. = FALSE)
      res <- "vep"
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


  # ==========================================================================
  # v4 setup: genotype-quality filter flags
  #   min_DP / min_GQ are thresholds; keep a record iff DP > min_DP AND GQ > min_GQ.
  #   Set either to NULL (or NA) to disable filtering on that field; set BOTH to
  #   NULL/NA to disable the genotype filter entirely (equivalent to the old
  #   filter_genotype = FALSE). Mirrors bcftools -e 'FORMAT/DP<=X | FORMAT/GQ<=Y'.
  # ==========================================================================
  filter_dp <- !is.null(min_DP) && !is.na(min_DP)
  filter_gq <- !is.null(min_GQ) && !is.na(min_GQ)
  if (filter_dp && !is.numeric(min_DP)) stop("min_DP must be numeric or NULL.")
  if (filter_gq && !is.numeric(min_GQ)) stop("min_GQ must be numeric or NULL.")

  # ==========================================================================
  # 0. Locate the input VCF(s)
  #    - folder mode  : process EVERY file matching `pattern` (multi-file merge).
  #    - vcf_path mode : process exactly that one file (backward compatible).
  # ==========================================================================
  if (is.null(vcf_path)) {
    if (!dir.exists(folder))
      stop(sprintf("Folder does not exist: %s", folder))
    hits <- list.files(folder, pattern = pattern, full.names = TRUE,
                       ignore.case = FALSE)
    if (length(hits) == 0L)
      stop(sprintf("No file matching '%s' found in: %s", pattern, folder))
    # deterministic order: by filename ascending so _01 precedes _02 precedes _03 ...
    hits <- hits[order(basename(hits))]
    vcf_paths <- hits
  } else {
    if (!file.exists(vcf_path))
      stop(sprintf("vcf_path does not exist: %s", vcf_path))
    vcf_paths <- vcf_path
  }

  if (verbose) {
    message(sprintf("Found %d file(s):", length(vcf_paths)))
    for (p in vcf_paths) message("  - ", basename(p))
  }
  # ==========================================================================
  # 0b. Auto-detect annotator per file + dispatch
  #     read.gvr() handles VEP-annotated VCFs (its original role); SnpEff
  #     batches are delegated to read.gvr.snpeff(). In a package, the sibling
  #     is always available; in standalone mode, the auto-source above loads it.
  # ==========================================================================
  detected <- vapply(vcf_paths, .detect_annotator, character(1L))
  na_mask  <- is.na(detected)
  if (any(na_mask)) {
    ## Surface up to 7 INFO tags from the first un-annotated file so the user
    ## can see what *was* annotated (e.g. raw GATK output vs. exotic annotator).
    ## .detect_annotator() returns NA with attr 'info_tags' attached.
    .unann_first <- vcf_paths[na_mask][1L]
    .tags_first  <- attr(.detect_annotator(.unann_first), "info_tags")
    .tag_str <- if (length(.tags_first))
      sprintf("INFO tags found in '%s': %s.\n",
              basename(.unann_first),
              paste(utils::head(.tags_first, 7L), collapse = ", "))
    else ""
    stop(sprintf(
      "read.gvr: file(s) have no VEP CSQ or SnpEff ANN INFO tag:\n%s\n%s%s",
      paste(sprintf("  - %s", basename(vcf_paths[na_mask])), collapse = "\n"),
      .tag_str,
      paste0("Annotate these VCFs with Ensembl VEP (writes ##INFO=<ID=CSQ>) ",
             "or SnpEff (writes ##INFO=<ID=ANN>) before passing to read.gvr() ",
             "/ read.gvr.snpeff().")
    ), call. = FALSE)
  }
  uniq_ann <- unique(detected)
  if (length(uniq_ann) > 1L) {
    by_ann <- split(basename(vcf_paths), detected)
    parts  <- vapply(names(by_ann), function(a)
      sprintf("  %s: %s", a, paste(by_ann[[a]], collapse = ", ")),
      character(1L))
    stop(sprintf("read.gvr: mixed annotators in batch -- refusing to merge:\n%s",
                 paste(parts, collapse = "\n")), call. = FALSE)
  }
  if (uniq_ann == "snpeff") {
    if (!exists("read.gvr.snpeff", mode = "function", inherits = TRUE)) {
      stop("read.gvr: detected SnpEff VCF(s) but 'read.gvr.snpeff' is not ",
           "defined. Source 'read.gvr.snpeff.R' (placed next to this file) ",
           "first.", call. = FALSE)
    }
    if (verbose) message("read.gvr: SnpEff-annotated input detected; delegating to read.gvr.snpeff().")
    return(read.gvr.snpeff(
      folder            = folder,
      vcf_path          = vcf_path,
      pattern           = pattern,
      write_tsv         = write_tsv,
      write_rds         = write_rds,
      write_xlsx        = write_xlsx,
      out_dir           = out_dir,
      out_prefix        = out_prefix,
      chunk_size        = chunk_size,
      ncbi_build        = ncbi_build,
      add_genotype      = add_genotype,
      strip_hgvs_prefix = strip_hgvs_prefix,
      dedup_columns     = dedup_columns,
      drop_empty_cols   = drop_empty_cols,
      add_abraom        = add_abraom,
      abraom_path       = abraom_path,
      abraom_url        = abraom_url,
      cache_dir         = cache_dir,
      min_DP            = min_DP,
      min_GQ            = min_GQ,
      genes             = genes,
      vc_nonSyn         = vc_nonSyn,
      ncores            = ncores,
      verbose           = verbose
    ))
  }
  # else uniq_ann == "vep": continue with the existing VEP body below.


  # ==========================================================================
  # 1. Local helper definitions (kept inside the function => fully self-contained)
  #    [UNCHANGED CORE ENGINE from v1]
  # ==========================================================================

  ## 1a. Ensembl consequence severity priority (vcf2maf-style; lower = more severe)
  effect_priority <- c(
    "transcript_ablation"=1,"exon_loss_variant"=2,"splice_donor_variant"=3,
    "splice_acceptor_variant"=3,"stop_gained"=4,"frameshift_variant"=5,
    "stop_lost"=6,"start_lost"=7,"initiator_codon_variant"=8,
    "transcript_amplification"=9,"protein_altering_variant"=10,
    "missense_variant"=11,"conservative_missense_variant"=11,
    "rare_amino_acid_variant"=11,"incomplete_terminal_codon_variant"=14,
    "splice_region_variant"=13,"splice_donor_5th_base_variant"=13,
    "splice_donor_region_variant"=13,"splice_polypyrimidine_tract_variant"=13,
    "stop_retained_variant"=15,"synonymous_variant"=15,"coding_sequence_variant"=16,
    "mature_miRNA_variant"=17,"5_prime_UTR_variant"=18,
    "5_prime_UTR_premature_start_codon_gain_variant"=18,"3_prime_UTR_variant"=19,
    "non_coding_transcript_exon_variant"=20,"non_coding_exon_variant"=20,
    "intron_variant"=21,"non_coding_transcript_variant"=22,"nc_transcript_variant"=22,
    "NMD_transcript_variant"=23,"upstream_gene_variant"=24,"downstream_gene_variant"=25,
    "TFBS_ablation"=26,"TFBS_amplification"=27,"TF_binding_site_variant"=28,
    "regulatory_region_ablation"=29,"regulatory_region_amplification"=30,
    "regulatory_region_variant"=31,"feature_elongation"=32,"feature_truncation"=33,
    "intergenic_variant"=34
  )

  # B2 + A2(memoized): single pass that returns BOTH the most-severe term AND its rank.
  # The computation (strsplit on '&' + effect_priority lookup) is a PURE function of the
  # raw Consequence string, drawn from only ~107 distinct values across ~148k allele-rows,
  # so results are cached in a per-call environment keyed by the exact string. On a hit the
  # stored list is returned; on a miss the value is computed by the unchanged inner function
  # and stored. Output is byte-identical to the uncached path by construction (same function,
  # same input). The thin most_severe_term()/consequence_rank() wrappers route through the
  # cached resolver, so all three call sites share one cache.
  .mstr_cache <- new.env(parent = emptyenv())
  .most_severe_term_and_rank_uncached <- function(consequence) {
    if (is.na(consequence) || consequence == "")
      return(list(term = NA_character_, rank = 999L))
    terms <- strsplit(consequence, "&", fixed = TRUE)[[1]]
    pr <- effect_priority[terms]; pr[is.na(pr)] <- 999L
    wm <- which.min(pr)
    list(term = terms[wm], rank = as.integer(pr[wm]))
  }
  most_severe_term_and_rank <- function(consequence) {
    key <- if (is.na(consequence)) "\001NA" else consequence
    hit <- .mstr_cache[[key]]
    if (!is.null(hit)) return(hit)
    val <- .most_severe_term_and_rank_uncached(consequence)
    assign(key, val, envir = .mstr_cache)
    val
  }
  most_severe_term <- function(consequence) most_severe_term_and_rank(consequence)$term
  consequence_rank <- function(consequence) most_severe_term_and_rank(consequence)$rank

  ## 1b. VEP consequence -> MAF Variant_Classification  (A2: memoized)
  ## var_class is a PURE function of (top_term, var_type); both are small-cardinality
  ## (19 distinct classes; a handful of var_types), so cache on the composite key.
  ## Inner function unchanged; cached value byte-identical.
  .v2m_cache <- new.env(parent = emptyenv())
  .vep_to_maf_class_uncached <- function(term, var_type) {
    m <- switch(term,
      "splice_acceptor_variant"="Splice_Site","splice_donor_variant"="Splice_Site",
      "transcript_ablation"="Splice_Site","exon_loss_variant"="Splice_Site",
      "stop_gained"="Nonsense_Mutation","stop_lost"="Nonstop_Mutation",
      "start_lost"="Translation_Start_Site","initiator_codon_variant"="Translation_Start_Site",
      "missense_variant"="Missense_Mutation","conservative_missense_variant"="Missense_Mutation",
      "rare_amino_acid_variant"="Missense_Mutation","protein_altering_variant"="Missense_Mutation",
      "transcript_amplification"="Intron","splice_region_variant"="Splice_Region",
      "splice_donor_5th_base_variant"="Splice_Region","splice_donor_region_variant"="Splice_Region",
      "splice_polypyrimidine_tract_variant"="Splice_Region","stop_retained_variant"="Silent",
      "synonymous_variant"="Silent","incomplete_terminal_codon_variant"="Silent",
      "coding_sequence_variant"="Missense_Mutation","mature_miRNA_variant"="RNA",
      "5_prime_UTR_variant"="5'UTR","5_prime_UTR_premature_start_codon_gain_variant"="5'UTR",
      "3_prime_UTR_variant"="3'UTR","non_coding_transcript_exon_variant"="RNA",
      "non_coding_exon_variant"="RNA","non_coding_transcript_variant"="RNA",
      "nc_transcript_variant"="RNA","NMD_transcript_variant"="Silent",
      "intron_variant"="Intron","upstream_gene_variant"="5'Flank",
      "downstream_gene_variant"="3'Flank","TFBS_ablation"="Targeted_Region",
      "TFBS_amplification"="Targeted_Region","TF_binding_site_variant"="IGR",
      "regulatory_region_ablation"="Targeted_Region","regulatory_region_amplification"="Targeted_Region",
      "regulatory_region_variant"="IGR","feature_elongation"="Targeted_Region",
      "feature_truncation"="Targeted_Region","intergenic_variant"="IGR",
      "Targeted_Region"
    )
    if (term == "frameshift_variant") m <- if (var_type == "DEL") "Frame_Shift_Del" else "Frame_Shift_Ins"
    if (term == "inframe_insertion")            m <- "In_Frame_Ins"
    if (term == "inframe_deletion")             m <- "In_Frame_Del"
    if (term == "disruptive_inframe_insertion") m <- "In_Frame_Ins"
    if (term == "disruptive_inframe_deletion")  m <- "In_Frame_Del"
    m
  }
  vep_to_maf_class <- function(term, var_type) {
    key <- paste0(if (is.na(term)) "\001NA" else term, "\002",
                  if (is.na(var_type)) "\001NA" else var_type)
    hit <- .v2m_cache[[key]]
    if (!is.null(hit)) return(hit)
    val <- .vep_to_maf_class_uncached(term, var_type)
    assign(key, val, envir = .v2m_cache)
    val
  }

  ## 1c. Amino-acid 3->1 and HGVSp_Short
  aa3to1 <- c(Ala="A",Arg="R",Asn="N",Asp="D",Cys="C",Gln="Q",Glu="E",Gly="G",
              His="H",Ile="I",Leu="L",Lys="K",Met="M",Phe="F",Pro="P",Ser="S",
              Thr="T",Trp="W",Tyr="Y",Val="V",Ter="*",Sec="U",Pyl="O",Asx="B",
              Glx="Z",Xaa="X",Xle="J")
  url_decode <- function(x) {
    if (is.na(x) || x == "") return(x)
    x <- gsub("%3D","=",x,fixed=TRUE); x <- gsub("%3B",";",x,fixed=TRUE)
    x <- gsub("%2C",",",x,fixed=TRUE); x <- gsub("%3A",":",x,fixed=TRUE); x
  }
  ## C2: single-pass 3-letter -> 1-letter amino-acid substitution. The 27 codes are
  ## mutually non-overlapping (no code is a substring of another) and the 1-letter
  ## outputs can never re-match a 3-letter code, so a single regex alternation pass is
  ## order-independent and byte-identical to the original 27 sequential gsub() calls
  ## (verified: 0 diffs over all 28,686 unique real HGVSp strings). aa3to1_pat is built
  ## once (all codes are literal alpha, regex-safe).
  aa3to1_pat <- paste(names(aa3to1), collapse = "|")
  make_hgvsp_short <- function(hgvsp) {
    if (is.na(hgvsp) || hgvsp == "") return("")
    s <- url_decode(hgvsp); s <- sub("^[^:]*:", "", s); s <- gsub("[()]", "", s)
    m <- gregexpr(aa3to1_pat, s, perl = TRUE)
    regmatches(s, m) <- lapply(regmatches(s, m), function(hit) unname(aa3to1[hit]))
    s
  }

  ## 1c-v2. Strip the leading "feature_id:" prefix from an HGVSc/HGVSp string.
  ##        "ENST00000831140.1:n.1889G>A" -> "n.1889G>A". Empty/NA pass through.
  ##        (URL-decoding is applied upstream; this only removes up to first ':'.)
  strip_feature_prefix <- function(x) {
    if (is.na(x) || x == "") return(x)
    sub("^[^:]*:", "", x)
  }

  ## 1d. Candidate VEP CSQ-allele encodings (priority order: anchor, LCP, bidir, raw)
  vep_allele_candidates <- function(ref, alt) {
    if (alt == "*") return("*")
    if (nchar(ref) == 1 && nchar(alt) == 1) return(alt)            # SNV
    rl <- nchar(ref); al <- nchar(alt); cand <- character(0)
    # (a) anchor-trim: drop exactly one leading base
    r1 <- substr(ref,2,rl); a1 <- substr(alt,2,al)
    cand <- c(cand, if (nchar(a1)==0) "-" else if (nchar(r1)==0) a1 else a1)
    # (b) longest-common-prefix trim
    k <- 0L; mx <- min(rl,al)
    while (k < mx && substr(ref,k+1,k+1) == substr(alt,k+1,k+1)) k <- k+1L
    rp <- substr(ref,k+1,rl); ap <- substr(alt,k+1,al)
    cand <- c(cand, if (nchar(ap)==0) "-" else if (nchar(rp)==0) ap else ap)
    # (c) bidirectional trim (repeat-aware): common suffix then common prefix
    rr <- ref; aa <- alt
    while (nchar(rr)>0 && nchar(aa)>0 &&
           substr(rr,nchar(rr),nchar(rr)) == substr(aa,nchar(aa),nchar(aa))) {
      rr <- substr(rr,1,nchar(rr)-1); aa <- substr(aa,1,nchar(aa)-1)
    }
    j <- 0L; mx2 <- min(nchar(rr),nchar(aa))
    while (j < mx2 && substr(rr,j+1,j+1) == substr(aa,j+1,j+1)) j <- j+1L
    rb <- substr(rr,j+1,nchar(rr)); ab <- substr(aa,j+1,nchar(aa))
    cand <- c(cand, if (nchar(ab)==0) "-" else if (nchar(rb)==0) ab else ab)
    # (d) raw ALT
    cand <- c(cand, alt)
    unique(cand[cand != ""])
  }

  ## 1e. MAF coordinate + allele conversion for one REF/ALT pair
  maf_coords <- function(pos, ref, alt) {
    pos <- as.integer(pos)
    if (alt == "*") return(list(var_type="DEL", start=pos, end=pos,
                                ref_allele=ref, tum_allele2="*"))
    rl <- nchar(ref); al <- nchar(alt)
    if (rl == al && rl == 1)
      return(list(var_type="SNP", start=pos, end=pos, ref_allele=ref, tum_allele2=alt))
    if (rl == al && rl > 1) {
      vt <- switch(as.character(rl), "2"="DNP","3"="TNP","ONP")
      return(list(var_type=vt, start=pos, end=pos+rl-1, ref_allele=ref, tum_allele2=alt))
    }
    if (al > rl) {                                   # insertion
      ins <- substr(alt, rl+1, al)
      return(list(var_type="INS", start=pos+rl-1, end=pos+rl,
                  ref_allele="-", tum_allele2=ins))
    } else {                                         # deletion
      del <- substr(ref, al+1, rl)
      return(list(var_type="DEL", start=pos+al, end=pos+rl-1,
                  ref_allele=del, tum_allele2="-"))
    }
  }

  ## 1f. INFO field accessor
  info_get <- function(info_str, key) {
    pat <- paste0("(^|;)", key, "=([^;]*)")
    m <- regmatches(info_str, regexec(pat, info_str))[[1]]
    if (length(m) >= 3) m[3] else NA_character_
  }
  ## C1: parse an INFO string ONCE into a named character vector, then look keys up
  ## by name. Equivalent to calling info_get() per key (verified 0 mismatches over
  ## 35,000 real comparisons): only "key=value" tokens are kept, so a missing key OR a
  ## bare flag returns NA exactly as info_get() does. ~5.7x faster than 6+ regex scans.
  info_parse <- function(info_str) {
    if (is.na(info_str) || info_str == "") return(setNames(character(0), character(0)))
    kv  <- strsplit(info_str, ";", fixed = TRUE)[[1]]
    eq  <- regexpr("=", kv, fixed = TRUE)
    has <- eq > 0L
    keys <- substr(kv, 1L, eq - 1L)            # only meaningful where has==TRUE
    vals <- substr(kv, eq + 1L, nchar(kv))
    setNames(vals[has], keys[has])
  }

  ## 1g. Zygosity-aware genotype allele CODES for a split row (Allele2 = this ALT)
  gt_codes_for_alt <- function(gt, ai) {
    ai <- as.integer(ai)
    gidx <- suppressWarnings(as.integer(strsplit(gt, "[/|]")[[1]]))
    gidx <- gidx[!is.na(gidx)]
    if (length(gidx) == 0L) return(list(c1 = NA_integer_, c2 = ai))
    others <- gidx[gidx != ai]
    if (length(others) == 0L) return(list(c1 = ai, c2 = ai))   # hom-alt
    partner <- if (0L %in% others) 0L else others[1]
    list(c1 = partner, c2 = ai)
  }

  ## 1h. Header parsers
  get_csq_fields <- function(path) {
    con <- gzfile(path, "r"); on.exit(close(con))
    repeat {
      line <- readLines(con, n = 1L)
      if (length(line) == 0L) stop("CSQ header not found")
      if (startsWith(line, "##INFO=<ID=CSQ")) {
        fmt <- sub('.*Format: ', '', line); fmt <- sub('">$', '', fmt)
        return(strsplit(fmt, "|", fixed = TRUE)[[1]])
      }
      if (startsWith(line, "#CHROM")) stop("Reached #CHROM before CSQ header")
    }
  }
  get_sample_name <- function(path) {
    con <- gzfile(path, "r"); on.exit(close(con))
    repeat {
      line <- readLines(con, n = 1L)
      if (length(line) == 0L) stop("#CHROM not found")
      if (startsWith(line, "#CHROM")) {
        cols <- strsplit(line, "\t", fixed = TRUE)[[1]]
        return(cols[length(cols)])
      }
    }
  }

  ## 1i. Convert one chunk (data.table of raw VCF rows) -> MAF rows
  ## C3 (conservative): the per-record control flow is UNCHANGED (the CSQ-block -> ALT
  ## assignment and per-ALT selection are genuinely irregular and were measured to be
  ## slower, not faster, when forced through flat strsplit + regroup). The one safe,
  ## byte-identical structural change kept here is pre-extracting the 10 VCF columns to
  ## plain atomic vectors ONCE per chunk, so the hot loop indexes vectors (chrom_v[r])
  ## instead of doing data.table `$`/is.factor dispatch on every field of every record
  ## (~4.7x faster for the indexing step in isolation). Output is identical.
  convert_chunk <- function(dt, csq_fields, sample_name) {
    n_csq <- length(csq_fields)
    ci <- function(name) match(name, csq_fields)
    P_Allele <- ci("Allele"); P_Cons <- ci("Consequence"); P_SYMBOL <- ci("SYMBOL")
    P_Feature <- ci("Feature"); P_HGVSc <- ci("HGVSc"); P_HGVSp <- ci("HGVSp")
    P_Existing <- ci("Existing_variation"); P_CANONICAL <- ci("CANONICAL")
    P_MANESEL <- ci("MANE_SELECT")

    # C3: hoist columns to plain vectors once (kills per-record `$`/is.factor overhead)
    CHROM_v <- dt$CHROM; POS_v <- dt$POS; ID_v <- dt$ID; REF_v <- dt$REF; ALT_v <- dt$ALT
    QUAL_v <- dt$QUAL; FILTER_v <- dt$FILTER; INFO_v <- dt$INFO
    FORMAT_v <- dt$FORMAT; SAMPLE_v <- dt$SAMPLE
    nrow_dt <- nrow(dt)

    # O2 (PERF): split the ALT column ONCE for the whole chunk instead of calling
    # strsplit(altf, ",") per record. base strsplit() over a character vector returns
    # a list-of-vectors; alts_all[[r]] is byte-identical to the old per-record
    # strsplit(ALT_v[r], ",", fixed=TRUE)[[1]]. ~30x faster on the isolated split
    # (multi-ALT is rare, so per-record calls were almost pure call-overhead).
    alts_all <- strsplit(ALT_v, ",", fixed = TRUE)

    # O1 (PERF): the per-record loop split FORMAT and SAMPLE on ":" every iteration and
    # then looked GT/AD/DP/GQ up BY NAME. In real VEP/GATK germline VCFs the FORMAT
    # string is almost always constant within a chunk (e.g. "GT:AD:DP:GQ:PL"). When it
    # IS constant we resolve the GT/AD/DP/GQ field POSITIONS once and pull them with a
    # single vectorized tstrsplit(SAMPLE, ":") over the whole column (~4x faster). If
    # FORMAT varies anywhere in the chunk we fall back to the exact per-record path, so
    # output is byte-identical in every case (verified by the validation gate).
    fmt_u <- unique(FORMAT_v)
    fmt_constant <- length(fmt_u) == 1L && !is.na(fmt_u[1L])
    if (fmt_constant) {
      fmt_keys0 <- strsplit(fmt_u[1L], ":", fixed = TRUE)[[1]]
      pos_GT <- match("GT", fmt_keys0); pos_AD <- match("AD", fmt_keys0)
      pos_DP <- match("DP", fmt_keys0); pos_GQ <- match("GQ", fmt_keys0)
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
    if (filter_dp || filter_gq) {
      keep <- rep(TRUE, nrow_dt)
      if (fmt_constant) {
        # Fast path: DP/GQ already extracted as chunk-level vectors
        if (filter_dp && !is.null(DP_col)) {
          dp_num <- suppressWarnings(as.numeric(DP_col))
          keep <- keep & (is.na(dp_num) | dp_num > min_DP)
        }
        if (filter_gq && !is.null(GQ_col)) {
          gq_num <- suppressWarnings(as.numeric(GQ_col))
          keep <- keep & (is.na(gq_num) | gq_num > min_GQ)
        }
      } else {
        # Slow path: FORMAT varies -- extract DP/GQ per record (same logic as the loop)
        for (r in seq_len(nrow_dt)) {
          fmt_keys <- strsplit(FORMAT_v[r], ":", fixed = TRUE)[[1]]
          smp_vals <- strsplit(SAMPLE_v[r], ":", fixed = TRUE)[[1]]
          names(smp_vals) <- fmt_keys[seq_along(smp_vals)]
          if (filter_dp) {
            sdp <- if ("DP" %in% names(smp_vals)) smp_vals[["DP"]] else NA_character_
            dp_num <- suppressWarnings(as.numeric(sdp))
            if (!is.na(dp_num) && dp_num <= min_DP) keep[r] <- FALSE
          }
          if (filter_gq) {
            sgq <- if ("GQ" %in% names(smp_vals)) smp_vals[["GQ"]] else NA_character_
            gq_num <- suppressWarnings(as.numeric(sgq))
            if (!is.na(gq_num) && gq_num <= min_GQ) keep[r] <- FALSE
          }
        }
      }
      n_dropped_dpgq <- sum(!keep)
      if (n_dropped_dpgq > 0L) {
        # Subset ALL chunk-level vectors to surviving records only
        idx <- which(keep)
        CHROM_v <- CHROM_v[idx]; POS_v <- POS_v[idx]; ID_v <- ID_v[idx]
        REF_v <- REF_v[idx]; ALT_v <- ALT_v[idx]; QUAL_v <- QUAL_v[idx]
        FILTER_v <- FILTER_v[idx]; INFO_v <- INFO_v[idx]
        FORMAT_v <- FORMAT_v[idx]; SAMPLE_v <- SAMPLE_v[idx]
        alts_all <- alts_all[idx]
        if (fmt_constant) {
          GT_col  <- if (!is.null(GT_col))  GT_col[idx]  else NULL
          AD_col  <- if (!is.null(AD_col))  AD_col[idx]  else NULL
          DP_col  <- if (!is.null(DP_col))  DP_col[idx]  else NULL
          GQ_col  <- if (!is.null(GQ_col))  GQ_col[idx]  else NULL
        }
        nrow_dt <- length(idx)
      }
    }


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
        gene_hit <- grepl(pat, INFO_v, ignore.case = TRUE)
        n_dropped_gene_rough <- sum(!gene_hit)
        if (n_dropped_gene_rough > 0L) {
          idx <- which(gene_hit)
          CHROM_v <- CHROM_v[idx]; POS_v <- POS_v[idx]; ID_v <- ID_v[idx]
          REF_v <- REF_v[idx]; ALT_v <- ALT_v[idx]; QUAL_v <- QUAL_v[idx]
          FILTER_v <- FILTER_v[idx]; INFO_v <- INFO_v[idx]
          FORMAT_v <- FORMAT_v[idx]; SAMPLE_v <- SAMPLE_v[idx]
          alts_all <- alts_all[idx]
          if (fmt_constant) {
            GT_col  <- if (!is.null(GT_col))  GT_col[idx]  else NULL
            AD_col  <- if (!is.null(AD_col))  AD_col[idx]  else NULL
            DP_col  <- if (!is.null(DP_col))  DP_col[idx]  else NULL
            GQ_col  <- if (!is.null(GQ_col))  GQ_col[idx]  else NULL
          }
          nrow_dt <- length(idx)
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
      vc_hit <- grepl(pat_vc, INFO_v, ignore.case = FALSE)  # VEP terms are lowercase
      n_dropped_vc_rough <- sum(!vc_hit)
      if (n_dropped_vc_rough > 0L) {
        idx <- which(vc_hit)
        CHROM_v <- CHROM_v[idx]; POS_v <- POS_v[idx]; ID_v <- ID_v[idx]
        REF_v <- REF_v[idx]; ALT_v <- ALT_v[idx]; QUAL_v <- QUAL_v[idx]
        FILTER_v <- FILTER_v[idx]; INFO_v <- INFO_v[idx]
        FORMAT_v <- FORMAT_v[idx]; SAMPLE_v <- SAMPLE_v[idx]
        alts_all <- alts_all[idx]
        if (fmt_constant) {
          GT_col  <- if (!is.null(GT_col))  GT_col[idx]  else NULL
          AD_col  <- if (!is.null(AD_col))  AD_col[idx]  else NULL
          DP_col  <- if (!is.null(DP_col))  DP_col[idx]  else NULL
          GQ_col  <- if (!is.null(GQ_col))  GQ_col[idx]  else NULL
        }
        nrow_dt <- length(idx)
      }
    }

    out <- vector("list", nrow_dt); oi <- 0L
    n_dropped <- n_dropped_dpgq + n_dropped_gene_rough + n_dropped_vc_rough
    for (r in seq_len(nrow_dt)) {
      chrom <- CHROM_v[r]; pos <- POS_v[r]; vid <- ID_v[r]
      ref <- REF_v[r]; altf <- ALT_v[r]; qual <- QUAL_v[r]; filt <- FILTER_v[r]
      info <- INFO_v[r]; fmt <- FORMAT_v[r]; smp <- SAMPLE_v[r]
      alts <- alts_all[[r]]                          # O2: pre-split once per chunk

      ip <- info_parse(info)                       # C1: parse INFO once, then index
      info_DP <- unname(ip["DP"]); info_AC <- unname(ip["AC"])
      info_AF <- unname(ip["AF"]); info_MQ <- unname(ip["MQ"])
      info_QD <- unname(ip["QD"]); info_CNN <- unname(ip["CNN_1D"])
      ac_vec <- if (!is.na(info_AC)) strsplit(info_AC,",",fixed=TRUE)[[1]] else NA
      af_vec <- if (!is.na(info_AF)) strsplit(info_AF,",",fixed=TRUE)[[1]] else NA

      csq_raw <- unname(ip["CSQ"])
      csq_blocks <- if (!is.na(csq_raw)) strsplit(csq_raw,",",fixed=TRUE)[[1]] else character(0)
      # A3 (PERF): split ALL CSQ blocks in ONE vectorized strsplit() call rather than
      # one strsplit() per block in an lapply. base strsplit() drops trailing empty
      # fields, so each result is padded back to n_csq with `length(f) <- n_csq`
      # (NA fill). Byte-identical to the per-block path (verified on 186k real blocks
      # + edge cases); ~22% faster on the isolated split. fixed=TRUE, no new deps.
      block_fields <- lapply(strsplit(csq_blocks, "|", fixed = TRUE), function(f) {
        length(f) <- n_csq; f
      })
      block_allele <- vapply(block_fields, function(f) {
        a <- f[P_Allele]; if (is.na(a)) "" else a }, character(1))

      # rank-aware CSQ block -> ALT owner assignment
      cand_by_alt <- lapply(alts, function(a) vep_allele_candidates(ref, a))
      block_owner <- integer(length(block_allele))
      if (length(block_allele) > 0L) {
        for (bi in seq_along(block_allele)) {
          bs <- block_allele[bi]; best_alt <- 0L; best_rank <- .Machine$integer.max
          for (ai2 in seq_along(alts)) {
            rk <- match(bs, cand_by_alt[[ai2]])
            if (!is.na(rk) && rk < best_rank) { best_rank <- rk; best_alt <- ai2 }
          }
          block_owner[bi] <- best_alt
        }
      }

      # O1: GT/AD/DP/GQ extraction. Fast path uses the chunk-level position split when
      # FORMAT is constant; otherwise the exact original per-record name-based path.
      # In BOTH paths, a missing field resolves to the SAME defaults the original used
      # (gt="./.", others NA_character_). For the fast path, NA from tstrsplit (SAMPLE
      # shorter than the position) is treated as "field absent" -> default, which is
      # identical to the original `seq_along(smp_vals)` name-assignment behaviour.
      if (fmt_constant) {
        gt0  <- if (is.null(GT_col)) NA_character_ else GT_col[r]
        gt   <- if (is.na(gt0)) "./." else gt0
        ad   <- if (is.null(AD_col)) NA_character_ else AD_col[r]
        sdp  <- if (is.null(DP_col)) NA_character_ else DP_col[r]
        gq   <- if (is.null(GQ_col)) NA_character_ else GQ_col[r]
      } else {
        fmt_keys <- strsplit(fmt,":",fixed=TRUE)[[1]]
        smp_vals <- strsplit(smp,":",fixed=TRUE)[[1]]
        names(smp_vals) <- fmt_keys[seq_along(smp_vals)]
        gt  <- if ("GT" %in% names(smp_vals)) smp_vals[["GT"]] else "./."
        ad  <- if ("AD" %in% names(smp_vals)) smp_vals[["AD"]] else NA_character_
        sdp <- if ("DP" %in% names(smp_vals)) smp_vals[["DP"]] else NA_character_
        gq  <- if ("GQ" %in% names(smp_vals)) smp_vals[["GQ"]] else NA_character_
      }
      ad_vec <- if (!is.na(ad)) strsplit(ad,",",fixed=TRUE)[[1]] else NA

      for (ai in seq_along(alts)) {
        alt <- alts[ai]
        sel <- which(block_owner == ai)
        coords <- maf_coords(pos, ref, alt); vt <- coords$var_type

        if (length(sel) > 0L) {
          ranks <- vapply(sel, function(k) consequence_rank(block_fields[[k]][P_Cons]), integer(1))
          canon <- vapply(sel, function(k) { v <- block_fields[[k]][P_CANONICAL]
            if (!is.na(v) && v=="YES") 0L else 1L }, integer(1))
          mane  <- vapply(sel, function(k) { v <- block_fields[[k]][P_MANESEL]
            if (!is.na(v) && v!="") 0L else 1L }, integer(1))
          feat  <- vapply(sel, function(k) { v <- block_fields[[k]][P_Feature]
            if (is.na(v)) "zzz" else v }, character(1))
          ord <- order(ranks, canon, mane, feat)
          chosen <- block_fields[[ sel[ord[1]] ]]
        } else {
          chosen <- rep(NA_character_, n_csq)
        }

        cons_str <- chosen[P_Cons]; top_term <- most_severe_term(cons_str)
        var_class <- if (is.na(top_term)) "Targeted_Region" else vep_to_maf_class(top_term, vt)
        symbol <- chosen[P_SYMBOL]
        hugo <- if (is.na(symbol) || symbol == "") "Unknown" else symbol

        gc <- gt_codes_for_alt(gt, ai)
        map_code <- function(code) {
          if (is.na(code)) return(".")
          if (code == 0L)  return(coords$ref_allele)
          if (code == ai)  return(coords$tum_allele2)
          if (code >= 1L && code <= length(alts)) return(maf_coords(pos, ref, alts[code])$tum_allele2)
          "."
        }
        t_allele1 <- map_code(gc$c1); t_allele2 <- map_code(gc$c2)

        t_ref_count <- if (length(ad_vec) >= 1 && !any(is.na(ad_vec))) ad_vec[1] else NA_character_
        t_alt_count <- if (length(ad_vec) >= (ai+1) && !any(is.na(ad_vec))) ad_vec[ai+1] else NA_character_

        exist <- chosen[P_Existing]; dbsnp <- NA_character_
        if (!is.na(exist) && exist != "") {
          rs <- grep("^rs", strsplit(exist,"&",fixed=TRUE)[[1]], value = TRUE)
          if (length(rs)) dbsnp <- rs[1]
        }
        if (is.na(dbsnp) && !is.na(vid) && vid != ".") {
          rs <- grep("^rs", strsplit(vid,";",fixed=TRUE)[[1]], value = TRUE)
          if (length(rs)) dbsnp <- rs[1]
        }

        hgvsp <- url_decode(chosen[P_HGVSp]); hgvsc <- url_decode(chosen[P_HGVSc])
        vep_vals <- as.list(chosen); names(vep_vals) <- csq_fields

        oi <- oi + 1L
        out[[oi]] <- c(
          list(
            Hugo_Symbol=hugo, Entrez_Gene_Id="0", Center=".", NCBI_Build=ncbi_build,
            Chromosome=chrom, Start_Position=coords$start, End_Position=coords$end,
            Strand="+", Variant_Classification=var_class, Variant_Type=vt,
            Reference_Allele=coords$ref_allele, Tumor_Seq_Allele1=t_allele1,
            Tumor_Seq_Allele2=t_allele2, dbSNP_RS=if (is.na(dbsnp)) "" else dbsnp,
            Tumor_Sample_Barcode=sample_name, Match_Norm_Seq_Allele1="",
            Match_Norm_Seq_Allele2="", HGVSc=if (is.na(hgvsc)) "" else hgvsc,
            HGVSp=if (is.na(hgvsp)) "" else hgvsp, HGVSp_Short=make_hgvsp_short(hgvsp),
            Transcript_ID={ v <- chosen[P_Feature]; if (is.na(v)) "" else v },
            Consequence=if (is.na(cons_str)) "" else cons_str,
            t_depth=if (is.na(sdp)) "" else sdp,
            t_ref_count=if (is.na(t_ref_count)) "" else t_ref_count,
            t_alt_count=if (is.na(t_alt_count)) "" else t_alt_count
          ),
          vep_vals,
          list(
            FILTER=filt,
            QUAL=qual,
            INFO_DP=if (is.na(info_DP)) "" else info_DP,
            INFO_AC=if (length(ac_vec)>=ai && !any(is.na(ac_vec))) ac_vec[ai] else if (!is.na(info_AC)) info_AC else "",
            INFO_AF=if (length(af_vec)>=ai && !any(is.na(af_vec))) af_vec[ai] else if (!is.na(info_AF)) info_AF else "",
            INFO_MQ=if (is.na(info_MQ)) "" else info_MQ,
            INFO_QD=if (is.na(info_QD)) "" else info_QD,
            CNN_1D=if (is.na(info_CNN)) "" else info_CNN,
            GT=gt, AD=if (is.na(ad)) "" else ad,
            sample_DP=if (is.na(sdp)) "" else sdp, GQ=if (is.na(gq)) "" else gq
          )
        )
      }
    }
    res <- data.table::rbindlist(out[seq_len(oi)], use.names = TRUE, fill = TRUE)
    data.table::setattr(res, "n_dropped", n_dropped)   # v4: carry filter drop count to caller
    res
  }

  ## 1j-v2. Convert ONE vcf file end-to-end (header parse + chunked streaming).
  ##        Returns a per-file MAF data.table. Progress printed if verbose.
  convert_one_vcf <- function(path, file_idx, file_total, file_total_records = NA_integer_) {
    csq_fields  <- get_csq_fields(path)
    sample_name <- get_sample_name(path)
    if (verbose)
      message(sprintf("[file %d/%d] Converting %s (sample: %s) | CSQ fields: %d",
                      file_idx, file_total, basename(path), sample_name, length(csq_fields)))

    con <- gzfile(path, "r")
    # Exception-safe cleanup: close on any early return/error so a partially-read
    # connection is never left for the GC to reclaim (which would emit
    # "closing unused connection N (<file>)" warnings). On the NORMAL path the explicit
    # close(con) below already ran, leaving `con` invalid; isOpen() then ERRORS with
    # "invalid connection", so the close attempt is wrapped in tryCatch and any
    # error/warning is swallowed (this on.exit is a safety net, not the primary close).
    on.exit(tryCatch(close(con), error = function(e) NULL, warning = function(w) NULL),
            add = TRUE)
    repeat { l <- readLines(con, 1L); if (length(l)==0L || startsWith(l, "#CHROM")) break }

    chunks <- list(); ci2 <- 0L; total_in <- 0L; n_drop_file <- 0L; t0 <- Sys.time()
    file_done <- 0L   # per-file running record counter (numerator of progress line)
    repeat {
      lines <- readLines(con, chunk_size)
      if (length(lines) == 0L) break
      ci2 <- ci2 + 1L
      mat <- data.table::tstrsplit(lines, "\t", fixed = TRUE)
      dtc <- data.table::data.table(CHROM=mat[[1]], POS=mat[[2]], ID=mat[[3]], REF=mat[[4]],
                        ALT=mat[[5]], QUAL=mat[[6]], FILTER=mat[[7]], INFO=mat[[8]],
                        FORMAT=mat[[9]], SAMPLE=mat[[10]])
      ck <- convert_chunk(dtc, csq_fields, sample_name)
      nd <- attr(ck, "n_dropped"); if (is.null(nd)) nd <- 0L
      n_drop_file <- n_drop_file + nd
      chunks[[ci2]] <- ck
      total_in <- total_in + nrow(dtc)
      file_done <- file_done + nrow(dtc)        # per-file running count (numerator)
      global_done <<- global_done + nrow(dtc)   # update shared global counter (drives global %)
      if (verbose) {
        # Per-file count (file_done / this file's total), GLOBAL percentage
        # (global_done / grand_total), GLOBAL elapsed (since whole-run start, t_global).
        el <- as.numeric(difftime(Sys.time(), t_global, units = "secs"))
        if (!is.na(file_total_records) && !is.na(grand_total) && grand_total > 0L) {
          message(sprintf("    %s/%s (%.0f%%) | %.0fs",
                          format(file_done, big.mark = ","),
                          format(file_total_records, big.mark = ","),
                          100 * global_done / grand_total, el))
        } else {
          message(sprintf("    %s records | %.0fs",
                          format(file_done, big.mark = ","), el))
        }
      }
    }
    close(con)

    maf_one <- data.table::rbindlist(chunks, use.names = TRUE, fill = TRUE)
    if (verbose) {
      # "records" = raw VCF variant records read for this file (pre DP/GQ filter,
      # pre multi-allelic split). The genotype-filter line below reports how many
      # of these were kept/dropped. (Final allele-row count is in the COMBINED /
      # Final Table Dimensions lines.)
      message(sprintf("[file %d/%d] done: %s records",
                      file_idx, file_total, format(total_in, big.mark = ",")))
      if (filter_dp || filter_gq) {
        kept <- total_in - n_drop_file
        crit <- paste(c(if (filter_dp) sprintf("DP>%g", min_DP),
                        if (filter_gq) sprintf("GQ>%g", min_GQ)), collapse = " & ")
        message(sprintf("    genotype filter (%s): kept %s / %s records (dropped %s; %.1f%%).",
                        crit,
                        format(kept, big.mark = ","), format(total_in, big.mark = ","),
                        format(n_drop_file, big.mark = ","),
                        if (total_in > 0L) 100 * n_drop_file / total_in else 0))
      }
    }
    attr(maf_one, "sample_name") <- sample_name
    attr(maf_one, "n_dropped")   <- n_drop_file
    maf_one
  }

  # ==========================================================================
  # 2. Convert every file and ROW-BIND into one combined MAF
  # ==========================================================================
  t_all <- Sys.time()

  ## 2a. Progress display in COUNT-ONLY mode (no pre-count pass).
  ##     PERF (v5): the previous implementation fully decompressed every gzip ONCE
  ##     just to count post-#CHROM variant lines (to drive a global % bar), then
  ##     convert_one_vcf decompressed each file AGAIN to convert it -> every file
  ##     was read twice. On the 2-file test set the redundant pre-count pass cost
  ##     ~34.5 s (~8% of runtime). Since the count is used ONLY for the cosmetic
  ##     progress percentage (it never touches the MAF), we drop it and let the
  ##     conversion report a cumulative running record count + elapsed time
  ##     instead of a %. This is byte-identical in OUTPUT; only verbose progress
  ##     text changes (and only when verbose=TRUE). convert_one_vcf already
  ##     handles a NA per-file total via its count-only branch.
  per_file_total <- rep(NA_integer_, length(vcf_paths))
  grand_total    <- NA_integer_

  if (verbose)
    message(sprintf("Converting %d file(s).", length(vcf_paths)))

  # global running counter + timer shared with convert_one_vcf (via <<-)
  global_done <- 0L
  t_global    <- t_all

  # v6 (PERF): optional multi-core conversion ACROSS FILES. Each VCF is converted
  # independently (its own MAF chunk) before the final rbindlist, so files are
  # embarrassingly parallel. ncores>1 forks parallel::mclapply over the file list;
  # mclapply preserves input order, so `per_file` (and thus the combined MAF) is
  # byte-identical to the sequential path regardless of ncores. Only effective with
  # >1 file and a fork-capable OS; otherwise we run the original sequential loop.
  # Default ncores=1L => exactly the previous behaviour. The shared global_done<<-
  # progress counter is fork-local (per child), so per-chunk % lines are suppressed
  # in workers to avoid garbled interleaved output; a clean per-file summary prints
  # on collection. parallel is base R (no new dependency).
  n_files  <- length(vcf_paths)
  want_par <- is.numeric(ncores) && length(ncores) == 1L && !is.na(ncores) &&
              ncores > 1L && n_files > 1L &&
              .Platform$OS.type == "unix" &&
              requireNamespace("parallel", quietly = TRUE)
  use_cores <- if (want_par)
                 max(1L, min(as.integer(ncores), parallel::detectCores(), n_files))
               else 1L

  if (want_par && use_cores > 1L) {
    if (verbose)
      message(sprintf("Parallel conversion: %d worker(s) over %d file(s) (mclapply).",
                      use_cores, n_files))
    per_file <- parallel::mclapply(
      seq_along(vcf_paths),
      function(i) suppressMessages(
        convert_one_vcf(vcf_paths[i], i, n_files, per_file_total[i])),
      mc.cores = use_cores, mc.preschedule = FALSE)
    # surface any worker error as a hard error (mclapply returns try-error objects)
    errs <- vapply(per_file, function(x) inherits(x, "try-error"), logical(1))
    if (any(errs)) {
      msgs <- vapply(per_file[errs], function(e) {
        cond <- attr(e, "condition")
        if (!is.null(cond)) conditionMessage(cond) else as.character(e)
      }, character(1))
      stop(sprintf("read.gvr: %d file(s) failed during parallel conversion:\n%s",
                   sum(errs), paste(unique(msgs), collapse = "\n")))
    }
    sample_names <- vapply(per_file, function(x) attr(x, "sample_name"), character(1))
    if (verbose) for (i in seq_along(vcf_paths))
      message(sprintf("[file %d/%d] done: %s (sample: %s, %s MAF rows)",
                      i, n_files, basename(vcf_paths[i]), sample_names[i],
                      format(nrow(per_file[[i]]), big.mark = ",")))
  } else {
    per_file <- vector("list", n_files)
    sample_names <- character(n_files)
    for (i in seq_along(vcf_paths)) {
      per_file[[i]] <- convert_one_vcf(vcf_paths[i], i, n_files,
                                       per_file_total[i])
      sample_names[i] <- attr(per_file[[i]], "sample_name")
    }
  }
  # collision guard: warn if two files map to the same sample barcode
  dup_samp <- sample_names[duplicated(sample_names)]
  if (length(dup_samp) > 0L)
    warning(sprintf("Duplicate Tumor_Sample_Barcode(s) across files: %s",
                    paste(unique(dup_samp), collapse=", ")))

  maf <- data.table::rbindlist(per_file, use.names = TRUE, fill = TRUE)
  if (verbose) message(sprintf("COMBINED: %d file(s) | %d MAF rows (%d cols) in %.0fs",
                               length(vcf_paths), nrow(maf), ncol(maf),
                               as.numeric(difftime(Sys.time(), t_all, units="secs"))))

  # ==========================================================================
  # 3. POST-PROCESSING (v2)
  # ==========================================================================

  ## 3a. Remove VERIFIED-IDENTICAL duplicate-named columns -------------------
  ##     For each repeated name: compare columns elementwise (NA==NA). If all
  ##     rows identical, keep the FIRST occurrence and drop the rest. If any
  ##     row differs, keep both and rename later copies "<name>.csq" + warn.
  ##     NOTE: run BEFORE HGVS stripping, so identity is checked on raw values
  ##     (otherwise stripping the core copy would make it differ from the CSQ
  ##     copy and the strict check would refuse to drop it).
  if (dedup_columns) {
    cn <- names(maf)
    dup_names <- unique(cn[duplicated(cn)])
    dropped <- character(0); renamed <- character(0)
    # Vectorized URL-decode (fast; mirrors the per-element url_decode helper).
    url_decode_vec <- function(x) {
      x <- gsub("%3D","=",x,fixed=TRUE); x <- gsub("%3B",";",x,fixed=TRUE)
      x <- gsub("%2C",",",x,fixed=TRUE); x <- gsub("%3A",":",x,fixed=TRUE); x
    }
    for (nm in dup_names) {
      idx <- which(cn == nm)
      keep <- idx[1]
      # Columns are "identical" iff, after URL-decoding and treating ""/NA as the
      # SAME 'missing' token, every row matches. This is the correct equivalence:
      # the MAF-core copy writes "" for missing while the raw CSQ copy writes NA,
      # so a naive (a==b) leaves NA at those rows and must not count as a difference.
      identical_to_keep <- function(j) {
        a <- url_decode_vec(maf[[keep]]); b <- url_decode_vec(maf[[j]])
        a[is.na(a)] <- ""        # normalize NA -> "" (missing)
        b[is.na(b)] <- ""
        all(a == b)              # no NA left, so all() is unambiguous
      }
      to_drop <- integer(0)
      for (j in idx[-1]) {
        if (identical_to_keep(j)) to_drop <- c(to_drop, j)
        else renamed <- c(renamed, nm)
      }
      if (length(to_drop) > 0L) {
        maf[, (to_drop) := NULL]
        dropped <- c(dropped, nm)
        cn <- names(maf)  # refresh after deletion
      }
    }
    # rename any surviving same-named-but-differing duplicates to avoid collision
    cn <- names(maf)
    if (any(duplicated(cn))) {
      for (nm in unique(cn[duplicated(cn)])) {
        idx <- which(cn == nm)
        for (k in seq_along(idx)[-1]) {
          new <- paste0(nm, ".csq", if (k > 2) k - 1 else "")
          data.table::setnames(maf, idx[k], new)
        }
        cn <- names(maf)
      }
    }
  }

  ## 3b. Strip Ensembl feature prefix from HGVSc / HGVSp ---------------------
  ##     "ENST...:n.1889G>A" -> "n.1889G>A". Applied to the surviving core cols.
  if (strip_hgvs_prefix) {
    # Vectorized strip: remove up to first ':'. Empty/NA pass through unchanged.
    strip_prefix_vec <- function(x) ifelse(is.na(x) | x == "", x, sub("^[^:]*:", "", x))
    for (col in c("HGVSc", "HGVSp")) {
      if (col %in% names(maf)) {
        data.table::set(maf, j = col, value = strip_prefix_vec(maf[[col]]))   # in-place, no copy
      }
    }
  }

  ## 3c. Add Genotype = Tumor_Seq_Allele1 / Tumor_Seq_Allele2 ----------------
  if (add_genotype) {
    data.table::setDT(maf)   # ensure a clean data.table (prior := NULL may have left a shallow copy)
    if (all(c("Tumor_Seq_Allele1","Tumor_Seq_Allele2") %in% names(maf))) {
      maf[, Genotype := paste(Tumor_Seq_Allele1, Tumor_Seq_Allele2, sep = "/")]
      # place Genotype directly after Tumor_Seq_Allele2
      nm <- names(maf)
      after <- which(nm == "Tumor_Seq_Allele2")
      neworder <- append(setdiff(nm, "Genotype"), "Genotype", after = after)
      data.table::setcolorder(maf, neworder)
    }
  }

  ## 3d. ABraOM Brazilian allele frequency (SABE-609-WES) --------------------
  ##     Adds ONE column `ABraOM_AF` = the ABraOM `Frequencies` value.
  ##
  ##     CRITICAL: the ABraOM 609-exome file is hg19/GRCh37 while these VCFs are
  ##     GRCh38. Positions are therefore NOT comparable (same rsID can sit >60 kb
  ##     apart between builds). We join on dbSNP rsID + alleles, which ARE
  ##     build-stable: MAF.dbSNP_RS == ABraOM.avsnp147 AND Ref/Alt agree.
  ##     (rsID alone is ambiguous at multi-allelic sites; alleles disambiguate.)
  ##     Blank "" where there is no rsID+allele match. Frequency only.
  if (add_abraom) {
    data.table::setDT(maf)
    maf[, ABraOM_AF := ""]                  # default blank (no-match convention)
    ok <- TRUE
    # Resolve the ABraOM file: explicit path, else cached copy, else download.
    apath <- abraom_path
    if (is.null(apath)) {
      cache_dir <- if (!is.null(cache_dir)) cache_dir
                   else tools::R_user_dir("germlinevaR", which = "cache")
      apath <- file.path(cache_dir, "ABRaOM_60plus_SABE_609_exomes_annotated.gz")
      # Known good size for the ABraOM reference (~50 MB gzipped). Used to detect
      # truncated partial downloads (e.g. from a timeout) that would crash fread().
      abraom_expected_bytes <- 50242984L
      need_dl <- !file.exists(apath) ||
                 file.info(apath)$size < abraom_expected_bytes * 0.9   # tolerate minor size drift
      if (need_dl) {
        # Delete any truncated partial file from a previous failed download
        if (file.exists(apath)) {
          if (verbose) message("ABraOM: cached file appears truncated; re-downloading ...")
          unlink(apath)
        }
        ok <- tryCatch({
          if (!dir.exists(cache_dir)) {
            dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
            # dir.create can silently fail on some mounts; verify and retry via shell
            if (!dir.exists(cache_dir)) {
              system2("mkdir", c("-p", shQuote(cache_dir)))
              if (!dir.exists(cache_dir))
                stop("cannot create cache directory: ", cache_dir,
                     ". Set cache_dir= to a writable path, or download the ABraOM file ",
                     "manually and pass it via abraom_path=.")
            }
          }
          if (verbose) message("ABraOM: downloading reference file (~50 MB, one-time cache) ...")
          utils::download.file(abraom_url, apath, mode = "wb", quiet = !verbose,
                               timeout = 600L)
          file.exists(apath) && file.info(apath)$size >= abraom_expected_bytes * 0.9
        }, error = function(e) {
          # Clean up partial file so the next attempt retries from scratch
          if (file.exists(apath)) unlink(apath)
          warning("read.gvr: ABraOM reference could not be downloaded (",
                  conditionMessage(e), "); 'ABraOM_AF' left blank. ",
                  "You can download it manually from:\n  ", abraom_url,
                  "\nand pass the local path via abraom_path=.", call. = FALSE)
          FALSE
        })
      }
    } else if (!file.exists(apath)) {
      warning("read.gvr: ABraOM reference path not found: ", apath,
              "; 'ABraOM_AF' left blank.", call. = FALSE)
      ok <- FALSE
    }

    if (ok) {
      ab <- tryCatch(suppressWarnings(data.table::fread(apath, sep = "\t", header = TRUE, quote = "",
                           showProgress = FALSE)), error = function(e) NULL)
      if (is.null(ab)) {
        warning("read.gvr: ABraOM reference unreadable; 'ABraOM_AF' left blank.",
                call. = FALSE)
      } else {
        data.table::setnames(ab, make.names(names(ab)))
        need <- c("avsnp147","Ref","Alt","Frequencies")
        if (!all(need %in% names(ab))) {
          warning("read.gvr: ABraOM file lacks expected columns ",
                  "(avsnp147/Ref/Alt/Frequencies); 'ABraOM_AF' left blank.",
                  call. = FALSE)
        } else {
          # Build keyed lookup: distinct (rsID, REF, ALT) -> Frequencies (as char).
          lut <- unique(ab[avsnp147 != "NA" & !is.na(avsnp147) & avsnp147 != "",
                           .(rs = avsnp147,
                             ref = toupper(Ref), alt = toupper(Alt),
                             af  = as.character(Frequencies))])
          data.table::setkey(lut, rs, ref, alt)
          # MAF join keys (uppercased; both sides use '-' for indel empty allele).
          maf[, `:=`(.rs  = dbSNP_RS,
                     .ref = toupper(Reference_Allele),
                     .alt = toupper(Tumor_Seq_Allele2))]
          hit <- lut[maf[, .(.rs, .ref, .alt)],
                     on = c(rs = ".rs", ref = ".ref", alt = ".alt"),
                     x.af]            # vector of matched AF, NA where no match
          maf[, ABraOM_AF := ifelse(is.na(hit), "", hit)]
          maf[, c(".rs",".ref",".alt") := NULL]
          # Place the new column with the other population-frequency columns:
          # right after MAX_AF_POPS / MAX_AF / the last gnomAD column (whichever
          # exists). This groups it with gnomAD/1000G AFs, NOT the INFO_* block.
          nm <- names(maf)
          anchors <- c(grep("^MAX_AF_POPS$", nm), grep("^MAX_AF$", nm),
                       grep("^gnomAD", nm))
          after <- if (length(anchors)) max(anchors) else length(nm) - 1L
          neworder <- append(setdiff(nm, "ABraOM_AF"),
                              "ABraOM_AF", after = after)
          data.table::setcolorder(maf, neworder)
          if (verbose) message(sprintf(
            "ABraOM: annotated %d/%d rows (%.1f%%) with ABraOM frequency (rsID+allele join).",
            sum(nzchar(maf$ABraOM_AF)), nrow(maf),
            100*mean(nzchar(maf$ABraOM_AF))))
        }
      }
    }
    # NOTE: when the reference is unavailable (download/path/read/column failure),
    # the specific warning is emitted at the point of failure above; `ABraOM_AF`
    # was initialized blank, so it is left blank and conversion continues.
  }

  ## 3d-bis. Gene-of-interest subset (opt-in) -------------------------------
  ##     Keep only rows whose Hugo_Symbol matches the user-supplied `genes`
  ##     vector (exact, case-insensitive). NULL = keep all genes. Applied as a
  ##     FINAL row filter, AFTER all annotation but BEFORE drop_empty_cols, so
  ##     emptiness is evaluated on the subsetted rows.
  if (!is.null(genes)) {
    data.table::setDT(maf)
    genes_chr <- as.character(genes)
    genes_chr <- genes_chr[!is.na(genes_chr) & nzchar(genes_chr)]
    want <- unique(toupper(trimws(genes_chr)))
    n_before <- nrow(maf)
    have <- toupper(trimws(as.character(maf$Hugo_Symbol)))
    maf <- maf[have %in% want]
    if (verbose) {
      found    <- intersect(want, unique(have))
      notfound <- setdiff(want, unique(have))
      message(sprintf(
        "gene subset: kept %s / %s rows across %d / %d requested gene(s)%s.",
        format(nrow(maf), big.mark = ","), format(n_before, big.mark = ","),
        length(found), length(want),
        if (length(notfound)) paste0(" (not found: ", paste(notfound, collapse = ", "), ")") else ""))
    }
  }

  ## 3d-ter. Variant-Classification non-synonymous filter (opt-in) -----------
  ##     Mirrors maftools' vc_nonSyn: keep only protein-altering classes.
  ##     FALSE (default) = keep all. TRUE = keep the 9 standard High/Moderate
  ##     classes. A custom character vector = keep those specific classes.
  ##     Applied AFTER gene subset, BEFORE drop_empty_cols.
  if (!identical(vc_nonSyn, FALSE)) {
    data.table::setDT(maf)
    vc_default <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site",
                    "Translation_Start_Site", "Nonsense_Mutation", "Nonstop_Mutation",
                    "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation")
    vc_keep <- if (isTRUE(vc_nonSyn)) vc_default else as.character(vc_nonSyn)
    vc_keep <- vc_keep[!is.na(vc_keep) & nzchar(vc_keep)]
    if (length(vc_keep) > 0L) {
      n_before <- nrow(maf)
      vc_col <- as.character(maf$Variant_Classification)
      # Remove rows with missing/blank Variant_Classification when filter is active
      maf <- maf[!is.na(vc_col) & nzchar(vc_col) & vc_col %in% vc_keep]
      if (verbose) {
        dropped <- n_before - nrow(maf)
        message(sprintf(
          "vc_nonSyn: kept %s / %s rows (%d classification(s): %s; removed %s silent/other)",
          format(nrow(maf), big.mark = ","), format(n_before, big.mark = ","),
          length(vc_keep), paste(vc_keep, collapse = ", "),
          format(dropped, big.mark = ",")))
      }
    }
  }

  ## 3e. Drop all-empty columns (opt-in) ------------------------------------
  ##     Remove any column whose values are ALL missing (NA or "") across the
  ##     entire combined table. Default FALSE = keep full schema.
  if (drop_empty_cols) {
    data.table::setDT(maf)
    is_all_empty <- function(col) {
      v <- maf[[col]]; vc <- as.character(v)
      all(is.na(v) | vc == "")
    }
    empty_cols <- names(maf)[vapply(names(maf), is_all_empty, logical(1))]
    if (length(empty_cols) > 0L) {
      maf[, (empty_cols) := NULL]
      if (verbose) message(sprintf("drop_empty_cols: removed %d all-empty column(s): %s",
                                   length(empty_cols), paste(empty_cols, collapse=", ")))
    } else if (verbose) {
      message("drop_empty_cols: no all-empty columns found.")
    }
  }

  if (verbose) message(sprintf("Final Table Dimensions: %d rows x %d columns.", nrow(maf), ncol(maf)))

  # Tag annotator BEFORE writes so saved TSV/RDS/XLSX carry the attribute too.
  # (R's saveRDS preserves attributes; for TSV/XLSX it's only for in-memory
  #  consumers -- but we still keep it consistent.)
  data.table::setattr(maf, "annotator", "vep")

  # ==========================================================================
  # 4. Optional file outputs
  # ==========================================================================
  if (write_tsv || write_rds || write_xlsx) {
    if (is.null(out_dir)) out_dir <- if (is.null(vcf_path)) folder else dirname(vcf_path)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (is.null(out_prefix)) {
      if (length(vcf_paths) > 1L) {
        out_prefix <- "combined.maf"                                  # multi-file
      } else {
        out_prefix <- sub("\\.vcf\\.gz$", "", basename(vcf_paths[1])) # single file
        out_prefix <- paste0(out_prefix, ".maf")
      }
    }
    if (write_tsv) {
      tsv_path <- file.path(out_dir, paste0(out_prefix, ".tsv"))
      data.table::fwrite(maf, tsv_path, sep = "\t", quote = FALSE, na = "")
      if (verbose) message("Wrote TSV: ", tsv_path)
    }
    if (write_rds) {
      # IMPORTANT: on S3-backed mounts (e.g. /mnt/results), R's file.copy() can
      # silently produce a 0-byte file. The ONLY reliable route is to saveRDS to
      # a local POSIX path (tempdir / /workspace) and then shell `cp` it over.
      rds_final <- file.path(out_dir, paste0(out_prefix, ".rds"))
      tmp_rds   <- file.path(tempdir(), paste0(out_prefix, ".rds"))
      saveRDS(maf, tmp_rds, compress = TRUE)
      system2("cp", c(shQuote(tmp_rds), shQuote(rds_final)))   # always shell-cp
      sz <- suppressWarnings(file.info(rds_final)$size)
      if (is.na(sz) || sz == 0)
        warning(sprintf("RDS write may have failed (0 bytes): %s", rds_final))
      else if (verbose)
        message(sprintf("Wrote RDS: %s (%.0f MB)", rds_final, sz/1e6))
    }
    if (write_xlsx) {
      # B1: optional Excel export of the FINAL MAF (one "MAF" sheet). Mirrors the
      # FUSE-safe openxlsx pattern used by gvr_summary: build the workbook, save to a
      # local temp file, then shell-cp to out_dir (openxlsx uses zip random-access
      # writes that can silently 0-byte on S3-backed mounts). Degrades gracefully:
      # if openxlsx is absent we warn and skip (TSV/RDS, if requested, still wrote).
      # NOTE: a germline MAF can be very large (hundreds of thousands of rows). Excel
      # handles it but the file is big and slow to open; write_rds / write_tsv remain
      # the better choice for large tables / downstream R use.
      if (!requireNamespace("openxlsx", quietly = TRUE)) {
        warning("read.gvr: 'openxlsx' not installed; skipping Excel export.")
      } else {
        xlsx_final <- file.path(out_dir, paste0(out_prefix, ".xlsx"))
        if (file.exists(xlsx_final) && verbose)
          message(sprintf("  Overwriting existing Excel: %s", xlsx_final))
        if (nrow(maf) > 1000000L)
          warning(sprintf("read.gvr: MAF has %s rows; Excel's per-sheet limit is 1,048,576 rows.",
                          format(nrow(maf), big.mark = ",")))
        wb <- openxlsx::createWorkbook()
        hs <- openxlsx::createStyle(textDecoration = "bold", halign = "center")
        openxlsx::addWorksheet(wb, "MAF")
        openxlsx::writeData(wb, "MAF", as.data.frame(maf), headerStyle = hs)
        openxlsx::freezePane(wb, "MAF", firstRow = TRUE)
        openxlsx::setColWidths(wb, "MAF", cols = seq_len(ncol(maf)), widths = "auto")
        tmp_xlsx <- file.path(tempdir(), paste0(out_prefix, ".xlsx"))
        wrote_ok <- tryCatch({ openxlsx::saveWorkbook(wb, tmp_xlsx, overwrite = TRUE); TRUE },
                             error = function(e) {
                               warning(sprintf("read.gvr: Excel write failed: %s", conditionMessage(e))); FALSE })
        if (wrote_ok) {
          system2("cp", c(shQuote(tmp_xlsx), shQuote(xlsx_final)))
          sz <- suppressWarnings(file.info(xlsx_final)$size)
          if (is.na(sz) || sz == 0) {
            warning(sprintf("read.gvr: copy to '%s' may have failed; Excel left at '%s'.",
                            xlsx_final, tmp_xlsx))
            xlsx_final <- tmp_xlsx
          }
          if (verbose) message(sprintf("Wrote Excel: %s", xlsx_final))
        }
      }
    }
  }

  # set a friendly key-less data.table and return
  data.table::setDT(maf)
  maf[]
}

# If sourced interactively this just defines read.gvr().
# Example (commented):
#   maf <- read.gvr("/mnt/user-uploads")                # merge all *_NN.vcf.gz
#   maf_pass <- maf[FILTER == "PASS"]
#   # genes of interest only, no DP/GQ filter:
#   maf_goi <- read.gvr("/mnt/user-uploads", min_DP = NULL, min_GQ = NULL,
#                       genes = c("MEN1","RET","CDKN1B"))

# NOTE: globalVariables() declarations for this package are consolidated in
# R/globals.R (one package-scoped block covering all functions).
