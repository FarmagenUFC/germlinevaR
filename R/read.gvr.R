#' Convert VEP-annotated germline VCF(s) to an MAF-like data.table
#'
#' @description
#' Converts VEP-annotated, single-sample germline VCFs (GATK HaplotypeCaller ->
#' CNN tranches -> Ensembl VEP, hg38) into an MAF-like table and returns it as
#' an in-memory `data.table` for downstream filtering ([gvr_filter()]) and
#' summarisation ([gvr_summary()]). In folder mode it finds every per-sample VCF,
#' converts each, and row-binds them into one combined gvr table. The conversion
#' uses base R and \pkg{data.table} only - no external annotation-package
#' dependency. This is the recommended entry point for all germline VCFs:
#' [read.gvr()] inspects each input's INFO tags and, when needed, delegates to
#' [read.gvr.snpeff()] for SnpEff-annotated VCFs or [read.gvr.dual()] for VCFs
#' carrying both VEP and SnpEff annotations, so a single call handles every
#' annotator combination.
#'
#' @details
#' Output and behaviour:
#' \itemize{
#'   \item Returns the final MAF-like `data.table`, one row per variant ALLELE
#'     (multi-allelic sites are split).
#'   \item A single most-severe transcript is chosen per allele (VEP severity ->
#'     `CANONICAL` -> `MANE_SELECT` -> transcript id).
#'   \item Columns include the MAF-like core fields, ALL VEP CSQ fields (read from the VCF
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
#'   Also used as the search root for `file=`.
#' @param vcf_path Character vector of one or more full paths to `.vcf.gz`
#'   files to convert. Use this to process a specific set of files outside
#'   the folder pattern. Mutually exclusive with `file=`. `NULL` (default)
#'   selects folder mode.
#' @param file Character vector of basenames (e.g.
#'   `c("S1.vep.vcf.gz", "S2.vep.vcf.gz")`) resolved against `folder=`.
#'   Use this to pick specific files from a folder that contains files you
#'   do NOT want to merge. Mutually exclusive with `vcf_path=`. `NULL`
#'   (default) selects either `vcf_path=` mode or folder-pattern mode.
#' @param pattern Regular expression identifying per-sample VCFs in folder mode.
#'   Default `"\\.vcf\\.gz$"` (matches any `.vcf.gz` file). The old default
#'   `"_\\d+(\\.(vep|snp[eE]ff))?\\.vcf\\.gz$"` (requires `_NN` suffix) is
#'   still available by passing it explicitly.
#' @param write_tsv Logical; if `TRUE`, also write the table as a TSV to `out_dir`.
#'   Default `FALSE`.
#' @param write_rds Logical; if `TRUE`, also write the table as an `.rds` to `out_dir`.
#'   Default `FALSE`.
#' @param write_xlsx Logical; if `TRUE`, also write the table as an `.xlsx` workbook
#'   (single `"gvr_table"` sheet) to `out_dir`. Requires the \pkg{openxlsx} package (a
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
#' @param ncbi_build Reference build label written into the table `NCBI_Build`
#'   column. Default `"auto"`: the function inspects the input VCF header
#'   (VEP `assembly=` / SnpEff `SnpEffCmd` database token / first `##contig=`
#'   length) and picks the canonical label among `"GRCh38"`, `"GRCh37"`,
#'   `"T2T-CHM13v2.0"`. When detection cannot decide, falls back to `"GRCh38"`
#'   with a verbose-mode message. Pass any literal (`"GRCh37"`, `"hg19"`,
#'   `"T2T-CHM13v2.0"`, internal lab codes) to override; the user-supplied
#'   value is written verbatim into `NCBI_Build` and the rest of the pipeline
#'   does not branch on its value. The ABraOM join uses dbSNP rsID + alleles
#'   and is build-stable. The aliases `"hg19"` and `"hg38"` are mapped to
#'   `"GRCh37"` / `"GRCh38"` for the mismatch check; a diagnostic warning
#'   fires when an explicit `ncbi_build` value (canonical or alias) disagrees
#'   with what auto-detection found at high confidence. Off-table user labels
#'   (e.g. internal lab codes) pass through silently.
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
#' @param panel Character vector of curated disease panel name(s) (e.g.
#'   `"breast cancer"`, `"hereditary prostate cancer"`, `"gist"`).
#'   Each name is resolved to a gene vector via [gvr_panel_genes()] and the
#'   union of all resolved genes is taken with `genes` (deduplicated,
#'   uppercased). Names are matched case-insensitively, trimmed, and `_` is
#'   treated as a space, so `"Breast_Cancer"`, `"breast cancer"`, and
#'   `" BREAST CANCER "` all resolve identically. A small alias table is
#'   also recognised (e.g. `"gastrointestinal stromal tumor"` -> `"gist"`,
#'   the typo `"pheocromocytoma"` -> `"pheochromocytoma"`). The registry
#'   currently ships 15 panels; see [gvr_list_panels()] for the full list.
#'   An unknown name raises an error listing the available panels. `NULL`
#'   (default) disables panel filtering; behaviour is then byte-identical to
#'   omitting the argument.
#' @param vc_nonSyn Logical or character vector. Controls which
#'   `Variant_Classification` values are retained (mirroring the
#'   convention of the `vc_nonSyn` argument). `FALSE` (default) keeps ALL variant classifications.
#'   `TRUE` keeps only protein-altering classes (High/Moderate VEP consequences):
#'   `"Frame_Shift_Del"`, `"Frame_Shift_Ins"`, `"Splice_Site"`,
#'   `"Translation_Start_Site"`, `"Nonsense_Mutation"`, `"Nonstop_Mutation"`,
#'   `"In_Frame_Del"`, `"In_Frame_Ins"`, `"Missense_Mutation"`. Alternatively,
#'   pass a custom character vector of classifications to keep. Rows with
#'   missing/blank `Variant_Classification` are always removed when this filter
#'   is active.
#' @param canonical_only Logical; when `TRUE` (default), drops table rows whose
#'   chosen VEP CSQ block has `CANONICAL != "YES"`. read.gvr() already prefers
#'   `CANONICAL=YES` when ranking CSQ blocks; this filter discards the
#'   fallback rows emitted when no canonical block exists for a given ALT.
#'   For SnpEff-annotated input, the ANN field has no CANONICAL flag —
#'   `canonical_only=TRUE` is ignored with a warning, and the result is the
#'   same as `canonical_only=FALSE`.
#' @param ncores Integer; number of worker processes for converting MULTIPLE input
#'   files in parallel via [parallel::mclapply()] (fork-based; Unix/macOS only).
#'   Default `1L` runs sequentially and is byte-identical to previous behaviour.
#'   Values `> 1` only help when more than one VCF is being read (each file is an
#'   independent task) and are clamped to `min(ncores, detectCores(), n_files)`. On
#'   non-fork platforms it falls back to sequential. A single file is unaffected.
#' @param verbose Logical; if `TRUE` (default) print per-file and per-chunk progress
#'   (file i/N, cumulative records, elapsed seconds).
#' @return An MAF-like `data.table`: one row per variant allele, with MAF-like core columns, all
#'   VEP CSQ fields, key GATK QC fields, `Tumor_Sample_Barcode`, and (when enabled) the
#'   `Genotype` and `ABraOM_AF` columns. TSV/RDS files are written as a side
#'   effect when `write_tsv`/`write_rds` is `TRUE`.
#'
#' @seealso [gvr_filter()] to filter the returned table, [gvr_summary()] to summarise it,
#'   [read.gvr.snpeff()] for SnpEff-annotated VCFs, [read.gvr.dual()] for VCFs with both
#'   VEP and SnpEff annotations.
#' @family germlinevaR
#' @author germlinevaR authors
#'
#' @examples
#' ## read.gvr() reads VEP/SnpEff-annotated VCFs and returns a parsed
#' ## data.table. Reading the bundled 62-variant fixture takes ~20s on
#' ## a typical CI worker, so the real read.gvr() calls are wrapped in
#' ## a donttest block below; a pre-parsed equivalent .rds is also bundled.
#'
#' ## The function is exported and callable:
#' is.function(read.gvr)
#'
#' ## Pre-parsed equivalent (instantaneous; same shape and content as the
#' ## result of read.gvr() on the same fixture, minus the ABraOM_AF column
#' ## which was added in a later code revision):
#' gvr <- readRDS(system.file("extdata", "example_gvr.rds",
#'     package = "germlinevaR"))
#' dim(gvr)
#'
#' \donttest{
#' ## Real read.gvr() call on the bundled VCF directory:
#' vcf_dir <- system.file("extdata", package = "germlinevaR")
#' gvr <- read.gvr(vcf_dir, verbose = FALSE)
#' dim(gvr)            # 62 rows x 116 columns
#'
#' ## Same call but with write-out to tempdir
#' out <- tempdir()
#' gvr2 <- read.gvr(vcf_dir, write_tsv = TRUE, write_rds = TRUE,
#'     out_dir = out, verbose = FALSE)
#' list.files(out, pattern = "\\.(tsv|rds)$")
#'
#' ## Single-file mode: full path to the bundled VCF
#' vcf_file <- system.file("extdata", "example.vep.vcf.gz",
#'     package = "germlinevaR")
#' gvr3 <- read.gvr(vcf_path = vcf_file, verbose = FALSE)
#' nrow(gvr3)
#' }
#' @importFrom data.table data.table as.data.table rbindlist fread fwrite setnames setcolorder setDT set setattr tstrsplit setkey :=
#' @importFrom stats setNames
#' @importFrom utils download.file
#' @importFrom openxlsx createWorkbook
#' @export
read.gvr <- function(folder = ".",
                     vcf_path   = NULL,
                     file       = NULL,   # vN+4: basename(s) inside `folder`; mutex with vcf_path
                     pattern    = "\\.vcf\\.gz$",   # vN+4: any .vcf.gz; old default was "_\\d+(\\.(vep|snp[eE]ff))?\\.vcf\\.gz$"
                     write_tsv  = FALSE,
                     write_rds  = FALSE,
                     write_xlsx = FALSE,   # v6: also write the table as .xlsx (one "gvr_table" sheet)
                     out_dir    = NULL,
                     out_prefix = NULL,
                     chunk_size = 25000L,
                     ncbi_build = "auto",
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
                     panel             = NULL,   # vN: curated disease gene panel(s); union'd with `genes`
                     vc_nonSyn         = FALSE,  # v8: keep only protein-altering Variant_Classification
                     canonical_only    = TRUE,   # vN+4: drop rows whose chosen CSQ block has CANONICAL != "YES"
                     ncores            = 1L,     # v6: parallel files (>1 forks mclapply; 1 = sequential, default)
                     verbose    = TRUE) {
    # ===========================================================================

    # Auto-source siblings if not already available
    # (in a package, they always are; in standalone mode, this provides them)
    if (!exists("read.gvr.snpeff", mode = "function", inherits = TRUE)) {
        sib <- .gvr_locate_sibling("read.gvr.snpeff.R")
        if (!is.null(sib)) {
            tryCatch(
                source(sib, local = FALSE, chdir = FALSE),
                error = function(e) warning(sprintf(
                    "read.gvr: failed to source sibling '%s': %s", sib, conditionMessage(e)))
            )
        }
    }
    if (!exists("read.gvr.dual", mode = "function", inherits = TRUE)) {
        sib <- .gvr_locate_sibling("read.gvr.dual.R")
        if (!is.null(sib)) {
            tryCatch(
                source(sib, local = FALSE, chdir = FALSE),
                error = function(e) warning(sprintf(
                    "read.gvr: failed to source sibling '%s': %s", sib, conditionMessage(e)))
            )
        }
    }

    # ==========================================================================
    # vN setup: disease-panel resolution (Phase N)
    #   If `panel` is supplied, resolve each name to its gene vector via
    #   .gvr_resolve_panels() and UNION with `genes`. The resulting vector is
    #   stored back in `genes` so the downstream filtering machinery (rough
    #   INFO prefilter + final exact Hugo_Symbol filter) needs no further
    #   changes. When `panel = NULL` this block is a strict no-op and `genes`
    #   is byte-identical to the caller-supplied value.
    # ==========================================================================
    if (!is.null(panel)) {
        if (!exists(".gvr_resolve_panels", mode = "function", inherits = TRUE)) {
            stop("read.gvr: panel resolution requires gvr_panels.R; please source/load the package.")
        }
        panel_genes <- .gvr_resolve_panels(panel)        # uppercased, deduplicated
        if (length(panel_genes) > 0L) {
            effective_genes <- unique(toupper(trimws(c(panel_genes, as.character(genes)))))
            effective_genes <- effective_genes[!is.na(effective_genes) & nzchar(effective_genes)]
            # Phase N+1: capture pre-union panel state for the post-filter
            # "panel subset:" coverage message (sibling of "gene subset:" below).
            # Only created when panel was supplied AND resolved to >=1 gene; the
            # post-filter block gates on exists(.panel_genes_for_summary).
            .panel_genes_for_summary  <- panel_genes
            .panel_extras_for_summary <- setdiff(
                toupper(trimws(as.character(genes))), panel_genes
            )
            .panel_extras_for_summary <- .panel_extras_for_summary[
                !is.na(.panel_extras_for_summary) & nzchar(.panel_extras_for_summary)
            ]
            if (verbose) {
                message(sprintf(
                    "panel: resolved %d panel(s) -> %d unique Hugo_Symbol(s)%s.",
                    length(unique(as.character(panel))),
                    length(effective_genes),
                    if (!is.null(genes) && length(as.character(genes)) > 0L)
                        sprintf(" (panel %d + extras %d)",
                            length(panel_genes),
                            length(setdiff(toupper(trimws(as.character(genes))), panel_genes)))
                    else ""))
            }
            genes <- effective_genes
        }
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
    #    Three input modes, precedence vcf_path > file > folder:
    #      - vcf_path : character vector of full paths to process (multi-file OK).
    #      - file     : character vector of basenames resolved against `folder`.
    #      - folder   : list.files(folder, pattern) -- the default fan-out mode.
    #    `vcf_path` and `file` are mutually exclusive.
    # ==========================================================================
    if (!is.null(vcf_path) && !is.null(file))
        stop("read.gvr: pass either `vcf_path=` (full paths) or `file=` ",
            "(basenames inside `folder=`), not both.", call. = FALSE)

    if (!is.null(vcf_path)) {
        vcf_path <- as.character(vcf_path)
        if (!length(vcf_path) || any(!nzchar(vcf_path)))
            stop("read.gvr: `vcf_path=` must be a non-empty character vector.",
                call. = FALSE)
        miss <- vcf_path[!file.exists(vcf_path)]
        if (length(miss))
            stop(sprintf("read.gvr: vcf_path file(s) do not exist:\n%s",
                paste0("  - ", miss, collapse = "\n")), call. = FALSE)
        vcf_paths <- vcf_path
    } else if (!is.null(file)) {
        file <- as.character(file)
        if (!length(file) || any(!nzchar(file)))
            stop("read.gvr: `file=` must be a non-empty character vector of basenames.",
                call. = FALSE)
        if (!dir.exists(folder))
            stop(sprintf("Folder does not exist: %s", folder), call. = FALSE)
        candidates <- file.path(folder, file)
        miss <- file[!file.exists(candidates)]
        if (length(miss))
            stop(sprintf("read.gvr: file(s) not found under folder='%s':\n%s",
                folder, paste0("  - ", miss, collapse = "\n")),
            call. = FALSE)
        vcf_paths <- candidates
    } else {
        if (!dir.exists(folder))
            stop(sprintf("Folder does not exist: %s", folder), call. = FALSE)
        hits <- list.files(folder, pattern = pattern, full.names = TRUE,
            ignore.case = FALSE)
        if (length(hits) == 0L)
            stop(sprintf("No file matching '%s' found in: %s", pattern, folder),
                call. = FALSE)
        # deterministic order: by filename ascending so _01 precedes _02 precedes _03 ...
        hits <- hits[order(basename(hits))]
        vcf_paths <- hits
    }

    # When called internally from read.gvr.dual() (via the package's force-annotator
    # flag consumed below), suppress the file-listing message: the outer read.gvr()
    # already printed it before routing here, so re-printing would just duplicate.
    if (verbose && !.gvr_force_active()) {
        message(sprintf("Found %d file(s):", length(vcf_paths)))
        for (p in vcf_paths) message("  - ", basename(p))
    }
    # ==========================================================================
    # 0b. Auto-detect annotator per file + dispatch
    #     read.gvr() handles VEP-annotated VCFs (its original role); SnpEff
    #     batches are delegated to read.gvr.snpeff(). In a package, the sibling
    #     is always available; in standalone mode, the auto-source above loads it.
    # ==========================================================================
    detected <- vapply(vcf_paths, .gvr_detect_annotator, character(1L))
    # Internal override: read.gvr.dual() sets a flag in the package's private
    # env via .gvr_set_force_annotator() before calling read.gvr(), then this
    # block consumes-and-clears the flag. The clear-on-read guarantees the
    # override cannot leak into a subsequent user-level read.gvr() call even
    # if read.gvr.dual() errors before returning. Accepts "vep" or "snpeff".
    .fa <- .gvr_consume_force_annotator()
    if (!is.null(.fa)) {
        if (!.fa %in% c("vep", "snpeff"))
            stop("read.gvr: internal force-annotator flag must be 'vep' or 'snpeff'", call. = FALSE)
        detected <- rep(.fa, length(detected))
    }
    na_mask  <- is.na(detected)
    if (any(na_mask)) {
        ## Surface up to 7 INFO tags from the first un-annotated file so the user
        ## can see what *was* annotated (e.g. raw GATK output vs. exotic annotator).
        ## .gvr_detect_annotator() returns NA with attr 'info_tags' attached.
        .unann_first <- vcf_paths[na_mask][1L]
        .tags_first  <- attr(.gvr_detect_annotator(.unann_first), "info_tags")
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
            file              = file,                  # vN+4
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
            genes             = genes,  # already unioned with resolved `panel` above
            panel             = panel,  # Phase N+1: pass through so sibling can run its
            # own "panel subset:" verbose block; sibling
            # re-resolves (idempotent: same registry, same
            # union math, returns the same effective_genes).
            vc_nonSyn         = vc_nonSyn,
            canonical_only    = canonical_only,        # vN+4
            ncores            = ncores,
            verbose           = verbose
        ))
    }
    if (uniq_ann == "dual") {
        if (!exists("read.gvr.dual", mode = "function", inherits = TRUE)) {
            stop("read.gvr: detected dual-annotated VCF(s) (VEP CSQ + SnpEff ANN) ",
                "but 'read.gvr.dual' is not defined. Source 'read.gvr.dual.R' ",
                "(placed next to this file) first.", call. = FALSE)
        }
        if (verbose) message("read.gvr: dual-annotated input detected (VEP CSQ + SnpEff ANN); delegating to read.gvr.dual().")
        return(read.gvr.dual(
            folder            = folder,
            vcf_path          = vcf_path,
            file              = file,
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
            panel             = panel,
            vc_nonSyn         = vc_nonSyn,
            canonical_only    = canonical_only,
            ncores            = ncores,
            verbose           = verbose
        ))
    }
    # else uniq_ann == "vep": continue with the existing VEP body below.


    # ==========================================================================
    # 0c. Resolve ncbi_build: auto-detect from VCF header(s), honour user override.
    #     (Phase N+2)
    #     - ncbi_build = "auto" (default): inspect first file's header via
    #       .gvr_detect_build(); use detected canonical label; fall back to "GRCh38"
    #       when detection is ambiguous; emit a verbose diagnostic either way.
    #     - ncbi_build = <any literal> (override): use the value verbatim; warn()
    #       loudly if auto-detection finds a CONFIDENT canonical label that
    #       disagrees with the user-supplied value (catches typos).
    # ==========================================================================
    {
        .bd_first <- vcf_paths[1L]
        .bd <- tryCatch(.gvr_detect_build(.bd_first),
            error = function(e) list(label = NA_character_,
                confidence = "none",
                signals = list(vep = NA_character_,
                    snpeff = NA_character_,
                    contig = NA_character_)))
        .bd_label <- .bd$label
        .bd_conf  <- .bd$confidence
        .bd_sigs  <- .bd$signals
        .bd_present <- vapply(.bd_sigs, function(x) !is.na(x), logical(1L))
        .bd_sig_str <- if (any(.bd_present))
            paste(names(.bd_sigs)[.bd_present], collapse = " + ") else "none"

        if (identical(ncbi_build, "auto")) {
            if (!is.na(.bd_label)) {
                effective_ncbi_build <- .bd_label
                if (verbose) message(sprintf(
                    "genome build: detected %s (%s agree, %s confidence).",
                    .bd_label, .bd_sig_str, .bd_conf))
            } else {
                effective_ncbi_build <- "GRCh38"
                if (verbose) message(
                    "genome build: could not detect from VCF header (signals: ",
                    .bd_sig_str, "); defaulting to 'GRCh38'. Pass ncbi_build= to override.")
            }
        } else {
            effective_ncbi_build <- ncbi_build
            # Map user value through the same alias table the detector uses, so we
            # treat 'hg19'/'hg38' as canonical equivalents and only warn when the
            # user supplied a build-aware label that contradicts a high-confidence
            # detection. Arbitrary lab codes (no alias hit) pass silently.
            .alias <- c(GRCh38 = "GRCh38", GRCh37 = "GRCh37",
                `T2T-CHM13v2.0` = "T2T-CHM13v2.0",
                hg38 = "GRCh38", hg19 = "GRCh37")
            .user_canon <- if (ncbi_build %in% names(.alias)) unname(.alias[ncbi_build]) else NA_character_
            if (!is.na(.bd_label) && identical(.bd_conf, "high") &&
                !is.na(.user_canon) && !identical(.user_canon, .bd_label)) {
                warning(sprintf(
                    "read.gvr: ncbi_build='%s' but VCF appears to be %s (%s agree). Continuing with user-supplied value.",
                    ncbi_build, .bd_label, .bd_sig_str), call. = FALSE)
            } else if (verbose) {
                message(sprintf(
                    "genome build: '%s' (user-supplied%s).",
                    ncbi_build,
                    if (!is.na(.bd_label) && !is.na(.user_canon) && identical(.user_canon, .bd_label))
                        "; agrees with auto-detection"
                    else if (!is.na(.bd_label) && is.na(.user_canon))
                        "; not in detection table"
                    else if (!is.na(.bd_label))
                        sprintf("; auto-detection found %s [%s]", .bd_label, .bd_conf)
                    else
                        "; auto-detection could not corroborate"))
            }
        }
        # The downstream row-construction site reads `ncbi_build` from lexical
        # scope; rebind so a single call site (`NCBI_Build = ncbi_build` below)
        # remains the canonical write -- no second NCBI_Build= reference needed.
        ncbi_build <- effective_ncbi_build
    }



    # ==========================================================================
    # Per-call memoisation caches for consequence-resolution helpers (Turn 5a).
    # The helpers themselves now live in R/read_gvr_consequence.R and accept
    # `cache` as an argument; we create fresh envs here so the cache lifetime
    # is scoped to a single read.gvr() invocation (identical to pre-Turn-5a
    # behaviour).
    # ==========================================================================
    .gvr_mstr_cache <- new.env(parent = emptyenv())
    .gvr_v2m_cache  <- new.env(parent = emptyenv())

    ## 1i. Convert one chunk (data.table of raw VCF rows) -> table rows
    ## C3 (conservative): the per-record control flow is UNCHANGED (the CSQ-block -> ALT
    ## assignment and per-ALT selection are genuinely irregular and were measured to be
    ## slower, not faster, when forced through flat strsplit + regroup). The one safe,
    ## byte-identical structural change kept here is pre-extracting the 10 VCF columns to
    ## plain atomic vectors ONCE per chunk, so the hot loop indexes vectors (chrom_v[r])
    ## instead of doing data.table `$`/is.factor dispatch on every field of every record
    ## (~4.7x faster for the indexing step in isolation). Output is identical.
    ## 1j. Convert ONE chunk of a VCF file to a MAF-like data.table.
    ##
    ## Turn-5b refactor: the per-chunk parsing pipeline (originally a single 454-line
    ## function body) was lifted into 5 themed package-internal helpers in
    ## `read_gvr_chunk.R`: setup_columns, resolve_format, filter_dpgq, filter_rough,
    ## build_record. The orchestrator here calls them in order, threading a `state`
    ## list of chunk-level vectors through each filter, then iterates over surviving
    ## records to produce per-ALT MAF rows. PERF tricks (O1/O2/O6/O7/O8/A3/P1/P2)
    ## and `# why:` comments are preserved verbatim inside the helpers. Numeric
    ## equality vs Turn-5a is gated by `data.table::all.equal()`; the hero PDF
    ## SHA256 gate is the final byte-identity check.
    ##
    ## Closure captures from `read.gvr()` scope:
    ##   - filter_dp, filter_gq, min_DP, min_GQ, genes, vc_nonSyn, canonical_only,
    ##     ncbi_build (the 8 user-facing args)
    ##   - .gvr_mstr_cache, .gvr_v2m_cache (per-call cache envs, defined at L583-584)
    convert_chunk <- function(dt, csq_fields, sample_name) {
        # 1) Lift VCF columns + cache CSQ field positions (C3 + O2 + ci()).
        state <- .gvr_chunk_setup_columns(dt, csq_fields)
        # 2) Resolve FORMAT/SAMPLE -> chunk-level GT/AD/DP/GQ vectors when constant (O1).
        state <- utils::modifyList(state, .gvr_chunk_resolve_format(state$FORMAT_v, state$SAMPLE_v))
        # 3) Vectorized DP/GQ pre-filter (O6).
        dpgq <- .gvr_chunk_filter_dpgq(state, filter_dp, filter_gq, min_DP, min_GQ)
        state <- dpgq$state
        # 4) Rough gene + vc_nonSyn pre-filters on raw INFO (O7 + O8).
        rough <- .gvr_chunk_filter_rough(state, genes, vc_nonSyn)
        state <- rough$state

        n_dropped <- dpgq$n_dropped_dpgq + rough$n_dropped_gene + rough$n_dropped_vc
        n_dropped_canonical <- 0L          # vN+4: rows skipped by canonical_only=TRUE
        out <- vector("list", state$nrow_dt)
        oi <- 0L

        # 5) Per-record loop. Scalar per-record extraction stays inline (it consumes
        # both the chunk state and the loop index r); the CSQ-block expansion + per-ALT
        # row build is delegated to .gvr_chunk_build_record().
        fmt_constant <- state$fmt_constant
        for (r in seq_len(state$nrow_dt)) {
            # O1: GT/AD/DP/GQ extraction. Fast path uses the chunk-level position split when
            # FORMAT is constant; otherwise the exact original per-record name-based path.
            # In BOTH paths, a missing field resolves to the SAME defaults the original used
            # (gt="./.", others NA_character_). For the fast path, NA from tstrsplit (SAMPLE
            # shorter than the position) is treated as "field absent" -> default, which is
            # identical to the original `seq_along(smp_vals)` name-assignment behaviour.
            if (fmt_constant) {
                gt0  <- if (is.null(state$GT_col)) NA_character_ else state$GT_col[r]
                gt   <- if (is.na(gt0)) "./." else gt0
                ad   <- if (is.null(state$AD_col)) NA_character_ else state$AD_col[r]
                sdp  <- if (is.null(state$DP_col)) NA_character_ else state$DP_col[r]
                gq   <- if (is.null(state$GQ_col)) NA_character_ else state$GQ_col[r]
            } else {
                fmt_keys <- strsplit(state$FORMAT_v[r], ":", fixed = TRUE)[[1]]
                smp_vals <- strsplit(state$SAMPLE_v[r], ":", fixed = TRUE)[[1]]
                names(smp_vals) <- fmt_keys[seq_along(smp_vals)]
                gt  <- if ("GT" %in% names(smp_vals)) smp_vals[["GT"]] else "./."
                ad  <- if ("AD" %in% names(smp_vals)) smp_vals[["AD"]] else NA_character_
                sdp <- if ("DP" %in% names(smp_vals)) smp_vals[["DP"]] else NA_character_
                gq  <- if ("GQ" %in% names(smp_vals)) smp_vals[["GQ"]] else NA_character_
            }
            ad_vec <- if (!is.na(ad)) strsplit(ad, ",", fixed = TRUE)[[1]] else NA

            record_ctx <- list(
                chrom = state$CHROM_v[r], pos = state$POS_v[r], vid = state$ID_v[r],
                ref = state$REF_v[r], altf = state$ALT_v[r], qual = state$QUAL_v[r],
                filt = state$FILTER_v[r], info = state$INFO_v[r],
                alts = state$alts_all[[r]],                          # O2: pre-split once per chunk
                gt = gt, ad = ad, sdp = sdp, gq = gq, ad_vec = ad_vec,
                sample_name = sample_name
            )
            rec <- .gvr_chunk_build_record(
                record_ctx, csq_fields, state$n_csq, state$P,
                canonical_only, ncbi_build,
                .gvr_mstr_cache, .gvr_v2m_cache
            )
            n_dropped_canonical <- n_dropped_canonical + rec$n_dropped_canonical
            if (length(rec$rows) > 0L) {
                for (row in rec$rows) {
                    oi <- oi + 1L
                    out[[oi]] <- row
                }
            }
        }

        res <- data.table::rbindlist(out[seq_len(oi)], use.names = TRUE, fill = TRUE)
        data.table::setattr(res, "n_dropped", n_dropped)   # v4: carry filter drop count to caller
        data.table::setattr(res, "n_dropped_canonical", n_dropped_canonical)  # vN+4
        res
    }

    ## 1j-v2. Convert ONE vcf file end-to-end (header parse + chunked streaming).
    ##        Returns a per-file MAF-like data.table. Progress printed if verbose.
    convert_one_vcf <- function(path, file_idx, file_total, file_total_records = NA_integer_) {
        csq_fields  <- .gvr_get_csq_fields(path)
        sample_name <- .gvr_get_sample_name(path)
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
        repeat {
            l <- readLines(con, 1L)
            if (length(l) == 0L || startsWith(l, "#CHROM")) break
        }

        chunks <- list()
        ci2 <- 0L
        total_in <- 0L
        n_drop_file <- 0L
        n_drop_canon_file <- 0L  # vN+4
        t0 <- Sys.time()
        file_done <- 0L   # per-file running record counter (numerator of progress line)
        repeat {
            lines <- readLines(con, chunk_size)
            if (length(lines) == 0L) break
            ci2 <- ci2 + 1L
            mat <- data.table::tstrsplit(lines, "\t", fixed = TRUE)
            dtc <- data.table::data.table(CHROM = mat[[1]], POS = mat[[2]], ID = mat[[3]], REF = mat[[4]],
                ALT = mat[[5]], QUAL = mat[[6]], FILTER = mat[[7]], INFO = mat[[8]],
                FORMAT = mat[[9]], SAMPLE = mat[[10]])
            ck <- convert_chunk(dtc, csq_fields, sample_name)
            nd <- attr(ck, "n_dropped")
            if (is.null(nd)) nd <- 0L
            n_drop_file <- n_drop_file + nd
            ndc <- attr(ck, "n_dropped_canonical")
            if (is.null(ndc)) ndc <- 0L  # vN+4
            n_drop_canon_file <- n_drop_canon_file + ndc
            chunks[[ci2]] <- ck
            total_in <- total_in + nrow(dtc)
            file_done <- file_done + nrow(dtc)        # per-file running count (numerator)
            .progress$global_done <- .progress$global_done + nrow(dtc)  # shared global counter (drives global %)
            if (verbose) {
                # Per-file count (file_done / this file's total), GLOBAL percentage
                # (global_done / grand_total), GLOBAL elapsed (since whole-run start, t_global).
                el <- as.numeric(difftime(Sys.time(), t_global, units = "secs"))
                if (!is.na(file_total_records) && !is.na(grand_total) && grand_total > 0L) {
                    message(sprintf("    %s/%s (%.0f%%) | %.0fs",
                        format(file_done, big.mark = ","),
                        format(file_total_records, big.mark = ","),
                        100 * .progress$global_done / grand_total, el))
                } else {
                    message(sprintf("    %s records | %.0fs",
                        format(file_done, big.mark = ","), el))
                }
            }
        }
        close(con)

        gvr_one <- data.table::rbindlist(chunks, use.names = TRUE, fill = TRUE)
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
            if (canonical_only && n_drop_canon_file > 0L) {
                # vN+4: emitted rows so far = nrow(gvr_one); dropped + emitted is the
                # pre-filter ALT-row count we would have emitted at this file.
                pre <- nrow(gvr_one) + n_drop_canon_file
                message(sprintf(
                    "    canonical_only: dropped %s of %s ALT-row(s) (%.1f%% non-canonical).",
                    format(n_drop_canon_file, big.mark = ","),
                    format(pre,              big.mark = ","),
                    if (pre > 0L) 100 * n_drop_canon_file / pre else 0))
            }
        }
        attr(gvr_one, "sample_name")          <- sample_name
        attr(gvr_one, "n_dropped")            <- n_drop_file
        attr(gvr_one, "n_dropped_canonical")  <- n_drop_canon_file    # vN+4
        gvr_one
    }

    # ==========================================================================
    # 2. Convert every file and ROW-BIND into one combined gvr table
    # ==========================================================================
    t_all <- Sys.time()

    ## 2a. Progress display in COUNT-ONLY mode (no pre-count pass).
    ##     PERF (v5): the previous implementation fully decompressed every gzip ONCE
    ##     just to count post-#CHROM variant lines (to drive a global % bar), then
    ##     convert_one_vcf decompressed each file AGAIN to convert it -> every file
    ##     was read twice. On the 2-file test set the redundant pre-count pass cost
    ##     ~34.5 s (~8% of runtime). Since the count is used ONLY for the cosmetic
    ##     progress percentage (it never touches the table), we drop it and let the
    ##     conversion report a cumulative running record count + elapsed time
    ##     instead of a %. This is byte-identical in OUTPUT; only verbose progress
    ##     text changes (and only when verbose=TRUE). convert_one_vcf already
    ##     handles a NA per-file total via its count-only branch.
    per_file_total <- rep(NA_integer_, length(vcf_paths))
    grand_total    <- NA_integer_

    if (verbose)
        message(sprintf("Converting %d file(s).", length(vcf_paths)))

    # global running counter + timer shared with convert_one_vcf via an explicit
    # environment (was: <<- on a parent-scope variable; BiocCheck-preferred idiom).
    # Fork-local in parallel mode (each mclapply child gets its own copy).
    .progress <- new.env(parent = emptyenv())
    .progress$global_done <- 0L
    t_global    <- t_all

    # v6 (PERF): optional multi-core conversion ACROSS FILES. Each VCF is converted
    # independently (its own table chunk) before the final rbindlist, so files are
    # embarrassingly parallel. ncores>1 forks parallel::mclapply over the file list;
    # mclapply preserves input order, so `per_file` (and thus the combined gvr table) is
    # byte-identical to the sequential path regardless of ncores. Only effective with
    # >1 file and a fork-capable OS; otherwise we run the original sequential loop.
    # Default ncores=1L => exactly the previous behaviour. The shared .progress env
    # is fork-local (per child), so per-chunk % lines are suppressed
    # in workers to avoid garbled interleaved output. In addition (vN+12), a small
    # side-fork prints a heartbeat line every `heartbeat_secs` while mclapply runs;
    # this gives visible liveness feedback ("workers still running"). The side fork
    # is killed unconditionally on mclapply return, and the per-file `done:` summary
    # is emitted in deterministic index order. parallel is base R (no new dependency).
    n_files  <- length(vcf_paths)
    want_par <- is.numeric(ncores) && length(ncores) == 1L && !is.na(ncores) &&
        ncores > 1L && n_files > 1L &&
        .Platform$OS.type == "unix" &&
        requireNamespace("parallel", quietly = TRUE)
    use_cores <- if (want_par)
        max(1L, min(as.integer(ncores), parallel::detectCores(), n_files))
    else 1L

    if (want_par && use_cores > 1L) {
        # vN+12 (parallel heartbeat): keep the original parallel::mclapply() call
        # for the actual work — it's correct, ordered, and reaps its children
        # cleanly — and add a side-fork "heartbeat" process that prints a single
        # liveness line every `heartbeat_secs` while mclapply runs. The side fork
        # is detached via parallel::mcparallel() and unconditionally killed via
        # tools::pskill() once mclapply returns; it has no dependence on the
        # work results, only on wall-clock. This gives the user visible feedback
        # ("R is alive, parallel conversion still running") without changing the
        # data path: the per-file `done:` summary lines still print AFTER mclapply
        # returns, in the original index order.
        heartbeat_secs <- 15L
        if (verbose)
            message(sprintf("Parallel conversion: %d worker(s) over %d file(s) (mclapply + heartbeat).",
                use_cores, n_files))

        hb_job <- NULL
        if (verbose) {
            t_hb_start <- Sys.time()
            n_hb       <- n_files
            cores_hb   <- use_cores
            hb_job <- parallel::mcparallel(
                {
                    # heartbeat fork: independent of work; prints lines, exits on kill.
                    # Note: messages from this fork go to the same stderr as the parent;
                    # we accept a small interleave risk for liveness signal only.
                    repeat {
                        Sys.sleep(heartbeat_secs)
                        el <- as.numeric(difftime(Sys.time(), t_hb_start, units = "secs"))
                        message(sprintf(
                            "  ... parallel conversion still running: %d worker(s), %.0fs elapsed",
                            cores_hb, el))
                    }
                },
                name = "read_gvr_heartbeat",
                silent = FALSE)
        }

        # actual work: same call as the previous version of read.gvr
        per_file <- tryCatch(
            parallel::mclapply(
                seq_along(vcf_paths),
                # why: suppressMessages on the worker call so per-chunk progress lines from convert_one_vcf() do not interleave with the parent heartbeat output.
                function(i) suppressMessages(
                    convert_one_vcf(vcf_paths[i], i, n_files, per_file_total[i])),
                mc.cores = use_cores, mc.preschedule = FALSE),
            finally = {
                # always reap the heartbeat fork; suppress expected warnings
                # ("1 parallel job did not deliver a result") since we just killed it.
                if (!is.null(hb_job)) {
                    try(tools::pskill(hb_job$pid), silent = TRUE)
                    # why: parallel::mccollect on an already-dead heartbeat fork can warn about missing children; the try() makes the cleanup best-effort.
                    suppressWarnings(try(
                        parallel::mccollect(hb_job, wait = FALSE, timeout = 1), silent = TRUE))
                }
        })

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
            message(sprintf("[file %d/%d] done: %s (sample: %s, %s rows)",
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
            paste(unique(dup_samp), collapse = ", ")))

    gvr <- data.table::rbindlist(per_file, use.names = TRUE, fill = TRUE)
    if (verbose) message(sprintf("COMBINED: %d file(s) | %d rows (%d cols) in %.0fs",
        length(vcf_paths), nrow(gvr), ncol(gvr),
        as.numeric(difftime(Sys.time(), t_all, units = "secs"))))
    # vN+4: aggregated canonical_only summary across all files.
    if (canonical_only) {
        n_can_drop_total <- sum(vapply(per_file, function(x) {
            v <- attr(x, "n_dropped_canonical")
            if (is.null(v)) 0L else as.integer(v)
        }, integer(1L)))
        if (verbose && n_can_drop_total > 0L) {
            pre_total <- nrow(gvr) + n_can_drop_total
            message(sprintf(
                "canonical_only: dropped %s of %s ALT-row(s) total (%.2f%% non-canonical). Pass canonical_only=FALSE to disable.",
                format(n_can_drop_total, big.mark = ","),
                format(pre_total,        big.mark = ","),
                if (pre_total > 0L) 100 * n_can_drop_total / pre_total else 0))
        }
        data.table::setattr(gvr, "n_dropped_canonical", n_can_drop_total)
    }

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
        cn <- names(gvr)
        dup_names <- unique(cn[duplicated(cn)])
        dropped <- character(0)
        renamed <- character(0)
        # Vectorized URL-decode (fast; mirrors the per-element .gvr_url_decode helper).
        url_decode_vec <- function(x) {
            x <- gsub("%3D", "=", x, fixed = TRUE)
            x <- gsub("%3B", ";", x, fixed = TRUE)
            x <- gsub("%2C", ",", x, fixed = TRUE)
            x <- gsub("%3A", ":", x, fixed = TRUE)
            x
        }
        for (nm in dup_names) {
            idx <- which(cn == nm)
            keep <- idx[1]
            # Columns are "identical" iff, after URL-decoding and treating ""/NA as the
            # SAME 'missing' token, every row matches. This is the correct equivalence:
            # the MAF-like-core copy writes "" for missing while the raw CSQ copy writes NA,
            # so a naive (a==b) leaves NA at those rows and must not count as a difference.
            identical_to_keep <- function(j) {
                a <- url_decode_vec(gvr[[keep]])
                b <- url_decode_vec(gvr[[j]])
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
                gvr[, (to_drop) := NULL]
                dropped <- c(dropped, nm)
                cn <- names(gvr)  # refresh after deletion
            }
        }
        # rename any surviving same-named-but-differing duplicates to avoid collision
        cn <- names(gvr)
        if (any(duplicated(cn))) {
            for (nm in unique(cn[duplicated(cn)])) {
                idx <- which(cn == nm)
                for (k in seq_along(idx)[-1]) {
                    new <- paste0(nm, ".csq", if (k > 2) k - 1 else "")
                    data.table::setnames(gvr, idx[k], new)
                }
                cn <- names(gvr)
            }
        }
    }

    ## 3b. Strip Ensembl feature prefix from HGVSc / HGVSp ---------------------
    ##     "ENST...:n.1889G>A" -> "n.1889G>A". Applied to the surviving core cols.
    if (strip_hgvs_prefix) {
        # Vectorized strip: remove up to first ':'. Empty/NA pass through unchanged.
        strip_prefix_vec <- function(x) ifelse(is.na(x) | x == "", x, sub("^[^:]*:", "", x))
        for (col in c("HGVSc", "HGVSp")) {
            if (col %in% names(gvr)) {
                data.table::set(gvr, j = col, value = strip_prefix_vec(gvr[[col]]))   # in-place, no copy
            }
        }
    }

    ## 3c. Add Genotype = Tumor_Seq_Allele1 / Tumor_Seq_Allele2 ----------------
    if (add_genotype) {
        data.table::setDT(gvr)   # ensure a clean data.table (prior := NULL may have left a shallow copy)
        if (all(c("Tumor_Seq_Allele1", "Tumor_Seq_Allele2") %in% names(gvr))) {
            gvr[, Genotype := paste(Tumor_Seq_Allele1, Tumor_Seq_Allele2, sep = "/")]
            # place Genotype directly after Tumor_Seq_Allele2
            nm <- names(gvr)
            after <- which(nm == "Tumor_Seq_Allele2")
            neworder <- append(setdiff(nm, "Genotype"), "Genotype", after = after)
            data.table::setcolorder(gvr, neworder)
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
        data.table::setDT(gvr)
        gvr[, ABraOM_AF := ""]                  # default blank (no-match convention)
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
                ok <- tryCatch(
                    {
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
                    },
                    error = function(e) {
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
            # why: data.table::fread() may warn about ragged ABraOM rows or quoted whitespace; we tolerate that, and a fatal error is caught by tryCatch.
            ab <- tryCatch(suppressWarnings(data.table::fread(apath, sep = "\t", header = TRUE, quote = "",
                showProgress = FALSE)), error = function(e) NULL)
            if (is.null(ab)) {
                warning("read.gvr: ABraOM reference unreadable; 'ABraOM_AF' left blank.",
                    call. = FALSE)
            } else {
                data.table::setnames(ab, make.names(names(ab)))
                need <- c("avsnp147", "Ref", "Alt", "Frequencies")
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
                    gvr[, `:=`(.rs  = dbSNP_RS,
                        .ref = toupper(Reference_Allele),
                        .alt = toupper(Tumor_Seq_Allele2))]
                    hit <- lut[gvr[, .(.rs, .ref, .alt)],
                        on = c(rs = ".rs", ref = ".ref", alt = ".alt"),
                        x.af]            # vector of matched AF, NA where no match
                    gvr[, ABraOM_AF := ifelse(is.na(hit), "", hit)]
                    gvr[, c(".rs", ".ref", ".alt") := NULL]
                    # Place the new column with the other population-frequency columns:
                    # right after MAX_AF_POPS / MAX_AF / the last gnomAD column (whichever
                    # exists). This groups it with gnomAD/1000G AFs, NOT the INFO_* block.
                    nm <- names(gvr)
                    anchors <- c(grep("^MAX_AF_POPS$", nm), grep("^MAX_AF$", nm),
                        grep("^gnomAD", nm))
                    after <- if (length(anchors)) max(anchors) else length(nm) - 1L
                    neworder <- append(setdiff(nm, "ABraOM_AF"),
                        "ABraOM_AF", after = after)
                    data.table::setcolorder(gvr, neworder)
                    if (verbose) message(sprintf(
                        "ABraOM: annotated %d/%d rows (%.1f%%) with ABraOM frequency (rsID+allele join).",
                        sum(nzchar(gvr$ABraOM_AF)), nrow(gvr),
                        100 * mean(nzchar(gvr$ABraOM_AF))))
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
        data.table::setDT(gvr)
        genes_chr <- as.character(genes)
        genes_chr <- genes_chr[!is.na(genes_chr) & nzchar(genes_chr)]
        want <- unique(toupper(trimws(genes_chr)))
        n_before <- nrow(gvr)
        have <- toupper(trimws(as.character(gvr$Hugo_Symbol)))
        gvr <- gvr[have %in% want]
        if (verbose) {
            found    <- intersect(want, unique(have))
            notfound <- setdiff(want, unique(have))
            message(sprintf(
                "gene subset: kept %s / %s rows across %d / %d requested gene(s)%s.",
                format(nrow(gvr), big.mark = ","), format(n_before, big.mark = ","),
                length(found), length(want),
                if (length(notfound)) paste0(" (not found: ", paste(notfound, collapse = ", "), ")") else ""))
        }
    }

    ## 3d-bis2. Phase N+1: post-filter panel coverage message -----------------
    ##     Sibling of "gene subset:" (3d-bis). Fires only when `panel` was
    ##     supplied AND resolved to >=1 gene (gated on
    ##     exists(.panel_genes_for_summary), which was set inside the vN setup
    ##     block at the top of the body). Reports how many of the ORIGINAL
    ##     (pre-union) panel genes are present in the post-filter result,
    ##     lists present + missing, and (when `panel` was combined with
    ##     explicit `genes` extras) names them inline. Observation point
    ##     matches "gene subset:" semantics (pre vc_nonSyn / drop_empty_cols).
    if (verbose && exists(".panel_genes_for_summary", inherits = FALSE)) {
        data.table::setDT(gvr)
        have_syms <- unique(toupper(trimws(as.character(gvr$Hugo_Symbol))))
        p_present <- intersect(.panel_genes_for_summary, have_syms)
        p_missing <- setdiff(.panel_genes_for_summary, have_syms)
        extras_chunk <- if (length(.panel_extras_for_summary) > 0L) {
            sprintf(" (+ %d extra gene(s) from `genes`: %s)",
                length(.panel_extras_for_summary),
                paste(.panel_extras_for_summary, collapse = ", "))
        } else ""
        missing_chunk <- if (length(p_missing) > 0L) {
            sprintf(" (missing: %s)", paste(p_missing, collapse = ", "))
        } else ""
        message(sprintf(
            "panel subset: %d / %d panel gene(s) present: %s%s%s.",
            length(p_present), length(.panel_genes_for_summary),
            if (length(p_present) > 0L) paste(p_present, collapse = ", ") else "(none)",
            missing_chunk, extras_chunk))
    }

    ## 3d-ter. Variant-Classification non-synonymous filter (opt-in) -----------
    ##     Mirrors the non-synonymous-only filter: keep only protein-altering classes.
    ##     FALSE (default) = keep all. TRUE = keep the 9 standard High/Moderate
    ##     classes. A custom character vector = keep those specific classes.
    ##     Applied AFTER gene subset, BEFORE drop_empty_cols.
    if (!identical(vc_nonSyn, FALSE)) {
        data.table::setDT(gvr)
        vc_default <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site",
            "Translation_Start_Site", "Nonsense_Mutation", "Nonstop_Mutation",
            "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation")
        vc_keep <- if (isTRUE(vc_nonSyn)) vc_default else as.character(vc_nonSyn)
        vc_keep <- vc_keep[!is.na(vc_keep) & nzchar(vc_keep)]
        if (length(vc_keep) > 0L) {
            n_before <- nrow(gvr)
            vc_col <- as.character(gvr$Variant_Classification)
            # Remove rows with missing/blank Variant_Classification when filter is active
            gvr <- gvr[!is.na(vc_col) & nzchar(vc_col) & vc_col %in% vc_keep]
            if (verbose) {
                dropped <- n_before - nrow(gvr)
                message(sprintf(
                    "vc_nonSyn: kept %s / %s rows (%d classification(s): %s; removed %s silent/other)",
                    format(nrow(gvr), big.mark = ","), format(n_before, big.mark = ","),
                    length(vc_keep), paste(vc_keep, collapse = ", "),
                    format(dropped, big.mark = ",")))
            }
        }
    }

    ## 3e. Drop all-empty columns (opt-in) ------------------------------------
    ##     Remove any column whose values are ALL missing (NA or "") across the
    ##     entire combined table. Default FALSE = keep full schema.
    if (drop_empty_cols) {
        data.table::setDT(gvr)
        is_all_empty <- function(col) {
            v <- gvr[[col]]
            vc <- as.character(v)
            all(is.na(v) | vc == "")
        }
        empty_cols <- names(gvr)[vapply(names(gvr), is_all_empty, logical(1))]
        if (length(empty_cols) > 0L) {
            gvr[, (empty_cols) := NULL]
            if (verbose) message(sprintf("drop_empty_cols: removed %d all-empty column(s): %s",
                length(empty_cols), paste(empty_cols, collapse = ", ")))
        } else if (verbose) {
            message("drop_empty_cols: no all-empty columns found.")
        }
    }

    if (verbose) message(sprintf("Final Table Dimensions: %d rows x %d columns.", nrow(gvr), ncol(gvr)))

    # Tag annotator BEFORE writes so saved TSV/RDS/XLSX carry the attribute too.
    # (R's saveRDS preserves attributes; for TSV/XLSX it's only for in-memory
    #  consumers -- but we still keep it consistent.)
    data.table::setattr(gvr, "annotator", "vep")

    # ==========================================================================
    # 4. Optional file outputs
    # ==========================================================================
    if (write_tsv || write_rds || write_xlsx) {
        # vN+4: handle vector vcf_path (use first), file= mode (folder), and folder mode (folder)
        if (is.null(out_dir)) {
            out_dir <- if (!is.null(vcf_path)) dirname(vcf_path[1L]) else folder
        }
        if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        if (is.null(out_prefix)) {
            if (length(vcf_paths) > 1L) {
                out_prefix <- "combined.gvr.tsv"                              # multi-file
            } else {
                out_prefix <- sub("\\.vcf\\.gz$", "", basename(vcf_paths[1])) # single file
                out_prefix <- paste0(out_prefix, ".gvr.tsv")
            }
        }
        if (write_tsv) {
            tsv_path <- file.path(out_dir, paste0(out_prefix, ".tsv"))
            data.table::fwrite(gvr, tsv_path, sep = "\t", quote = FALSE, na = "")
            if (verbose) message("Wrote TSV: ", tsv_path)
        }
        if (write_rds) {
            # IMPORTANT: on S3-backed mounts (e.g. /mnt/results), R's file.copy() can
            # silently produce a 0-byte file. The ONLY reliable route is to saveRDS to
            # a local POSIX path (tempdir / /workspace) and then shell `cp` it over.
            rds_final <- file.path(out_dir, paste0(out_prefix, ".rds"))
            tmp_rds   <- file.path(tempdir(), paste0(out_prefix, ".rds"))
            saveRDS(gvr, tmp_rds, compress = TRUE)
            system2("cp", c(shQuote(tmp_rds), shQuote(rds_final)))   # always shell-cp
            # why: file.info()$size may warn / return NA right after a shell cp on S3 FUSE; the is.na()||size==0 check below handles that.
            sz <- suppressWarnings(file.info(rds_final)$size)
            if (is.na(sz) || sz == 0)
                warning(sprintf("RDS write may have failed (0 bytes): %s", rds_final))
            else if (verbose)
                message(sprintf("Wrote RDS: %s (%.0f MB)", rds_final, sz / 1e6))
        }
        if (write_xlsx) {
            # B1: optional Excel export of the FINAL gvr table (one "gvr_table" sheet). Mirrors the
            # FUSE-safe openxlsx pattern used by gvr_summary: build the workbook, save to a
            # local temp file, then shell-cp to out_dir (openxlsx uses zip random-access
            # writes that can silently 0-byte on S3-backed mounts). Degrades gracefully:
            # if openxlsx is absent we warn and skip (TSV/RDS, if requested, still wrote).
            # NOTE: a germline gvr table can be very large (hundreds of thousands of rows). Excel
            # handles it but the file is big and slow to open; write_rds / write_tsv remain
            # the better choice for large tables / downstream R use.
            if (!requireNamespace("openxlsx", quietly = TRUE)) {
                warning("read.gvr: 'openxlsx' not installed; skipping Excel export.")
            } else {
                xlsx_final <- file.path(out_dir, paste0(out_prefix, ".xlsx"))
                if (file.exists(xlsx_final) && verbose)
                    message(sprintf("  Overwriting existing Excel: %s", xlsx_final))
                if (nrow(gvr) > 1000000L)
                    warning(sprintf("read.gvr: table has %s rows; Excel's per-sheet limit is 1,048,576 rows.",
                        format(nrow(gvr), big.mark = ",")))
                wb <- openxlsx::createWorkbook()
                hs <- openxlsx::createStyle(textDecoration = "bold", halign = "center")
                openxlsx::addWorksheet(wb, "gvr_table")
                openxlsx::writeData(wb, "gvr_table", as.data.frame(gvr), headerStyle = hs)
                openxlsx::freezePane(wb, "gvr_table", firstRow = TRUE)
                openxlsx::setColWidths(wb, "gvr_table", cols = seq_len(ncol(gvr)), widths = "auto")
                tmp_xlsx <- file.path(tempdir(), paste0(out_prefix, ".xlsx"))
                wrote_ok <- tryCatch(
                    {
                        openxlsx::saveWorkbook(wb, tmp_xlsx, overwrite = TRUE)
                        TRUE
                    },
                    error = function(e) {
                        warning(sprintf("read.gvr: Excel write failed: %s", conditionMessage(e)))
                        FALSE
                    })
                if (wrote_ok) {
                    system2("cp", c(shQuote(tmp_xlsx), shQuote(xlsx_final)))
                    # why: same as above — file.info()$size on a freshly-copied XLSX on a slow FS may warn; the is.na()||size==0 fallback below handles it.
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
    data.table::setDT(gvr)
    gvr[]
}

# If sourced interactively this just defines read.gvr().
# Example (commented):
#   gvr <- read.gvr("/mnt/user-uploads")                # merge all *_NN.vcf.gz
#   gvr_pass <- gvr[FILTER == "PASS"]
#   # genes of interest only, no DP/GQ filter:
#   gvr_goi <- read.gvr("/mnt/user-uploads", min_DP = NULL, min_GQ = NULL,
#                       genes = c("MEN1","RET","CDKN1B"))

# NOTE: globalVariables() declarations for this package are consolidated in
# R/globals.R (one package-scoped block covering all functions).
