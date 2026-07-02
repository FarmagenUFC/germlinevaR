# =============================================================================
# read.gvr.snpeff.R
# -----------------------------------------------------------------------------
# SnpEff-annotated VCF reader for the germlinevaR toolchain. Sibling to
# read.gvr() (VEP).
#
# Public function: read.gvr.snpeff(folder, vcf_path, ...) -- same signature
# and defaults as read.gvr(). Returns a data.table with the canonical table
# schema and attr(., 'annotator') = 'snpeff'.
#
# This file defines ONE public function. All helpers (.AF_FIELDS_SNPEFF,
# .CLINSIG_FIELDS_SNPEFF, .detect_annotator, get_ann_fields,
# .snpeff_strip_allele, .CSQ_FIELDS_VEP_CANONICAL, convert_chunk_snpeff)
# are nested inside read.gvr.snpeff() for R-package hygiene.
# =============================================================================


#' Convert SnpEff-annotated germline VCF(s) to a tabular variant data.table
#'
#' @description
#' Converts SnpEff-annotated, single-sample germline VCFs (GATK HaplotypeCaller ->
#' CNN tranches -> SnpEff, hg38) into a tabular variant table and returns it as
#' an in-memory `data.table` for downstream filtering ([gvr_filter()]) and
#' summarisation ([gvr_summary()]). In folder mode it finds every per-sample VCF,
#' converts each, and row-binds them into one combined gvr table. Emits the SAME
#' canonical 80-field schema as [read.gvr()] so downstream code is annotator-
#' agnostic. The conversion uses base R + \pkg{data.table} only - no external
#' annotation-package dependency.
#'
#' This function is the SnpEff sibling of [read.gvr()] (VEP-annotated VCFs).
#' Both readers are usually invoked through [read.gvr()], which auto-routes
#' SnpEff inputs to this function via the nested `.detect_annotator()` helper.
#' Calling `read.gvr.snpeff()` directly is supported for fully SnpEff-only
#' pipelines.
#'
#' @details
#' Output and behaviour:
#' \itemize{
#'   \item Returns the final tabular variant `data.table`, one row per variant ALLELE
#'     (multi-allelic sites are split). Multiple `ANN` annotation blocks per
#'     allele are reduced to one most-severe transcript per allele.
#'   \item Columns include the canonical MAF-style core fields (Hugo_Symbol, Variant_Classification, Start_Position, Reference_Allele, Tumor_Seq_Allele2, Tumor_Sample_Barcode, HGVSp_Short, IMPACT), the canonical 80 VEP CSQ field
#'     names (populated from the equivalent SnpEff ANN fields where available;
#'     blank for fields SnpEff does not supply), and key GATK QC fields.
#'     `FILTER` is retained as a column and ALL variants (PASS and non-PASS) are
#'     kept.
#'   \item `Tumor_Seq_Allele1`/`Tumor_Seq_Allele2` are zygosity-aware (vcf2maf-
#'     style), and an optional `Genotype` column (`Tumor_Seq_Allele1/Tumor_Seq_Allele2`,
#'     e.g. `"T/C"`) is added next to the alleles.
#'   \item Each variant keeps its source sample in `Tumor_Sample_Barcode`.
#'   \item Absent values are written as the empty string `""` (not `NA`);
#'     downstream [gvr_filter()] / [gvr_summary()] treat `NA` and `""` identically
#'     as "missing".
#'   \item Output is tagged with `attr(gvr, "annotator") = "snpeff"` so downstream
#'     code can distinguish the source.
#' }
#'
#' SnpEff-specific notes:
#' \itemize{
#'   \item The nested `.detect_annotator()` helper inspects `##INFO=<ID=*` header
#'     lines and returns `"snpeff"` when `##INFO=<ID=ANN` is found (and `"vep"`
#'     when `##INFO=<ID=CSQ` is found). If both are present, VEP takes priority.
#'   \item The nested `get_ann_fields()` helper parses SnpEff's ANN INFO header
#'     line, extracting the `'F1 | F2 | ... | FN'` field-name vector from inside
#'     the single-quoted `Description=` block.
#'   \item The nested `.snpeff_strip_allele()` helper handles SnpEff's two
#'     non-standard allele forms (Cingolani spec): cancer-somatic-vs-germline
#'     (`"G-C"` -> `"G"`) and compound (`"C-chr1:123456_A>T"` -> `"C"`).
#'     Standard ALT alleles (no dash) pass through unchanged.
#'   \item Population AF columns (`gnomADe_AF`, etc.) are populated from the first
#'     matching SnpEff ANN-side INFO key in `.AF_FIELDS_SNPEFF` (case-sensitive,
#'     first match wins). Bare `AF` is intentionally excluded - in the VCF spec
#'     it is the per-ALT caller frequency, not a population frequency.
#' }
#'
#' Processing options (same as [read.gvr()]):
#' \itemize{
#'   \item MULTI-FILE: in folder mode, every file matching `pattern` (default
#'     `"*_NN.vcf.gz"`) is converted and row-bound.
#'   \item HGVS CLEANUP (`strip_hgvs_prefix`): strips the Ensembl feature prefix
#'     from `HGVSc`/`HGVSp` (e.g. `"ENST00000831140.1:n.1889G>A"` ->
#'     `"n.1889G>A"`).
#'   \item DEDUP (`dedup_columns`): removes duplicate-named columns ONLY when
#'     their values are byte-for-byte identical across all rows (otherwise keeps
#'     + warns).
#'   \item ABraOM (`add_abraom`): joins the Brazilian ABraOM SABE-609 allele
#'     frequency as the `ABraOM_AF` column (downloaded/cached from `abraom_url`).
#'   \item GENOTYPE-QUALITY FILTER (`min_DP`/`min_GQ`): keeps a record iff
#'     `DP > min_DP` AND `GQ > min_GQ`; mirrors
#'     `bcftools view -e 'FORMAT/DP<=X | FORMAT/GQ<=Y'`. Set either to `NULL` to
#'     disable that field; set both to `NULL` to disable the genotype filter
#'     entirely.
#'   \item GENE SUBSET (`genes`): restrict to a set of `Hugo_Symbol`s (exact,
#'     case-insensitive).
#' }
#'
#' @param folder Directory to scan in folder mode; every file matching `pattern` is
#'   converted and row-bound. Default `"."`. Ignored when `vcf_path` is supplied.
#'   Also used as the search root for `file=`.
#' @param vcf_path Character vector of one or more full paths to `.vcf.gz`
#'   files to convert. Mutually exclusive with `file=`. `NULL` (default)
#'   selects folder mode.
#' @param file Character vector of basenames (e.g.
#'   `c("S1.snpeff.vcf.gz", "S2.snpeff.vcf.gz")`) resolved against `folder=`.
#'   Use this to pick specific files from a folder that contains files you
#'   do NOT want to merge. Mutually exclusive with `vcf_path=`. `NULL`
#'   (default) selects either `vcf_path=` mode or folder-pattern mode.
#' @param pattern Regular expression identifying per-sample VCFs in folder mode.
#'   Default `"\\.vcf\\.gz$"` (matches any `.vcf.gz` file). The old default
#'   `"_\\d+(\\.(vep|snp[eE]ff))?\\.vcf\\.gz$"` is still available by passing
#'   it explicitly.
#' @param write_tsv Logical; if `TRUE`, also write the table as a TSV to `out_dir`.
#'   Default `FALSE`.
#' @param write_rds Logical; if `TRUE`, also write the table as an `.rds` to `out_dir`.
#'   Default `FALSE`.
#' @param write_xlsx Logical; if `TRUE`, also write the table as an `.xlsx` workbook
#'   (single `"gvr_table"` sheet) to `out_dir`. Requires the \pkg{openxlsx} package (a
#'   `Suggests` dependency); if it is not installed the export is skipped with a
#'   warning. Default `FALSE`. Excel caps a sheet at 1,048,576 rows.
#' @param out_dir Output directory for written TSV/RDS/XLSX. `NULL` (default) uses
#'   the input location/working directory. Only used when
#'   `write_tsv`/`write_rds`/`write_xlsx` is `TRUE`.
#' @param out_prefix Filename prefix for written outputs. `NULL` (default) derives
#'   one from the input.
#' @param chunk_size Integer; number of VCF records processed per chunk (controls
#'   peak memory and progress granularity). Default `25000L`.
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
#' @param dedup_columns Logical; if `TRUE` (default) drop duplicate-named columns
#'   when byte-for-byte identical (otherwise keep and warn).
#' @param drop_empty_cols Logical; if `TRUE`, drop columns that are entirely
#'   `NA`/blank. Default `FALSE`.
#' @param add_abraom Logical; if `TRUE` (default) join the ABraOM SABE-609 allele
#'   frequency as `ABraOM_AF`.
#' @param abraom_path Path to a local ABraOM annotation file. `NULL` (default) uses
#'   an auto-managed cache (see `cache_dir`), downloading from `abraom_url` if
#'   needed.
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
#' @param hpo Character vector of Human Phenotype Ontology term id(s)
#'   (e.g. `"HP:0003002"`, `c("HP:0003002", "HP:0025022")`). Each id is
#'   resolved to its gene vector via [gvr_hpo_genes()] using the HPO
#'   `phenotype_to_genes.txt` annotation file, and the union of all
#'   resolved genes is taken with `genes` and any panel genes
#'   (deduplicated, uppercased). Lenient input forms are accepted and
#'   normalised: `"HP:0003002"`, `"hp:0003002"`, `"hp:3002"`, `"3002"`,
#'   and `"0003002"` all coerce to canonical `"HP:0003002"`. Only the
#'   exact terms are used; descendants in the ontology are NOT expanded.
#'   The annotation file is fetched once per session and cached under
#'   [tools::R_user_dir()] `"germlinevaR"` `"cache"`; the cache is
#'   auto-refreshed after 30 days. For offline / hermetic runs, point
#'   `options(gvr.hpo_path = )` at a local copy. `NULL` (default)
#'   disables HPO filtering; behaviour is then byte-identical to omitting
#'   the argument.
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
#' @param canonical_only Logical; accepted for API symmetry with [read.gvr()].
#'   SnpEff's ANN field has no `CANONICAL` flag, so this filter cannot be
#'   applied. When `TRUE` (the default), `read.gvr.snpeff()` emits a one-time
#'   warning and returns the unfiltered table (same result as
#'   `canonical_only = FALSE`).
#' @param ncores Integer; number of worker processes for converting MULTIPLE input
#'   files in parallel via [parallel::mclapply()] (fork-based; Unix/macOS only).
#'   Default `1L` runs sequentially and is byte-identical to previous behaviour.
#'   Values `> 1` only help when more than one VCF is being read (each file is an
#'   independent task) and are clamped to `min(ncores, detectCores(), n_files)`. On
#'   non-fork platforms it falls back to sequential. A single file is unaffected.
#' @param normalize_alleles Logical; if `TRUE` (default, since 0.99.2) apply
#'   bcftools-norm-style trimming of common REF/ALT prefix and suffix nucleotides
#'   before deriving the trimmed (Start_Position, Reference_Allele, Tumor_Seq_Allele2) coordinates. See [read.gvr()] for full rationale. Set
#'   `FALSE` to reproduce pre-0.99.2 coords for reproducibility with an older
#'   analysis; not recommended for new research.
#' @param verbose Logical; if `TRUE` (default) print per-file and per-chunk progress
#'   (file i/N, cumulative records, elapsed seconds).
#'
#' @return A tabular variant `data.table`: one row per variant allele, with the canonical variant table columns,
#'   the canonical 80 VEP CSQ field names (populated from SnpEff ANN where
#'   available; blank otherwise), key GATK QC fields, `Tumor_Sample_Barcode`, and
#'   (when enabled) the `Genotype` and `ABraOM_AF` columns.
#'   `attr(., "annotator") = "snpeff"`. TSV/RDS/XLSX files are written as a side
#'   effect when `write_tsv`/`write_rds`/`write_xlsx` is `TRUE`.
#'
#' @seealso [read.gvr()] for the VEP-annotated sibling reader (and the recommended
#'   entrypoint, which auto-routes SnpEff input to this function),
#'   [gvr_filter()] to filter the returned table, [gvr_summary()] to summarise it.
#' @family germlinevaR
#' @author germlinevaR authors
#'
#' @examples
#' ## The function signature is exported and callable:
#' is.function(read.gvr.snpeff)
#'
#' \donttest{
#' ## read.gvr.snpeff() expects VCFs annotated by SnpEff (ANN/LOF/NMD
#' ## INFO fields). The shipped example VCF is VEP-annotated, so a real
#' ## SnpEff example needs your own VCFs; we therefore guard each call on
#' ## the path existing, so the example skips cleanly on machines without
#' ## the data.
#' snpeff_dir <- "/path/to/snpeff-vcfs/"
#' if (dir.exists(snpeff_dir)) {
#'     gvr <- read.gvr.snpeff(snpeff_dir)
#'
#'     ## Or via the auto-router when the VCF header declares SnpEff fields:
#'     gvr <- read.gvr(snpeff_dir)
#' }
#'
#' ## Single-file mode: full path
#' snpeff_vcf <- "/path/to/SAMPLE_01.snpeff.vcf.gz"
#' if (file.exists(snpeff_vcf)) {
#'     gvr <- read.gvr.snpeff(vcf_path = snpeff_vcf)
#' }
#' }
#'
#' @importFrom data.table data.table as.data.table rbindlist fread fwrite setnames setcolorder setDT set setattr tstrsplit setkey :=
#' @importFrom stats setNames
#' @importFrom utils download.file
#' @importFrom openxlsx createWorkbook
#' @export
read.gvr.snpeff <- function(folder = ".",
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
                            hpo               = NULL,   # vN+2: HPO term id(s); union'd with `genes` and `panel`
                            vc_nonSyn         = FALSE,  # v8: keep only protein-altering Variant_Classification
                            canonical_only    = TRUE,   # vN+4: API symmetry; SnpEff ANN has no CANONICAL -> warn and ignore
                            ncores            = 1L,     # v6: parallel files (>1 forks mclapply; 1 = sequential, default)
                            normalize_alleles = TRUE,   # v0.99.2: bcftools-norm-style trim of common prefix/suffix nt before deriving the (Start_Position, Reference_Allele, Tumor_Seq_Allele2) coordinates (recommended); FALSE reproduces pre-0.99.2 coords for reproducibility with older analyses
                            verbose    = TRUE) {
    # ===========================================================================
    # Nested helpers and constants (previously top-level)
    # ---------------------------------------------------------------------------
    # .AF_FIELDS_SNPEFF, .CLINSIG_FIELDS_SNPEFF, .detect_annotator,
    # get_ann_fields, .snpeff_strip_allele, .CSQ_FIELDS_VEP_CANONICAL,
    # convert_chunk_snpeff. Nested for package hygiene; no semantic change.
    # ===========================================================================


    # -----------------------------------------------------------------------------
    # Canonical field-name lists used for INFO scanning when populating
    # gnomADe_AF and CLIN_SIG in SnpEff output. First case-sensitive match wins.
    # If a real-world file uses an exotic name, the column ends up empty (no
    # fabrication).
    # NOTE: bare 'AF' is intentionally excluded -- in VCF spec it is the per-ALT
    # caller allele frequency from the variant caller (GATK/etc.), not a
    # population frequency. Population AF always appears under a namespaced key.
    # -----------------------------------------------------------------------------
    .AF_FIELDS_SNPEFF <- c(
        "gnomADe_AF", "gnomADg_AF", "gnomAD_AF",
        "gnomAD_exomes_AF", "gnomAD_genomes_AF",
        "AF_popmax", "AF_nfe",
        "ExAC_AF", "1000Gp3_AF", "TOPMED_AF"
    )


    .CLINSIG_FIELDS_SNPEFF <- c(
        "CLIN_SIG", "CLNSIG", "ClinVar_CLNSIG",
        "ClinVar_CLNSIGCONF", "CLNSIGCONF"
    )


    # -----------------------------------------------------------------------------
    # .detect_annotator(path)
    # Reads the VCF header lines (until first non-`##` line) and decides whether
    # the file is VEP- or SnpEff-annotated.
    #
    # Returns one of:
    #   "vep"             -- found `##INFO=<ID=CSQ`
    #   "snpeff"          -- found `##INFO=<ID=ANN` (and no CSQ)
    #   NA_character_     -- neither present
    #
    # Always attaches attr(<result>, "info_tags") with all `##INFO=<ID=Xxx` tag
    # names found, so the un-annotated error path can list them in the message.
    # -----------------------------------------------------------------------------
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


    # ===========================================================================
    # Nested helper: .detect_build()  (Phase N+2; byte-identical to the helper
    # in read.gvr.R)
    # ---------------------------------------------------------------------------
    # Reads a VCF header and returns a list:
    #   list(label = <canonical or NA_character_>,
    #        confidence = c("high", "low", "none"),
    #        signals = list(vep = ..., snpeff = ..., contig = ...))
    # Canonical labels: "GRCh38", "GRCh37", "T2T-CHM13v2.0".
    # Sources combined: VEP `##VEP=` `assembly=` attr, SnpEff `##SnpEffCmd=`
    # database token, and the FIRST `##contig=<ID=...>` whose ID matches a
    # standard autosome -- length is compared against a small built-in table.
    # confidence='high' when >=2 signals agree on a canonical label;
    # confidence='low' when exactly one signal is present; 'none' otherwise.
    # ===========================================================================
    .detect_build <- function(path) {
        # Canonical contig lookup. chr1 length is the most discriminating signal:
        #   GRCh38         chr1 = 248,956,422  (NCBI GRCh38.p14)
        #   GRCh37         chr1 = 249,250,621  (NCBI GRCh37.p13)
        #   T2T-CHM13v2.0  chr1 = 248,387,328  (T2T consortium)
        contig_lookup <- list(
            `1` = list(`248956422` = "GRCh38",
                `249250621` = "GRCh37",
                `248387328` = "T2T-CHM13v2.0"),
            `M` = list(`16569` = "GRCh38",  `16571` = "GRCh37")  # chrM tiebreaker
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

            # VEP signal: ##VEP="v113.0" ... assembly="GRCh38.p14" ...
            if (is.na(vep_sig) && startsWith(line, "##VEP=")) {
                m <- regmatches(line, regexec("assembly=\"(GRCh3[78])", line))[[1]]
                if (length(m) == 2L) vep_sig <- m[2]
            }

            # SnpEff signal: ##SnpEffCmd="SnpEff ... GRCh38.99 input.vcf"
            if (is.na(snpeff_sig) && startsWith(line, "##SnpEffCmd=")) {
                # database token is the first GRCh38.x / GRCh37.x / hg38 / hg19 occurrence
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

            # Contig signal: first ##contig=<ID=(chr)?(1|M),length=N>
            if (!first_std_contig_seen && startsWith(line, "##contig=<ID=")) {
                m <- regmatches(line, regexec(
                    "^##contig=<ID=(?:chr)?([0-9XYM]+),length=([0-9]+)", line))[[1]]
                if (length(m) == 3L) {
                    first_std_contig_seen <- TRUE
                    chrom_key <- m[2]
                    len_key   <- m[3]
                    if (!is.null(contig_lookup[[chrom_key]]) &&
                        !is.null(contig_lookup[[chrom_key]][[len_key]])) {
                        contig_sig <- contig_lookup[[chrom_key]][[len_key]]
                    }
                    # else: first std contig found but length doesn't match any known build
                }
            }
        }

        # Combine signals: 2-of-3 agreement => high confidence;
        # 1 signal => low confidence; 0 or all-disagree => none.
        sigs <- c(vep = vep_sig, snpeff = snpeff_sig, contig = contig_sig)
        sigs_present <- sigs[!is.na(sigs)]
        if (length(sigs_present) == 0L) {
            return(list(label = NA_character_, confidence = "none",
                signals = as.list(sigs)))
        }
        tbl <- sort(table(sigs_present), decreasing = TRUE)
        top <- names(tbl)[1L]
        top_count <- as.integer(tbl[1L])
        if (top_count >= 2L) {
            return(list(label = top, confidence = "high",
                signals = as.list(sigs)))
        }
        if (length(sigs_present) == 1L) {
            return(list(label = top, confidence = "low",
                signals = as.list(sigs)))
        }
        # All 2-3 signals present but no two agree -> conflict
        list(label = NA_character_, confidence = "none",
            signals = as.list(sigs))
    }



    # -----------------------------------------------------------------------------
    # get_ann_fields(path)
    # Parses SnpEff's ANN INFO header into a character vector of field names.
    #
    # Canonical SnpEff format (verified against pcingola.github.io):
    #   ##INFO=<ID=ANN,Number=.,Type=String,
    #          Description="Functional annotations: 'F1 | F2 | F3 | ... | FN'">
    #
    # Implementation finds the first and last single-quote in the line, then
    # splits the substring between them on " | " with whitespace trimmed.
    # Robust to any preceding wording before the single-quoted block.
    # -----------------------------------------------------------------------------
    get_ann_fields <- function(path) {
        con <- gzfile(path, "r")
        on.exit(close(con))
        repeat {
            line <- readLines(con, n = 1L, warn = FALSE)
            if (length(line) == 0L) break
            if (!startsWith(line, "##")) break
            if (startsWith(line, "##INFO=<ID=ANN,")) {
                qpos <- as.integer(gregexpr("'", line, fixed = TRUE)[[1L]])
                if (length(qpos) < 2L || qpos[1L] < 0L) {
                    stop(sprintf(
                        "get_ann_fields: '%s' has ##INFO=<ID=ANN but the Description= block has no single-quoted field list. Cannot parse SnpEff schema.",
                        basename(path)
                    ))
                }
                inner <- substr(line, qpos[1L] + 1L, qpos[length(qpos)] - 1L)
                fields <- trimws(strsplit(inner, "|", fixed = TRUE)[[1L]])
                # Drop trailing empties from a stray trailing pipe
                fields <- fields[nzchar(fields)]
                return(fields)
            }
        }
        stop(sprintf("get_ann_fields: no ##INFO=<ID=ANN,...> header in '%s'.",
            basename(path)))
    }



    # -----------------------------------------------------------------------------
    # .snpeff_strip_allele(x)
    # Handles SnpEff's two non-standard allele forms (Cingolani spec):
    #   - cancer somatic vs germline: "G-C" -> "G"
    #   - compound variant:           "C-chr1:123456_A>T" -> "C"
    # Standard ALT (no dash) is returned unchanged. Empty input -> empty output.
    # Vectorized over `x`.
    # -----------------------------------------------------------------------------
    .snpeff_strip_allele <- function(x) {
        out <- vapply(strsplit(x, "-", fixed = TRUE),
            function(parts) if (length(parts)) parts[1L] else "",
            character(1L))
        out[!nzchar(x)] <- ""
        out
    }



    # -----------------------------------------------------------------------------
    # .CSQ_FIELDS_VEP_CANONICAL
    # Canonical 80-field VEP CSQ column-name vector emitted by read.gvr() (VEP).
    # read.gvr.snpeff() emits the same 80 names so downstream code is annotator-
    # agnostic. Only the subset that SnpEff ANN supplies is populated; the rest
    # are emitted as "" (per Round 1+2 user decision).
    #
    # Sourced from S1.vep_02.vcf.gz CSQ header (VEP v113.0, verified prior session).
    # -----------------------------------------------------------------------------
    .CSQ_FIELDS_VEP_CANONICAL <- c(
        "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type",
        "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp",
        "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
        "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS",
        "VARIANT_CLASS", "SYMBOL_SOURCE", "HGNC_ID", "CANONICAL", "MANE",
        "MANE_SELECT", "MANE_PLUS_CLINICAL", "TSL", "APPRIS", "CCDS", "ENSP",
        "SWISSPROT", "TREMBL", "UNIPARC", "UNIPROT_ISOFORM", "GENE_PHENO",
        "SIFT", "PolyPhen", "DOMAINS", "miRNA", "HGVS_OFFSET", "AF",
        "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF",
        "gnomADe_AF", "gnomADe_AFR_AF", "gnomADe_AMR_AF", "gnomADe_ASJ_AF",
        "gnomADe_EAS_AF", "gnomADe_FIN_AF", "gnomADe_MID_AF", "gnomADe_NFE_AF",
        "gnomADe_REMAINING_AF", "gnomADe_SAS_AF",
        "gnomADg_AF", "gnomADg_AFR_AF", "gnomADg_AMI_AF", "gnomADg_AMR_AF",
        "gnomADg_ASJ_AF", "gnomADg_EAS_AF", "gnomADg_FIN_AF", "gnomADg_MID_AF",
        "gnomADg_NFE_AF", "gnomADg_REMAINING_AF", "gnomADg_SAS_AF",
        "MAX_AF", "MAX_AF_POPS", "CLIN_SIG", "SOMATIC", "PHENO", "PUBMED",
        "MOTIF_NAME", "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE",
        "TRANSCRIPTION_FACTORS"
    )



    # -----------------------------------------------------------------------------
    # convert_chunk_snpeff(dt, ann_fields, sample_name, fns, opts)
    #
    # Per-chunk SnpEff -> table row decoder. Mirrors the data flow of VEP
    # convert_chunk() in read.gvr.R but parses SnpEff ANN blocks instead of CSQ.
    #
    # Arguments
    # ---------
    # dt           data.table with VCF columns: CHROM POS ID REF ALT QUAL FILTER
    #              INFO FORMAT SAMPLE. One row per VCF record (not per ALT, not
    #              per ANN block).
    # ann_fields   character vector of ANN field names, from get_ann_fields().
    # sample_name  string used to fill Tumor_Sample_Barcode.
    # fns          list of helper functions: info_parse, gvr_coords,
    #              vep_to_gvr_class, url_decode, make_hgvsp_short,
    #              gt_codes_for_alt. (closure-free design: passed in explicitly.)
    # opts         list of chunk-state variables:
    #                ncbi_build, dp_min, gq_min,
    #                do_dpgq_filter, do_gene_filter, do_vc_rough,
    #                genes_chr (or NULL),
    #                csq_fields_out (the 80-name vector this function emits)
    #
    # Returns a data.table with one row per (record, ALT) where any ANN block
    # matched the ALT (after .snpeff_strip_allele). Carries attr(., "n_dropped")
    # with the rough-filter drop count.
    # -----------------------------------------------------------------------------
    convert_chunk_snpeff <- function(dt, ann_fields, sample_name, fns, opts) {
        n_ann <- length(ann_fields)
        ai_ <- function(name) match(name, ann_fields)
        # ANN field positions (canonical 16 per Cingolani spec; absent ones -> NA)
        A_Allele       <- ai_("Allele")
        A_Annotation   <- ai_("Annotation")
        A_Impact       <- ai_("Annotation_Impact")
        A_Gene_Name    <- ai_("Gene_Name")
        A_Gene_ID      <- ai_("Gene_ID")
        A_Feature_Type <- ai_("Feature_Type")
        A_Feature_ID   <- ai_("Feature_ID")
        A_BioType      <- ai_("Transcript_BioType")
        A_Rank         <- ai_("Rank")
        A_HGVSc        <- ai_("HGVS.c")
        A_HGVSp        <- ai_("HGVS.p")
        A_cDNA         <- ai_("cDNA.pos / cDNA.length")
        A_CDS          <- ai_("CDS.pos / CDS.length")
        A_AA           <- ai_("AA.pos / AA.length")
        A_Distance     <- ai_("Distance")
        # ERRORS / WARNINGS / INFO field is at position 16 by convention but
        # we don't emit it; ignore.

        csq_fields_out <- opts$csq_fields_out
        n_csq_out <- length(csq_fields_out)

        # Position lookups in the 80-name canonical output vector for the few
        # ANN-derived slots. Cache them once so the inner loop is index-based.
        C_Allele     <- match("Allele",       csq_fields_out)
        C_Consequence <- match("Consequence",  csq_fields_out)
        C_IMPACT     <- match("IMPACT",       csq_fields_out)
        C_SYMBOL     <- match("SYMBOL",       csq_fields_out)
        C_Gene       <- match("Gene",         csq_fields_out)
        C_Feature_t  <- match("Feature_type", csq_fields_out)
        C_Feature    <- match("Feature",      csq_fields_out)
        C_BIOTYPE    <- match("BIOTYPE",      csq_fields_out)
        C_EXON       <- match("EXON",         csq_fields_out)
        C_HGVSc      <- match("HGVSc",        csq_fields_out)
        C_HGVSp      <- match("HGVSp",        csq_fields_out)
        C_cDNA       <- match("cDNA_position", csq_fields_out)
        C_CDS        <- match("CDS_position", csq_fields_out)
        C_Protein    <- match("Protein_position", csq_fields_out)
        C_DISTANCE   <- match("DISTANCE",     csq_fields_out)
        C_CLIN_SIG   <- match("CLIN_SIG",     csq_fields_out)
        C_gnomADe_AF <- match("gnomADe_AF",   csq_fields_out)

        # Hoist helpers from fns into local names for inner-loop speed
        info_parse        <- fns$info_parse
        gvr_coords        <- fns$gvr_coords
        vep_to_gvr_class  <- fns$vep_to_gvr_class
        url_decode        <- fns$url_decode
        make_hgvsp_short  <- fns$make_hgvsp_short
        gt_codes_for_alt  <- fns$gt_codes_for_alt

        ncbi_build      <- opts$ncbi_build
        dp_min          <- opts$dp_min
        gq_min          <- opts$gq_min
        do_dpgq_filter  <- isTRUE(opts$do_dpgq_filter)
        do_gene_filter  <- isTRUE(opts$do_gene_filter)
        do_vc_rough     <- isTRUE(opts$do_vc_rough)
        genes_chr       <- opts$genes_chr

        # C3: hoist columns to plain vectors once (no per-record `$` overhead)
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

        # O2: split ALT once per chunk
        alts_all <- strsplit(ALT_v, ",", fixed = TRUE)

        # O1: FORMAT-constant fast path -- pre-extract GT/AD/DP/GQ for whole chunk
        GT_col <- AD_col <- DP_col <- GQ_col <- NULL
        fmt_constant <- length(unique(FORMAT_v)) == 1L
        if (fmt_constant && nrow_dt > 0L) {
            fmt_keys <- strsplit(FORMAT_v[1L], ":", fixed = TRUE)[[1L]]
            p_GT <- match("GT", fmt_keys)
            p_AD <- match("AD", fmt_keys)
            p_DP <- match("DP", fmt_keys)
            p_GQ <- match("GQ", fmt_keys)
            splits <- data.table::tstrsplit(SAMPLE_v, ":", fixed = TRUE)
            if (!is.na(p_GT)) GT_col <- splits[[p_GT]]
            if (!is.na(p_AD)) AD_col <- splits[[p_AD]]
            if (!is.na(p_DP)) DP_col <- splits[[p_DP]]
            if (!is.na(p_GQ)) GQ_col <- splits[[p_GQ]]
        }

        # O6: DP/GQ pre-filter (vectorized, BEFORE main loop). Missing/non-numeric
        # values are treated as "pass" -- same convention as VEP convert_chunk.
        n_dropped_dpgq <- 0L
        if (do_dpgq_filter && nrow_dt > 0L) {
            if (fmt_constant) {
                dp_chr <- if (is.null(DP_col)) rep(NA_character_, nrow_dt) else DP_col
                gq_chr <- if (is.null(GQ_col)) rep(NA_character_, nrow_dt) else GQ_col
            } else {
                # Per-record: parse FORMAT for each record (slow path)
                dp_chr <- character(nrow_dt)
                gq_chr <- character(nrow_dt)
                for (r in seq_len(nrow_dt)) {
                    fk <- strsplit(FORMAT_v[r], ":", fixed = TRUE)[[1L]]
                    sv <- strsplit(SAMPLE_v[r], ":", fixed = TRUE)[[1L]]
                    names(sv) <- fk[seq_along(sv)]
                    dp_chr[r] <- if ("DP" %in% names(sv)) sv[["DP"]] else NA_character_
                    gq_chr[r] <- if ("GQ" %in% names(sv)) sv[["GQ"]] else NA_character_
                }
            }
            # why: as.integer() on the DP character field which may contain '' for missing; NAs are tolerated by the is.na() | dp_num>=dp_min filter.
            dp_num <- suppressWarnings(as.integer(dp_chr))
            # why: as.integer() on the GQ character field which may contain ''; NAs handled by the same is.na()|>=gq_min pattern.
            gq_num <- suppressWarnings(as.integer(gq_chr))
            pass_dp <- is.na(dp_num) | dp_num >= dp_min
            pass_gq <- is.na(gq_num) | gq_num >= gq_min
            pass <- pass_dp & pass_gq
            n_dropped_dpgq <- sum(!pass)
            if (n_dropped_dpgq > 0L) {
                idx <- which(pass)
                CHROM_v <- CHROM_v[idx]
                POS_v <- POS_v[idx]
                ID_v <- ID_v[idx]
                REF_v <- REF_v[idx]
                ALT_v <- ALT_v[idx]
                QUAL_v <- QUAL_v[idx]
                FILTER_v <- FILTER_v[idx]
                INFO_v <- INFO_v[idx]
                FORMAT_v <- FORMAT_v[idx]
                SAMPLE_v <- SAMPLE_v[idx]
                alts_all <- alts_all[idx]
                if (fmt_constant) {
                    GT_col <- if (!is.null(GT_col)) GT_col[idx] else NULL
                    AD_col <- if (!is.null(AD_col)) AD_col[idx] else NULL
                    DP_col <- if (!is.null(DP_col)) DP_col[idx] else NULL
                    GQ_col <- if (!is.null(GQ_col)) GQ_col[idx] else NULL
                }
                nrow_dt <- length(idx)
            }
        }

        # O7: gene rough filter via pipe-delimited match against INFO.
        # ANN blocks are pipe-separated and Gene_Name lives at position 4. Bare
        # substring grep would false-positive (e.g. "RET" matching "RETAIN"); the
        # |GENE_NAME| pattern guarantees a field-boundary match. This is more
        # permissive than VEP's because ANN field separators in INFO are simpler.
        n_dropped_gene_rough <- 0L
        if (do_gene_filter && !is.null(genes_chr) && length(genes_chr) > 0L && nrow_dt > 0L) {
            pat <- paste0(paste0("\\|", genes_chr, "\\|"), collapse = "|")
            gene_hit <- grepl(pat, INFO_v)
            n_dropped_gene_rough <- sum(!gene_hit)
            if (n_dropped_gene_rough > 0L) {
                idx <- which(gene_hit)
                CHROM_v <- CHROM_v[idx]
                POS_v <- POS_v[idx]
                ID_v <- ID_v[idx]
                REF_v <- REF_v[idx]
                ALT_v <- ALT_v[idx]
                QUAL_v <- QUAL_v[idx]
                FILTER_v <- FILTER_v[idx]
                INFO_v <- INFO_v[idx]
                FORMAT_v <- FORMAT_v[idx]
                SAMPLE_v <- SAMPLE_v[idx]
                alts_all <- alts_all[idx]
                if (fmt_constant) {
                    GT_col <- if (!is.null(GT_col)) GT_col[idx] else NULL
                    AD_col <- if (!is.null(AD_col)) AD_col[idx] else NULL
                    DP_col <- if (!is.null(DP_col)) DP_col[idx] else NULL
                    GQ_col <- if (!is.null(GQ_col)) GQ_col[idx] else NULL
                }
                nrow_dt <- length(idx)
            }
        }

        # O8: vc_nonSyn rough filter via SO terms. Identical term list to VEP path
        # because SnpEff's Annotation field uses the same Sequence Ontology
        # vocabulary (verified at pcingola.github.io). Compound terms can be joined
        # by '&' just as in VEP CSQ.
        n_dropped_vc_rough <- 0L
        if (do_vc_rough && nrow_dt > 0L) {
            nonSyn_so_terms <- c(
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
            pat_vc <- paste0(nonSyn_so_terms, collapse = "|")
            vc_hit <- grepl(pat_vc, INFO_v, ignore.case = FALSE)
            n_dropped_vc_rough <- sum(!vc_hit)
            if (n_dropped_vc_rough > 0L) {
                idx <- which(vc_hit)
                CHROM_v <- CHROM_v[idx]
                POS_v <- POS_v[idx]
                ID_v <- ID_v[idx]
                REF_v <- REF_v[idx]
                ALT_v <- ALT_v[idx]
                QUAL_v <- QUAL_v[idx]
                FILTER_v <- FILTER_v[idx]
                INFO_v <- INFO_v[idx]
                FORMAT_v <- FORMAT_v[idx]
                SAMPLE_v <- SAMPLE_v[idx]
                alts_all <- alts_all[idx]
                if (fmt_constant) {
                    GT_col <- if (!is.null(GT_col)) GT_col[idx] else NULL
                    AD_col <- if (!is.null(AD_col)) AD_col[idx] else NULL
                    DP_col <- if (!is.null(DP_col)) DP_col[idx] else NULL
                    GQ_col <- if (!is.null(GQ_col)) GQ_col[idx] else NULL
                }
                nrow_dt <- length(idx)
            }
        }

        out <- vector("list", nrow_dt)
        oi <- 0L
        n_dropped <- n_dropped_dpgq + n_dropped_gene_rough + n_dropped_vc_rough

        for (r in seq_len(nrow_dt)) {
            chrom <- CHROM_v[r]
            pos <- POS_v[r]
            vid <- ID_v[r]
            ref <- REF_v[r]
            altf <- ALT_v[r]
            qual <- QUAL_v[r]
            filt <- FILTER_v[r]
            info <- INFO_v[r]
            fmt <- FORMAT_v[r]
            smp <- SAMPLE_v[r]
            alts <- alts_all[[r]]

            ip <- info_parse(info)
            info_DP <- unname(ip["DP"])
            info_AC <- unname(ip["AC"])
            info_AF <- unname(ip["AF"])
            info_MQ <- unname(ip["MQ"])
            info_QD <- unname(ip["QD"])
            info_CNN <- unname(ip["CNN_1D"])
            ac_vec <- if (!is.na(info_AC)) strsplit(info_AC, ",", fixed = TRUE)[[1L]] else NA
            af_vec <- if (!is.na(info_AF)) strsplit(info_AF, ",", fixed = TRUE)[[1L]] else NA

            # ANN parsing: blocks separated by ',', fields by '|'.
            ann_raw <- unname(ip["ANN"])
            ann_blocks <- if (!is.na(ann_raw)) strsplit(ann_raw, ",", fixed = TRUE)[[1L]] else character(0L)
            # Pad each block to n_ann (base strsplit drops trailing empties)
            block_fields <- lapply(strsplit(ann_blocks, "|", fixed = TRUE), function(f) {
                length(f) <- n_ann
                f
            })
            # Per-block allele string, stripped of SnpEff's cancer/compound suffix
            block_allele_raw <- vapply(block_fields, function(f) {
                a <- f[A_Allele]
                if (is.na(a)) "" else a
            }, character(1L))
            block_allele <- .snpeff_strip_allele(block_allele_raw)

            # AF / CLIN_SIG INFO scan via canonical name list (first hit wins)
            info_gnomADe_AF <- ""
            for (nm in .AF_FIELDS_SNPEFF) {
                v <- ip[nm]
                if (!is.na(v) && nzchar(v)) {
                    info_gnomADe_AF <- unname(v)
                    break
                }
            }
            info_CLIN_SIG <- ""
            for (nm in .CLINSIG_FIELDS_SNPEFF) {
                v <- ip[nm]
                if (!is.na(v) && nzchar(v)) {
                    info_CLIN_SIG <- unname(v)
                    break
                }
            }

            # dbSNP_RS: SnpEff has no Existing_variation analog, so we look only in
            # the VCF ID column.
            dbsnp <- NA_character_
            if (!is.na(vid) && vid != ".") {
                rs <- grep("^rs", strsplit(vid, ";", fixed = TRUE)[[1L]], value = TRUE)
                if (length(rs)) dbsnp <- rs[1L]
            }

            # GT/AD/DP/GQ extraction -- identical logic to VEP path
            if (fmt_constant) {
                gt0  <- if (is.null(GT_col)) NA_character_ else GT_col[r]
                gt   <- if (is.na(gt0)) "./." else gt0
                ad   <- if (is.null(AD_col)) NA_character_ else AD_col[r]
                sdp  <- if (is.null(DP_col)) NA_character_ else DP_col[r]
                gq   <- if (is.null(GQ_col)) NA_character_ else GQ_col[r]
            } else {
                fmt_keys <- strsplit(fmt, ":", fixed = TRUE)[[1L]]
                smp_vals <- strsplit(smp, ":", fixed = TRUE)[[1L]]
                names(smp_vals) <- fmt_keys[seq_along(smp_vals)]
                gt  <- if ("GT" %in% names(smp_vals)) smp_vals[["GT"]] else "./."
                ad  <- if ("AD" %in% names(smp_vals)) smp_vals[["AD"]] else NA_character_
                sdp <- if ("DP" %in% names(smp_vals)) smp_vals[["DP"]] else NA_character_
                gq  <- if ("GQ" %in% names(smp_vals)) smp_vals[["GQ"]] else NA_character_
            }
            ad_vec <- if (!is.na(ad)) strsplit(ad, ",", fixed = TRUE)[[1L]] else NA

            # Per-ALT loop. For each ALT, pick the FIRST ANN block matching that ALT
            # (per Round 2 user decision: trust SnpEff's most-deleterious-first sort).
            for (ai in seq_along(alts)) {
                alt <- alts[ai]
                sel <- which(block_allele == alt)
                coords <- gvr_coords(pos, ref, alt)
                vt <- coords$var_type

                if (length(sel) > 0L) {
                    chosen <- block_fields[[sel[1L]]]
                } else {
                    chosen <- rep(NA_character_, n_ann)
                }

                # Derive ANN-sourced values
                ann_allele <- if (!is.na(A_Allele)) chosen[A_Allele]       else NA_character_
                ann_cons   <- if (!is.na(A_Annotation)) chosen[A_Annotation]   else NA_character_
                ann_impact <- if (!is.na(A_Impact)) chosen[A_Impact]       else NA_character_
                ann_gname  <- if (!is.na(A_Gene_Name)) chosen[A_Gene_Name]    else NA_character_
                ann_gid    <- if (!is.na(A_Gene_ID)) chosen[A_Gene_ID]      else NA_character_
                ann_ftype  <- if (!is.na(A_Feature_Type)) chosen[A_Feature_Type] else NA_character_
                ann_fid    <- if (!is.na(A_Feature_ID)) chosen[A_Feature_ID]   else NA_character_
                ann_bio    <- if (!is.na(A_BioType)) chosen[A_BioType]      else NA_character_
                ann_rank   <- if (!is.na(A_Rank)) chosen[A_Rank]         else NA_character_
                ann_hgvsc  <- if (!is.na(A_HGVSc)) chosen[A_HGVSc]        else NA_character_
                ann_hgvsp  <- if (!is.na(A_HGVSp)) chosen[A_HGVSp]        else NA_character_
                ann_cDNA   <- if (!is.na(A_cDNA)) chosen[A_cDNA]         else NA_character_
                ann_CDS    <- if (!is.na(A_CDS)) chosen[A_CDS]          else NA_character_
                ann_AA     <- if (!is.na(A_AA)) chosen[A_AA]           else NA_character_
                ann_dist   <- if (!is.na(A_Distance)) chosen[A_Distance]     else NA_character_

                # SnpEff Annotation field uses '&' to join compound consequences (same
                # convention as VEP). Pick the first term for Variant_Classification.
                top_term <- if (is.na(ann_cons) || ann_cons == "") NA_character_ else {
                    strsplit(ann_cons, "&", fixed = TRUE)[[1L]][1L]
                }
                var_class <- if (is.na(top_term)) "Targeted_Region" else vep_to_gvr_class(top_term, vt)
                hugo <- if (is.na(ann_gname) || ann_gname == "") "Unknown" else ann_gname

                gc <- gt_codes_for_alt(gt, ai)
                map_code <- function(code) {
                    if (is.na(code)) return(".")
                    if (code == 0L) return(coords$ref_allele)
                    if (code == ai) return(coords$tum_allele2)
                    if (code >= 1L && code <= length(alts)) return(gvr_coords(pos, ref, alts[code])$tum_allele2)
                    "."
                }
                t_allele1 <- map_code(gc$c1)
                t_allele2 <- map_code(gc$c2)

                t_ref_count <- if (length(ad_vec) >= 1 && !any(is.na(ad_vec))) ad_vec[1L] else NA_character_
                t_alt_count <- if (length(ad_vec) >= (ai + 1L) && !any(is.na(ad_vec))) ad_vec[ai + 1L] else NA_character_

                hgvsp_dec <- url_decode(ann_hgvsp)
                hgvsc_dec <- url_decode(ann_hgvsc)

                # Build the 80-name csq_fields_out vector for this record. Most slots
                # are "" (SnpEff doesn't provide them). The few that we DO have go in
                # their canonical positions.
                vep_vals <- as.list(rep("", n_csq_out))
                names(vep_vals) <- csq_fields_out
                # Populate the SnpEff-derivable slots
                if (!is.na(C_Allele)     && !is.na(ann_allele)) vep_vals[[C_Allele]]      <- ann_allele
                if (!is.na(C_Consequence) && !is.na(ann_cons)) vep_vals[[C_Consequence]] <- ann_cons
                if (!is.na(C_IMPACT)     && !is.na(ann_impact)) vep_vals[[C_IMPACT]]      <- ann_impact
                if (!is.na(C_SYMBOL)     && !is.na(ann_gname)) vep_vals[[C_SYMBOL]]      <- ann_gname
                if (!is.na(C_Gene)       && !is.na(ann_gid)) vep_vals[[C_Gene]]        <- ann_gid
                if (!is.na(C_Feature_t)  && !is.na(ann_ftype)) vep_vals[[C_Feature_t]]   <- ann_ftype
                if (!is.na(C_Feature)    && !is.na(ann_fid)) vep_vals[[C_Feature]]     <- ann_fid
                if (!is.na(C_BIOTYPE)    && !is.na(ann_bio)) vep_vals[[C_BIOTYPE]]     <- ann_bio
                if (!is.na(C_EXON)       && !is.na(ann_rank)) vep_vals[[C_EXON]]        <- ann_rank
                if (!is.na(C_HGVSc)      && !is.na(hgvsc_dec)) vep_vals[[C_HGVSc]]       <- hgvsc_dec
                if (!is.na(C_HGVSp)      && !is.na(hgvsp_dec)) vep_vals[[C_HGVSp]]       <- hgvsp_dec
                if (!is.na(C_cDNA)       && !is.na(ann_cDNA)) vep_vals[[C_cDNA]]        <- ann_cDNA
                if (!is.na(C_CDS)        && !is.na(ann_CDS)) vep_vals[[C_CDS]]         <- ann_CDS
                if (!is.na(C_Protein)    && !is.na(ann_AA)) vep_vals[[C_Protein]]     <- ann_AA
                if (!is.na(C_DISTANCE)   && !is.na(ann_dist)) vep_vals[[C_DISTANCE]]    <- ann_dist
                if (!is.na(C_CLIN_SIG)   && nzchar(info_CLIN_SIG)) vep_vals[[C_CLIN_SIG]]   <- info_CLIN_SIG
                if (!is.na(C_gnomADe_AF) && nzchar(info_gnomADe_AF)) vep_vals[[C_gnomADe_AF]] <- info_gnomADe_AF

                oi <- oi + 1L
                out[[oi]] <- c(
                    list(
                        Hugo_Symbol = hugo, Entrez_Gene_Id = "0", Center = ".", NCBI_Build = ncbi_build,
                        Chromosome = chrom, Start_Position = coords$start, End_Position = coords$end,
                        Strand = "+", Variant_Classification = var_class, Variant_Type = vt,
                        Reference_Allele = coords$ref_allele, Tumor_Seq_Allele1 = t_allele1,
                        Tumor_Seq_Allele2 = t_allele2, dbSNP_RS = if (is.na(dbsnp)) "" else dbsnp,
                        Tumor_Sample_Barcode = sample_name, Match_Norm_Seq_Allele1 = "",
                        Match_Norm_Seq_Allele2 = "",
                        HGVSc = if (is.na(hgvsc_dec)) "" else hgvsc_dec,
                        HGVSp = if (is.na(hgvsp_dec)) "" else hgvsp_dec,
                        HGVSp_Short = make_hgvsp_short(hgvsp_dec),
                        Transcript_ID = if (is.na(ann_fid)) "" else ann_fid,
                        Consequence = if (is.na(ann_cons)) "" else ann_cons,
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
        }
        res <- data.table::rbindlist(out[seq_len(oi)], use.names = TRUE, fill = TRUE)
        data.table::setattr(res, "n_dropped", n_dropped)
        res
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
            stop("read.gvr.snpeff: panel resolution requires gvr_panels.R; please source/load the package.")
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
    # vN+2 setup: HPO phenotype-term resolution (Phase N+2)
    #   Byte-identical to the read.gvr() HPO block. If `hpo` is supplied,
    #   resolve each HPO id to its gene vector via gvr_hpo_genes() and UNION
    #   with the existing `genes` (which may already include panel members
    #   from the vN block above). Runs AFTER panel so an HPO term augments,
    #   never replaces, a panel-supplied gene list. When `hpo = NULL` this
    #   block is a strict no-op and `genes` is byte-identical to what the
    #   vN block emitted.
    # ==========================================================================
    if (!is.null(hpo)) {
        if (!exists("gvr_hpo_genes", mode = "function", inherits = TRUE)) {
            stop("read.gvr.snpeff: HPO resolution requires gvr_hpo.R; please source/load the package.")
        }
        # Testing hook: options(gvr.hpo_path=...) allows hermetic test runs.
        .hpo_path_opt <- getOption("gvr.hpo_path", NULL)
        hpo_genes <- gvr_hpo_genes(hpo,
                                   hpo_path = .hpo_path_opt,
                                   verbose  = verbose)         # uppercased, deduped
        if (length(hpo_genes) > 0L) {
            effective_genes <- unique(toupper(trimws(c(hpo_genes, as.character(genes)))))
            effective_genes <- effective_genes[!is.na(effective_genes) & nzchar(effective_genes)]
            # Capture pre-union HPO state for the post-filter "hpo subset:"
            # coverage message (sibling of "panel subset:").
            .hpo_terms_for_summary   <- unique(trimws(as.character(hpo)))
            .hpo_genes_for_summary   <- hpo_genes
            .hpo_extras_for_summary  <- setdiff(
                toupper(trimws(as.character(genes))), hpo_genes
            )
            .hpo_extras_for_summary <- .hpo_extras_for_summary[
                !is.na(.hpo_extras_for_summary) & nzchar(.hpo_extras_for_summary)
            ]
            if (verbose) {
                message(sprintf(
                    "hpo: resolved %d HPO term(s) -> %d unique Hugo_Symbol(s)%s.",
                    length(.hpo_terms_for_summary),
                    length(effective_genes),
                    if (!is.null(genes) && length(as.character(genes)) > 0L)
                        sprintf(" (hpo %d + extras %d)",
                            length(hpo_genes),
                            length(setdiff(toupper(trimws(as.character(genes))), hpo_genes)))
                    else ""))
            }
            genes <- effective_genes
        } else if (verbose) {
            message("hpo: no genes resolved for supplied HPO term(s); ",
                    "no additional gene restriction applied.")
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
    #    Three input modes, precedence vcf_path > file > folder (vN+4):
    #      - vcf_path : character vector of full paths to process (multi-file OK).
    #      - file     : character vector of basenames resolved against `folder`.
    #      - folder   : list.files(folder, pattern) -- the default fan-out mode.
    #    `vcf_path` and `file` are mutually exclusive.
    # ==========================================================================
    if (!is.null(vcf_path) && !is.null(file))
        stop("read.gvr.snpeff: pass either `vcf_path=` (full paths) or `file=` ",
            "(basenames inside `folder=`), not both.", call. = FALSE)

    if (!is.null(vcf_path)) {
        vcf_path <- as.character(vcf_path)
        if (!length(vcf_path) || any(!nzchar(vcf_path)))
            stop("read.gvr.snpeff: `vcf_path=` must be a non-empty character vector.",
                call. = FALSE)
        miss <- vcf_path[!file.exists(vcf_path)]
        if (length(miss))
            stop(sprintf("read.gvr.snpeff: vcf_path file(s) do not exist:\n%s",
                paste0("  - ", miss, collapse = "\n")), call. = FALSE)
        vcf_paths <- vcf_path
    } else if (!is.null(file)) {
        file <- as.character(file)
        if (!length(file) || any(!nzchar(file)))
            stop("read.gvr.snpeff: `file=` must be a non-empty character vector of basenames.",
                call. = FALSE)
        if (!dir.exists(folder))
            stop(sprintf("Folder does not exist: %s", folder), call. = FALSE)
        candidates <- file.path(folder, file)
        miss <- file[!file.exists(candidates)]
        if (length(miss))
            stop(sprintf("read.gvr.snpeff: file(s) not found under folder='%s':\n%s",
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

    # vN+4: warn loudly once if canonical_only=TRUE was passed to snpeff -- ANN
    # has no CANONICAL flag, so the filter is silently no-op. Emit and continue.
    if (isTRUE(canonical_only)) {
        warning("read.gvr.snpeff: SnpEff ANN field has no CANONICAL flag; ",
            "canonical_only=TRUE has no effect. Result is identical to ",
            "canonical_only=FALSE.", call. = FALSE)
    }

    if (verbose) {
        message(sprintf("Found %d file(s):", length(vcf_paths)))
        for (p in vcf_paths) message("  - ", basename(p))
    }

    # ==========================================================================
    # 0b. Preflight annotator check (symmetric to read.gvr()).
    #     Three abort branches, in order:
    #       (i)   any file lacks both CSQ and ANN  -> un-annotated input
    #       (ii)  both VEP and SnpEff present     -> mixed annotators in batch
    #       (iii) all files are VEP-annotated     -> wrong route, use read.gvr()
    #     All ANN-only batches fall through unchanged into the chunk loop below.
    # ==========================================================================
    if (!exists(".detect_annotator", mode = "function")) {
        stop("read.gvr.snpeff: '.detect_annotator' helper missing; source ",
            "'read.gvr.snpeff.R' before calling.", call. = FALSE)
    }
    .detected_se <- vapply(vcf_paths, .detect_annotator, character(1L))
    ## (i) Un-annotated files: instructional message with INFO-tag listing.
    .na_mask_se <- is.na(.detected_se)
    if (any(.na_mask_se)) {
        .unann_first <- vcf_paths[.na_mask_se][1L]
        .tags_first  <- attr(.detect_annotator(.unann_first), "info_tags")
        .tag_str <- if (length(.tags_first))
            sprintf("INFO tags found in '%s': %s.\n",
                basename(.unann_first),
                paste(utils::head(.tags_first, 7L), collapse = ", "))
        else ""
        stop(sprintf(
            "read.gvr.snpeff: file(s) have no SnpEff ANN INFO tag:\n%s\n%s%s",
            paste(sprintf("  - %s", basename(vcf_paths[.na_mask_se])), collapse = "\n"),
            .tag_str,
            paste0("Annotate these VCFs with SnpEff (writes ##INFO=<ID=ANN>) ",
                "before passing to read.gvr.snpeff(), or use read.gvr() for ",
                "VEP-annotated input.")
        ), call. = FALSE)
    }
    ## (ii) Mixed batch (both 'vep' and 'snpeff' present): same wording as read.gvr().
    .uniq_se <- unique(.detected_se)
    if (length(.uniq_se) > 1L) {
        .by_ann <- split(basename(vcf_paths), .detected_se)
        .parts  <- vapply(names(.by_ann), function(a)
            sprintf("  %s: %s", a, paste(.by_ann[[a]], collapse = ", ")),
        character(1L))
        stop(sprintf(
            "read.gvr.snpeff: mixed annotators in batch -- refusing to merge:\n%s\n%s",
            paste(.parts, collapse = "\n"),
            "Use read.gvr() for VEP files, or process each annotator separately."
        ), call. = FALSE)
    }
    ## (iii) Pure VEP batch reached the SnpEff entry point: wrong route.
    if (.uniq_se == "vep") {
        stop(sprintf(
            "read.gvr.snpeff: file(s) are VEP-annotated, not SnpEff:\n%s\nUse read.gvr() for VEP files, or process each annotator separately.",
            paste(sprintf("  - %s", basename(vcf_paths)), collapse = "\n")
        ), call. = FALSE)
    }
    ## Fall-through: every file detects as 'snpeff'; continue to chunk loop.

    # ==========================================================================
    # 0c. Resolve ncbi_build: auto-detect from VCF header(s), honour user override.
    #     (Phase N+2; byte-identical logic to read.gvr())
    # ==========================================================================
    {
        .bd_first <- vcf_paths[1L]
        .bd <- tryCatch(.detect_build(.bd_first),
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
                    "read.gvr.snpeff: ncbi_build='%s' but VCF appears to be %s (%s agree). Continuing with user-supplied value.",
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
        ncbi_build <- effective_ncbi_build
    }

    # ==========================================================================
    # 1. Local helper definitions (kept inside the function => fully self-contained)
    #    [UNCHANGED CORE ENGINE from v1]
    # ==========================================================================
    ## 1b. VEP consequence -> Variant_Classification  (A2: memoized)
    ## var_class is a PURE function of (top_term, var_type); both are small-cardinality
    ## (19 distinct classes; a handful of var_types), so cache on the composite key.
    ## Inner function unchanged; cached value byte-identical.
    .v2m_cache <- new.env(parent = emptyenv())
    .vep_to_gvr_class_uncached <- function(term, var_type) {
        m <- switch(term,
            "splice_acceptor_variant" = "Splice_Site", "splice_donor_variant" = "Splice_Site",
            "transcript_ablation" = "Splice_Site", "exon_loss_variant" = "Splice_Site",
            "stop_gained" = "Nonsense_Mutation", "stop_lost" = "Nonstop_Mutation",
            "start_lost" = "Translation_Start_Site", "initiator_codon_variant" = "Translation_Start_Site",
            "missense_variant" = "Missense_Mutation", "conservative_missense_variant" = "Missense_Mutation",
            "rare_amino_acid_variant" = "Missense_Mutation", "protein_altering_variant" = "Missense_Mutation",
            "transcript_amplification" = "Intron", "splice_region_variant" = "Splice_Region",
            "splice_donor_5th_base_variant" = "Splice_Region", "splice_donor_region_variant" = "Splice_Region",
            "splice_polypyrimidine_tract_variant" = "Splice_Region", "stop_retained_variant" = "Silent",
            "synonymous_variant" = "Silent", "incomplete_terminal_codon_variant" = "Silent",
            "coding_sequence_variant" = "Missense_Mutation", "mature_miRNA_variant" = "RNA",
            "5_prime_UTR_variant" = "5'UTR", "5_prime_UTR_premature_start_codon_gain_variant" = "5'UTR",
            "3_prime_UTR_variant" = "3'UTR", "non_coding_transcript_exon_variant" = "RNA",
            "non_coding_exon_variant" = "RNA", "non_coding_transcript_variant" = "RNA",
            "nc_transcript_variant" = "RNA", "NMD_transcript_variant" = "Silent",
            "intron_variant" = "Intron", "upstream_gene_variant" = "5'Flank",
            "downstream_gene_variant" = "3'Flank", "TFBS_ablation" = "Targeted_Region",
            "TFBS_amplification" = "Targeted_Region", "TF_binding_site_variant" = "IGR",
            "regulatory_region_ablation" = "Targeted_Region", "regulatory_region_amplification" = "Targeted_Region",
            "regulatory_region_variant" = "IGR", "feature_elongation" = "Targeted_Region",
            "feature_truncation" = "Targeted_Region", "intergenic_variant" = "IGR",
            "Targeted_Region"
        )
        if (term == "frameshift_variant") m <- if (var_type == "DEL") "Frame_Shift_Del" else "Frame_Shift_Ins"
        if (term == "inframe_insertion") m <- "In_Frame_Ins"
        if (term == "inframe_deletion") m <- "In_Frame_Del"
        if (term == "disruptive_inframe_insertion") m <- "In_Frame_Ins"
        if (term == "disruptive_inframe_deletion") m <- "In_Frame_Del"
        m
    }
    vep_to_gvr_class <- function(term, var_type) {
        key <- paste0(if (is.na(term)) "\001NA" else term, "\002",
            if (is.na(var_type)) "\001NA" else var_type)
        hit <- .v2m_cache[[key]]
        if (!is.null(hit)) return(hit)
        val <- .vep_to_gvr_class_uncached(term, var_type)
        assign(key, val, envir = .v2m_cache)
        val
    }

    ## 1c. Amino-acid 3->1 and HGVSp_Short
    aa3to1 <- c(Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C", Gln = "Q", Glu = "E", Gly = "G",
        His = "H", Ile = "I", Leu = "L", Lys = "K", Met = "M", Phe = "F", Pro = "P", Ser = "S",
        Thr = "T", Trp = "W", Tyr = "Y", Val = "V", Ter = "*", Sec = "U", Pyl = "O", Asx = "B",
        Glx = "Z", Xaa = "X", Xle = "J")
    url_decode <- function(x) {
        if (is.na(x) || x == "") return(x)
        x <- gsub("%3D", "=", x, fixed = TRUE)
        x <- gsub("%3B", ";", x, fixed = TRUE)
        x <- gsub("%2C", ",", x, fixed = TRUE)
        x <- gsub("%3A", ":", x, fixed = TRUE)
        x
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
        s <- url_decode(hgvsp)
        s <- sub("^[^:]*:", "", s)
        s <- gsub("[()]", "", s)
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

    ## 1e. Coordinate + allele conversion for one REF/ALT pair
    # v0.99.2: Delegate to the shared .gvr_coords() helper (formerly a verbatim
    # inline copy). This ensures the SnpEff-only reader honors the
    # `normalize_alleles` argument (bcftools-norm-style prefix/suffix trimming
    # of REF/ALT before deriving the trimmed coords), which prevents distinct multi-ALT
    # records from collapsing to the same (chrom, start, ref, alt) join key.
    # Closes over `normalize_alleles` from read.gvr.snpeff()'s signature.
    gvr_coords <- function(pos, ref, alt) {
        .gvr_coords(pos, ref, alt, normalize_alleles = normalize_alleles)
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
        # why: as.integer() on a split genotype string; '.' and other non-numeric tokens become NA and are dropped by !is.na() below.
        gidx <- suppressWarnings(as.integer(strsplit(gt, "[/|]")[[1]]))
        gidx <- gidx[!is.na(gidx)]
        if (length(gidx) == 0L) return(list(c1 = NA_integer_, c2 = ai))
        others <- gidx[gidx != ai]
        if (length(others) == 0L) return(list(c1 = ai, c2 = ai))   # hom-alt
        partner <- if (0L %in% others) 0L else others[1]
        list(c1 = partner, c2 = ai)
    }

    ## 1h. Header parsers (SnpEff path: get_sample_name verbatim; get_ann_fields
    ##     is defined at top level in read.gvr.snpeff.R).
    get_sample_name <- function(path) {
        con <- gzfile(path, "r")
        on.exit(close(con))
        repeat {
            line <- readLines(con, n = 1L)
            if (length(line) == 0L) stop("#CHROM not found")
            if (startsWith(line, "#CHROM")) {
                cols <- strsplit(line, "\t", fixed = TRUE)[[1]]
                return(cols[length(cols)])
            }
        }
    }

    ## 1j. Convert ONE snpeff-annotated vcf file end-to-end. Returns a per-file variant data.table
    ##     data.table. Mirrors the VEP convert_one_vcf in read.gvr.R; the only
    ##     differences are (a) get_ann_fields instead of get_csq_fields, (b)
    ##     convert_chunk_snpeff instead of convert_chunk, and (c) the verbose
    ##     message wording ("ANN fields" vs "CSQ fields"). Progress wiring,
    ##     parallel safety, and exception handling are identical.
    convert_one_vcf <- function(path, file_idx, file_total, file_total_records = NA_integer_) {
        ann_fields  <- get_ann_fields(path)
        sample_name <- get_sample_name(path)
        if (verbose)
            message(sprintf("[file %d/%d] Converting %s (sample: %s) | ANN fields: %d",
                file_idx, file_total, basename(path), sample_name, length(ann_fields)))

        # Pack helpers + chunk-state for the closure-free convert_chunk_snpeff()
        fns_pack <- list(
            info_parse        = info_parse,
            gvr_coords        = gvr_coords,
            vep_to_gvr_class  = vep_to_gvr_class,
            url_decode        = url_decode,
            make_hgvsp_short  = make_hgvsp_short,
            gt_codes_for_alt  = gt_codes_for_alt
        )
        opts_pack <- list(
            ncbi_build      = ncbi_build,
            dp_min          = if (filter_dp) min_DP else 0L,
            gq_min          = if (filter_gq) min_GQ else 0L,
            do_dpgq_filter  = filter_dp || filter_gq,
            do_gene_filter  = !is.null(genes) && length(genes) > 0L,
            do_vc_rough     = isTRUE(vc_nonSyn),
            genes_chr       = genes,
            csq_fields_out  = .CSQ_FIELDS_VEP_CANONICAL
        )

        con <- gzfile(path, "r")
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
        t0 <- Sys.time()
        file_done <- 0L
        repeat {
            lines <- readLines(con, chunk_size)
            if (length(lines) == 0L) break
            ci2 <- ci2 + 1L
            mat <- data.table::tstrsplit(lines, "\t", fixed = TRUE)
            dtc <- data.table::data.table(CHROM = mat[[1]], POS = mat[[2]], ID = mat[[3]], REF = mat[[4]],
                ALT = mat[[5]], QUAL = mat[[6]], FILTER = mat[[7]], INFO = mat[[8]],
                FORMAT = mat[[9]], SAMPLE = mat[[10]])
            ck <- convert_chunk_snpeff(dtc, ann_fields, sample_name, fns_pack, opts_pack)
            nd <- attr(ck, "n_dropped")
            if (is.null(nd)) nd <- 0L
            n_drop_file <- n_drop_file + nd
            chunks[[ci2]] <- ck
            total_in <- total_in + nrow(dtc)
            file_done <- file_done + nrow(dtc)
            .progress$global_done <- .progress$global_done + nrow(dtc)
            if (verbose) {
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
        attr(gvr_one, "sample_name") <- sample_name
        attr(gvr_one, "n_dropped")   <- n_drop_file
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
            # why: suppressMessages on the worker call so per-chunk progress lines from convert_one_vcf() do not interleave with the parent heartbeat output.
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
        # Vectorized URL-decode (fast; mirrors the per-element url_decode helper).
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
            # the core-columns copy writes "" for missing while the raw CSQ copy writes NA,
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

    ## 3d-bis3. Phase N+2: post-filter HPO coverage message -------------------
    ##     Sibling of "panel subset:" (3d-bis2). Fires only when `hpo` was
    ##     supplied AND resolved to >=1 gene (gated on
    ##     exists(.hpo_genes_for_summary), which was set inside the vN+2
    ##     setup block at the top of the body). Reports how many of the
    ##     ORIGINAL (pre-union) HPO-derived genes are present in the
    ##     post-filter result, names the requested terms, lists present +
    ##     missing, and (when `hpo` was combined with explicit `genes`
    ##     extras) names them inline. Observation point matches "gene
    ##     subset:" semantics (pre vc_nonSyn / drop_empty_cols).
    if (verbose && exists(".hpo_genes_for_summary", inherits = FALSE)) {
        data.table::setDT(gvr)
        have_syms <- unique(toupper(trimws(as.character(gvr$Hugo_Symbol))))
        h_present <- intersect(.hpo_genes_for_summary, have_syms)
        h_missing <- setdiff(.hpo_genes_for_summary, have_syms)
        extras_chunk <- if (length(.hpo_extras_for_summary) > 0L) {
            sprintf(" (+ %d extra gene(s) from `genes`: %s)",
                length(.hpo_extras_for_summary),
                paste(.hpo_extras_for_summary, collapse = ", "))
        } else ""
        missing_chunk <- if (length(h_missing) > 0L) {
            sprintf(" (missing: %s)", paste(h_missing, collapse = ", "))
        } else ""
        message(sprintf(
            "hpo subset: %d / %d HPO gene(s) present for %s: %s%s%s.",
            length(h_present), length(.hpo_genes_for_summary),
            paste(.hpo_terms_for_summary, collapse = ", "),
            if (length(h_present) > 0L) paste(h_present, collapse = ", ") else "(none)",
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
    # (R's saveRDS preserves attributes; CSV/TSV is just columns so the attribute
    #  is for in-memory consumers -- but we still keep it consistent.)
    data.table::setattr(gvr, "annotator", "snpeff")

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
