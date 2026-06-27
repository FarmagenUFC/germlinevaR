# Convert VEP-annotated germline VCF(s) to an MAF-like data.table

Converts VEP-annotated, single-sample germline VCFs (GATK
HaplotypeCaller -\> CNN tranches -\> Ensembl VEP, hg38) into an MAF-like
table and returns it as an in-memory `data.table` for downstream
filtering
([`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md))
and summarisation
([`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md)).
In folder mode it finds every per-sample VCF, converts each, and
row-binds them into one combined gvr table. The conversion uses base R
and data.table only - no external annotation-package dependency. This is
the recommended entry point for all germline VCFs: `read.gvr()` inspects
each input's INFO tags and, when needed, delegates to
[`read.gvr.snpeff()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.snpeff.md)
for SnpEff-annotated VCFs or
[`read.gvr.dual()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.dual.md)
for VCFs carrying both VEP and SnpEff annotations, so a single call
handles every annotator combination.

## Usage

``` r
read.gvr(
  folder = ".",
  vcf_path = NULL,
  file = NULL,
  pattern = "\\.vcf\\.gz$",
  write_tsv = FALSE,
  write_rds = FALSE,
  write_xlsx = FALSE,
  out_dir = NULL,
  out_prefix = NULL,
  chunk_size = 25000L,
  ncbi_build = "auto",
  add_genotype = TRUE,
  strip_hgvs_prefix = TRUE,
  dedup_columns = TRUE,
  drop_empty_cols = FALSE,
  add_abraom = TRUE,
  abraom_path = NULL,
  abraom_url =
    "https://abraom.ib.usp.br/download/ABRaOM_60+_SABE_609_exomes_annotated.gz",
  cache_dir = NULL,
  min_DP = 10,
  min_GQ = 30,
  genes = NULL,
  panel = NULL,
  vc_nonSyn = FALSE,
  canonical_only = TRUE,
  ncores = 1L,
  verbose = TRUE
)
```

## Arguments

- folder:

  Directory to scan in folder mode; every file matching `pattern` is
  converted and row-bound. Default `"."`. Ignored when `vcf_path` is
  supplied. Also used as the search root for `file=`.

- vcf_path:

  Character vector of one or more full paths to `.vcf.gz` files to
  convert. Use this to process a specific set of files outside the
  folder pattern. Mutually exclusive with `file=`. `NULL` (default)
  selects folder mode.

- file:

  Character vector of basenames (e.g.
  `c("S1.vep.vcf.gz", "S2.vep.vcf.gz")`) resolved against `folder=`. Use
  this to pick specific files from a folder that contains files you do
  NOT want to merge. Mutually exclusive with `vcf_path=`. `NULL`
  (default) selects either `vcf_path=` mode or folder-pattern mode.

- pattern:

  Regular expression identifying per-sample VCFs in folder mode. Default
  `"\\.vcf\\.gz$"` (matches any `.vcf.gz` file). The old default
  `"_\\d+(\\.(vep|snp[eE]ff))?\\.vcf\\.gz$"` (requires `_NN` suffix) is
  still available by passing it explicitly.

- write_tsv:

  Logical; if `TRUE`, also write the table as a TSV to `out_dir`.
  Default `FALSE`.

- write_rds:

  Logical; if `TRUE`, also write the table as an `.rds` to `out_dir`.
  Default `FALSE`.

- write_xlsx:

  Logical; if `TRUE`, also write the table as an `.xlsx` workbook
  (single `"gvr_table"` sheet) to `out_dir`. Requires the openxlsx
  package (a `Suggests` dependency); if it is not installed the export
  is skipped with a warning. Default `FALSE`. Note: germline MAFs can be
  very large; Excel handles them but the file is big and slow to open,
  so `write_rds`/`write_tsv` are better for large tables. Excel also
  caps a sheet at 1,048,576 rows.

- out_dir:

  Output directory for written TSV/RDS/XLSX. `NULL` (default) uses the
  input location/working directory. Only used when
  `write_tsv`/`write_rds`/`write_xlsx` is `TRUE`.

- out_prefix:

  Filename prefix for written outputs. `NULL` (default) derives one from
  the input.

- chunk_size:

  Integer; number of VCF records processed per chunk (controls peak
  memory and progress granularity). Default `25000L`.

- ncbi_build:

  Reference build label written into the table `NCBI_Build` column.
  Default `"auto"`: the function inspects the input VCF header (VEP
  `assembly=` / SnpEff `SnpEffCmd` database token / first `##contig=`
  length) and picks the canonical label among `"GRCh38"`, `"GRCh37"`,
  `"T2T-CHM13v2.0"`. When detection cannot decide, falls back to
  `"GRCh38"` with a verbose-mode message. Pass any literal (`"GRCh37"`,
  `"hg19"`, `"T2T-CHM13v2.0"`, internal lab codes) to override; the
  user-supplied value is written verbatim into `NCBI_Build` and the rest
  of the pipeline does not branch on its value. The ABraOM join uses
  dbSNP rsID + alleles and is build-stable. The aliases `"hg19"` and
  `"hg38"` are mapped to `"GRCh37"` / `"GRCh38"` for the mismatch check;
  a diagnostic warning fires when an explicit `ncbi_build` value
  (canonical or alias) disagrees with what auto-detection found at high
  confidence. Off-table user labels (e.g. internal lab codes) pass
  through silently.

- add_genotype:

  Logical; if `TRUE` (default) add the `Genotype` column.

- strip_hgvs_prefix:

  Logical; if `TRUE` (default) strip the Ensembl feature prefix from
  `HGVSc`/`HGVSp`.

- dedup_columns:

  Logical; if `TRUE` (default) drop duplicate-named columns when
  byte-for-byte identical (otherwise keep and warn).

- drop_empty_cols:

  Logical; if `TRUE`, drop columns that are entirely `NA`/blank. Default
  `FALSE`.

- add_abraom:

  Logical; if `TRUE` (default) join the ABraOM SABE-609 allele frequency
  as `ABraOM_AF`.

- abraom_path:

  Path to a local ABraOM annotation file. `NULL` (default) uses an
  auto-managed cache (see `cache_dir`), downloading from `abraom_url` if
  needed.

- abraom_url:

  URL of the ABraOM SABE-609 annotated release used when `abraom_path`
  is `NULL`.

- cache_dir:

  Directory for the ABraOM reference cache (used only when `abraom_path`
  is `NULL`). `NULL` (default) uses
  `tools::R_user_dir("germlinevaR", "cache")`, which resolves to a
  platform-appropriate directory (e.g. `~/.cache/R/germlinevaR` on
  Linux). The directory is created on first download. Set explicitly if
  you prefer a custom location.

- min_DP:

  Numeric; keep only records with `DP > min_DP`. `NULL` (or `NA`)
  disables the depth filter. Default `10`.

- min_GQ:

  Numeric; keep only records with `GQ > min_GQ`. `NULL` (or `NA`)
  disables the genotype-quality filter. Default `30`.

- genes:

  Character vector of `Hugo_Symbol`s to keep (exact, case-insensitive),
  or `NULL` (default) to keep all genes.

- panel:

  Character vector of curated disease panel name(s) (e.g.
  `"breast cancer"`, `"hereditary prostate cancer"`, `"gist"`). Each
  name is resolved to a gene vector via
  [`gvr_panel_genes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_panel_genes.md)
  and the union of all resolved genes is taken with `genes`
  (deduplicated, uppercased). Names are matched case-insensitively,
  trimmed, and `_` is treated as a space, so `"Breast_Cancer"`,
  `"breast cancer"`, and `" BREAST CANCER "` all resolve identically. A
  small alias table is also recognised (e.g.
  `"gastrointestinal stromal tumor"` -\> `"gist"`, the typo
  `"pheocromocytoma"` -\> `"pheochromocytoma"`). The registry currently
  ships 15 panels; see
  [`gvr_list_panels()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_list_panels.md)
  for the full list. An unknown name raises an error listing the
  available panels. `NULL` (default) disables panel filtering; behaviour
  is then byte-identical to omitting the argument.

- vc_nonSyn:

  Logical or character vector. Controls which `Variant_Classification`
  values are retained (mirroring the convention of the `vc_nonSyn`
  argument). `FALSE` (default) keeps ALL variant classifications. `TRUE`
  keeps only protein-altering classes (High/Moderate VEP consequences):
  `"Frame_Shift_Del"`, `"Frame_Shift_Ins"`, `"Splice_Site"`,
  `"Translation_Start_Site"`, `"Nonsense_Mutation"`,
  `"Nonstop_Mutation"`, `"In_Frame_Del"`, `"In_Frame_Ins"`,
  `"Missense_Mutation"`. Alternatively, pass a custom character vector
  of classifications to keep. Rows with missing/blank
  `Variant_Classification` are always removed when this filter is
  active.

- canonical_only:

  Logical; when `TRUE` (default), drops table rows whose chosen VEP CSQ
  block has `CANONICAL != "YES"`. read.gvr() already prefers
  `CANONICAL=YES` when ranking CSQ blocks; this filter discards the
  fallback rows emitted when no canonical block exists for a given ALT.
  For SnpEff-annotated input, the ANN field has no CANONICAL flag —
  `canonical_only=TRUE` is ignored with a warning, and the result is the
  same as `canonical_only=FALSE`.

- ncores:

  Integer; number of worker processes for converting MULTIPLE input
  files in parallel via
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  (fork-based; Unix/macOS only). Default `1L` runs sequentially and is
  byte-identical to previous behaviour. Values `> 1` only help when more
  than one VCF is being read (each file is an independent task) and are
  clamped to `min(ncores, detectCores(), n_files)`. On non-fork
  platforms it falls back to sequential. A single file is unaffected.

- verbose:

  Logical; if `TRUE` (default) print per-file and per-chunk progress
  (file i/N, cumulative records, elapsed seconds).

## Value

An MAF-like `data.table`: one row per variant allele, with MAF-like core
columns, all VEP CSQ fields, key GATK QC fields, `Tumor_Sample_Barcode`,
and (when enabled) the `Genotype` and `ABraOM_AF` columns. TSV/RDS files
are written as a side effect when `write_tsv`/`write_rds` is `TRUE`.

## Details

Output and behaviour:

- Returns the final MAF-like `data.table`, one row per variant ALLELE
  (multi-allelic sites are split).

- A single most-severe transcript is chosen per allele (VEP severity -\>
  `CANONICAL` -\> `MANE_SELECT` -\> transcript id).

- Columns include the MAF-like core fields, ALL VEP CSQ fields (read
  from the VCF header), and key GATK QC fields. `FILTER` is retained as
  a column and ALL variants (PASS and non-PASS) are kept.

- `Tumor_Seq_Allele1`/`Tumor_Seq_Allele2` are zygosity-aware
  (vcf2maf-style), and an optional `Genotype` column
  (`Tumor_Seq_Allele1/Tumor_Seq_Allele2`, e.g. `"T/C"`) is added next to
  the alleles.

- Each variant keeps its source sample in `Tumor_Sample_Barcode`.

- Absent values are written as the empty string `""` (not `NA`);
  downstream
  [`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md)
  /
  [`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md)
  treat `NA` and `""` identically as "missing".

Processing options:

- MULTI-FILE: in folder mode, every file matching `pattern` (default
  `"*_NN.vcf.gz"`) is converted and row-bound.

- HGVS CLEANUP (`strip_hgvs_prefix`): strips the Ensembl feature prefix
  from `HGVSc`/`HGVSp` (e.g. `"ENST00000831140.1:n.1889G>A"` -\>
  `"n.1889G>A"`).

- DEDUP (`dedup_columns`): removes duplicate-named columns ONLY when
  their values are byte-for-byte identical across all rows (otherwise
  keeps + warns).

- ABraOM (`add_abraom`): joins the Brazilian ABraOM SABE-609 allele
  frequency as the `ABraOM_AF` column (downloaded/cached from
  `abraom_url`).

- GENOTYPE-QUALITY FILTER (`min_DP`/`min_GQ`): keeps a record iff
  `DP > min_DP` AND `GQ > min_GQ`; mirrors
  `bcftools view -e 'FORMAT/DP<=X | FORMAT/GQ<=Y'`. Set either to `NULL`
  to disable that field; set both to `NULL` to disable the genotype
  filter entirely.

- GENE SUBSET (`genes`): restrict to a set of `Hugo_Symbol`s (exact,
  case-insensitive).

## See also

[`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md)
to filter the returned table,
[`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md)
to summarise it,
[`read.gvr.snpeff()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.snpeff.md)
for SnpEff-annotated VCFs,
[`read.gvr.dual()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.dual.md)
for VCFs with both VEP and SnpEff annotations.

Other germlinevaR:
[`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md),
[`gvr_list_panels()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_list_panels.md),
[`gvr_panel_genes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_panel_genes.md),
[`gvr_plot()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_plot.md),
[`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md),
[`read.gvr.snpeff()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.snpeff.md)

## Author

germlinevaR authors

## Examples

``` r
## read.gvr() reads VEP/SnpEff-annotated VCFs and returns a parsed
## data.table. Reading the bundled 62-variant fixture takes ~20s on
## a typical CI worker, so the real read.gvr() calls are wrapped in
## a donttest block below; a pre-parsed equivalent .rds is also bundled.

## The function is exported and callable:
is.function(read.gvr)
#> [1] TRUE

## Pre-parsed equivalent (instantaneous; same shape and content as the
## result of read.gvr() on the same fixture, minus the ABraOM_AF column
## which was added in a later code revision):
gvr <- readRDS(system.file("extdata", "example_gvr.rds",
                           package = "germlinevaR"))
dim(gvr)
#> [1]  62 115

# \donttest{
  ## Real read.gvr() call on the bundled VCF directory:
  vcf_dir <- system.file("extdata", package = "germlinevaR")
  gvr <- read.gvr(vcf_dir, verbose = FALSE)
#> Warning: read.gvr: ABraOM reference unreadable; 'ABraOM_AF' left blank.
  dim(gvr)            # 62 rows x 116 columns
#> [1]  62 116

  ## Same call but with write-out to tempdir
  out <- tempdir()
  gvr2 <- read.gvr(vcf_dir, write_tsv = TRUE, write_rds = TRUE,
                   out_dir = out, verbose = FALSE)
#> Warning: read.gvr: ABraOM reference unreadable; 'ABraOM_AF' left blank.
  list.files(out, pattern = "\\.(tsv|rds)$")
#> [1] "example.vep.gvr.tsv.rds" "example.vep.gvr.tsv.tsv"

  ## Single-file mode: full path to the bundled VCF
  vcf_file <- system.file("extdata", "example.vep.vcf.gz",
                          package = "germlinevaR")
  gvr3 <- read.gvr(vcf_path = vcf_file, verbose = FALSE)
#> Warning: read.gvr: ABraOM reference unreadable; 'ABraOM_AF' left blank.
  nrow(gvr3)
#> [1] 62
# }
```
