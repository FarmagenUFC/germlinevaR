# Convert SnpEff-annotated germline VCF(s) to an MAF-like data.table

Converts SnpEff-annotated, single-sample germline VCFs (GATK
HaplotypeCaller -\> CNN tranches -\> SnpEff, hg38) into an MAF-like
table and returns it as an in-memory `data.table` for downstream
filtering
([`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md))
and summarisation
([`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md)).
In folder mode it finds every per-sample VCF, converts each, and
row-binds them into one combined gvr table. Emits the SAME canonical
80-field schema as
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
so downstream code is annotator- agnostic. The conversion uses base R +
data.table only - no external annotation-package dependency.

This function is the SnpEff sibling of
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
(VEP-annotated VCFs). Both readers are usually invoked through
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md),
which auto-routes SnpEff inputs to this function via the nested
`.detect_annotator()` helper. Calling `read.gvr.snpeff()` directly is
supported for fully SnpEff-only pipelines.

## Usage

``` r
read.gvr.snpeff(
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
  convert. Mutually exclusive with `file=`. `NULL` (default) selects
  folder mode.

- file:

  Character vector of basenames (e.g.
  `c("S1.snpeff.vcf.gz", "S2.snpeff.vcf.gz")`) resolved against
  `folder=`. Use this to pick specific files from a folder that contains
  files you do NOT want to merge. Mutually exclusive with `vcf_path=`.
  `NULL` (default) selects either `vcf_path=` mode or folder-pattern
  mode.

- pattern:

  Regular expression identifying per-sample VCFs in folder mode. Default
  `"\\.vcf\\.gz$"` (matches any `.vcf.gz` file). The old default
  `"_\\d+(\\.(vep|snp[eE]ff))?\\.vcf\\.gz$"` is still available by
  passing it explicitly.

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
  is skipped with a warning. Default `FALSE`. Excel caps a sheet at
  1,048,576 rows.

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

  Logical; accepted for API symmetry with
  [`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md).
  SnpEff's ANN field has no `CANONICAL` flag, so this filter cannot be
  applied. When `TRUE` (the default), `read.gvr.snpeff()` emits a
  one-time warning and returns the unfiltered table (same result as
  `canonical_only = FALSE`).

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
columns, the canonical 80 VEP CSQ field names (populated from SnpEff ANN
where available; blank otherwise), key GATK QC fields,
`Tumor_Sample_Barcode`, and (when enabled) the `Genotype` and
`ABraOM_AF` columns. `attr(., "annotator") = "snpeff"`. TSV/RDS/XLSX
files are written as a side effect when
`write_tsv`/`write_rds`/`write_xlsx` is `TRUE`.

## Details

Output and behaviour:

- Returns the final MAF-like `data.table`, one row per variant ALLELE
  (multi-allelic sites are split). Multiple `ANN` annotation blocks per
  allele are reduced to one most-severe transcript per allele.

- Columns include the MAF-like core fields, the canonical 80 VEP CSQ
  field names (populated from the equivalent SnpEff ANN fields where
  available; blank for fields SnpEff does not supply), and key GATK QC
  fields. `FILTER` is retained as a column and ALL variants (PASS and
  non-PASS) are kept.

- `Tumor_Seq_Allele1`/`Tumor_Seq_Allele2` are zygosity-aware (vcf2maf-
  style), and an optional `Genotype` column
  (`Tumor_Seq_Allele1/Tumor_Seq_Allele2`, e.g. `"T/C"`) is added next to
  the alleles.

- Each variant keeps its source sample in `Tumor_Sample_Barcode`.

- Absent values are written as the empty string `""` (not `NA`);
  downstream
  [`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md)
  /
  [`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md)
  treat `NA` and `""` identically as "missing".

- Output is tagged with `attr(gvr, "annotator") = "snpeff"` so
  downstream code can distinguish the source.

SnpEff-specific notes:

- The nested `.detect_annotator()` helper inspects `##INFO=<ID=*` header
  lines and returns `"snpeff"` when `##INFO=<ID=ANN` is found (and
  `"vep"` when `##INFO=<ID=CSQ` is found). If both are present, VEP
  takes priority.

- The nested `get_ann_fields()` helper parses SnpEff's ANN INFO header
  line, extracting the `'F1 | F2 | ... | FN'` field-name vector from
  inside the single-quoted `Description=` block.

- The nested `.snpeff_strip_allele()` helper handles SnpEff's two
  non-standard allele forms (Cingolani spec): cancer-somatic-vs-germline
  (`"G-C"` -\> `"G"`) and compound (`"C-chr1:123456_A>T"` -\> `"C"`).
  Standard ALT alleles (no dash) pass through unchanged.

- Population AF columns (`gnomADe_AF`, etc.) are populated from the
  first matching SnpEff ANN-side INFO key in `.AF_FIELDS_SNPEFF`
  (case-sensitive, first match wins). Bare `AF` is intentionally
  excluded - in the VCF spec it is the per-ALT caller frequency, not a
  population frequency.

Processing options (same as
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)):

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

[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
for the VEP-annotated sibling reader (and the recommended entrypoint,
which auto-routes SnpEff input to this function),
[`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md)
to filter the returned table,
[`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md)
to summarise it.

Other germlinevaR:
[`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md),
[`gvr_list_panels()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_list_panels.md),
[`gvr_panel_genes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_panel_genes.md),
[`gvr_plot()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_plot.md),
[`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md),
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)

## Author

germlinevaR authors

## Examples

``` r
## The shipped example is VEP-annotated; read.gvr.snpeff() is shown here
## with a minimal demonstration using the auto-router on the VEP example.
## For SnpEff VCFs, call read.gvr.snpeff() the same way as read.gvr().
gvr_list_panels()   # confirm package is loaded
#>  [1] "acromegaly"                          
#>  [2] "breast cancer"                       
#>  [3] "breast cancer somatic"               
#>  [4] "familial adenomatous polyposis"      
#>  [5] "gist"                                
#>  [6] "hereditary breast and ovarian cancer"
#>  [7] "hereditary cancer"                   
#>  [8] "hereditary colorectal cancer"        
#>  [9] "hereditary gastric cancer"           
#> [10] "hereditary melanoma cancer"          
#> [11] "hereditary prostate cancer"          
#> [12] "li-fraumeni syndrome"                
#> [13] "lynch syndrome"                      
#> [14] "men1"                                
#> [15] "pheochromocytoma"                    

if (FALSE) { # \dontrun{
## Folder mode: merge ALL *_NN.snpeff.vcf.gz into one gvr table
gvr <- read.gvr.snpeff("/path/to/folder")

## Or use read.gvr() to auto-route (SnpEff inputs are dispatched here):
gvr <- read.gvr("/path/to/folder")

## Single-file mode: full path
gvr <- read.gvr.snpeff(vcf_path = "/path/to/SAMPLE_01.snpeff.vcf.gz")

## Multi-file mode by full path
gvr <- read.gvr.snpeff(
  vcf_path = c("/p/S1.snpeff.vcf.gz", "/p/S2.snpeff.vcf.gz"))

## Pick basenames from a folder (merges these two but ignores other .vcf.gz)
gvr <- read.gvr.snpeff(folder = "/p",
                       file   = c("S1.snpeff.vcf.gz", "S2.snpeff.vcf.gz"))

## Disable the DP/GQ genotype filter entirely
gvr <- read.gvr.snpeff("/path/to/folder", min_DP = NULL, min_GQ = NULL)

## Restrict to genes of interest
gvr <- read.gvr.snpeff("/path/to/folder",
                       genes = c("MEN1", "RET", "CDKN1B", "CDC73"))

## Or use a curated disease panel:
gvr_list_panels()
gvr <- read.gvr.snpeff("/path/to/folder", panel = "breast cancer")

## Multiple panels are unioned (deduplicated):
gvr <- read.gvr.snpeff("/path/to/folder",
                       panel = c("breast cancer", "hereditary prostate cancer"))
} # }
```
