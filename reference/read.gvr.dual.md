# Convert dual-annotated germline VCF(s) (VEP + SnpEff) to an MAF-like data.table

Converts single-sample germline VCFs that carry **both** Ensembl VEP
`CSQ` and SnpEff `ANN` annotation INFO fields on the same records
(typical workflow: SnpEff -\> VEP, or vice versa). Returns an MAF-like
`data.table` with one row per ALT allele, using VEP's transcript pick as
the spine and adding SnpEff-derived comparison columns and LoF/NMD
predictions.

This function is the dual-annotator sibling of
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
(VEP-only) and
[`read.gvr.snpeff()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.snpeff.md)
(SnpEff-only). All three are normally invoked through
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md),
which auto-routes inputs based on the INFO tags found in the VCF header.

## Usage

``` r
read.gvr.dual(
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
  hpo = NULL,
  vc_nonSyn = FALSE,
  canonical_only = TRUE,
  ncores = 1L,
  normalize_alleles = TRUE,
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

- hpo:

  Character vector of Human Phenotype Ontology (HPO) term identifier(s),
  e.g. `"HP:0003002"` (Breast carcinoma) or
  `c("HP:0003002", "HP:0025022")`. Each term is resolved to a gene
  vector via
  [`gvr_hpo_genes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_hpo_genes.md)
  and the union of all resolved genes is added to any genes already
  produced by `genes` and `panel` (deduplicated, uppercased). Lenient
  input is accepted: `"HP:0003002"`, `"hp:0003002"`, `"3002"`, and
  `"0003002"` all normalise to the canonical `"HP:0003002"` form. Only
  exact-term associations are used (no ontology-descendant expansion).
  The HPO phenotype-to-gene table is downloaded and cached under
  `tools::R_user_dir("germlinevaR", "cache")` and auto-refreshed after
  30 days; see
  [`gvr_hpo_genes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_hpo_genes.md)
  for offline / air-gapped usage via `hpo_path=`. `NULL` (default)
  disables HPO filtering; behaviour is then byte-identical to omitting
  the argument.

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

- normalize_alleles:

  Logical; if `TRUE` (default, since 0.99.2) apply bcftools-norm-style
  trimming of common REF/ALT prefix and suffix nucleotides before
  deriving MAF-like coords (`Start_Position`, `Reference_Allele`,
  `Tumor_Seq_Allele2`). This is the recommended behaviour: it puts each
  variant on its unique minimal representation and prevents distinct
  multi-ALT records from collapsing to the same MAF key (which could
  previously drop or scramble annotations on the SnpEff side of the dual
  reader; VEP-only reads were unaffected). Set `FALSE` to reproduce
  pre-0.99.2 coords for reproducibility with an older analysis; note
  this is not recommended for new research.

- verbose:

  Logical; if `TRUE` (default) print per-file and per-chunk progress
  (file i/N, cumulative records, elapsed seconds).

## Value

A `data.table` with one row per ALT allele, ~124 columns. Carries
`attr(., "annotator") = "dual"`.

## Details

Field-level priority: **VEP wins.** The canonical table columns and all
80 (or 81 with `FREQS`) VEP-style CSQ columns hold VEP's pick
exclusively. SnpEff data is added in additional columns:

- `snpeff_consequence`, `snpeff_impact`, `snpeff_gene`, `snpeff_hgvsc` –
  the SnpEff ANN block matching the same ALT allele and gene as VEP's
  pick (most-deleterious block when multiple match). Empty when SnpEff
  has no ANN block for that allele.

- `LOF_Gene`, `LOF_Pct_Transcripts` – SnpEff's loss-of-function
  prediction (gene flagged, fraction of transcripts affected) from the
  INFO `LOF=` field. Empty when no LoF call.

- `NMD_Gene`, `NMD_Pct_Transcripts` – SnpEff's nonsense-mediated- decay
  prediction from the INFO `NMD=` field. Empty when no NMD call.

Transcript pick: VEP drives the row spine (one row per ALT allele). If
VEP was run with `--per_gene` (one transcript per gene per allele), the
chosen transcript is VEP's most-severe across genes for that allele;
otherwise the standard VEP most-severe-block ranking applies. SnpEff
annotations are attached *to* that row from the ANN block matching the
same allele and (preferentially) the same `Gene_Name`.

Schema additions vs
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
(VEP-only):

- `FREQS` – 81st canonical CSQ field, populated when VEP was run with
  `--everything` (v113+). Names the gnomAD population where MAX_AF was
  observed (e.g. "gnomADg_ASJ"). Empty when the CSQ header does not
  include this field.

- 4 LOF/NMD columns (above)

- 4 snpeff\_\* parallel columns (above)

Total: ~124 columns (vs ~116 for
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
and
[`read.gvr.snpeff()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.snpeff.md)).

All other behaviour, defaults, arguments, and filtering options are
identical to
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md);
see that function's documentation for the full option reference.

## See also

[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md),
[`read.gvr.snpeff()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.snpeff.md)

## Examples

``` r
## The function signature is exported and callable:
is.function(read.gvr.dual)
#> [1] TRUE

# \donttest{
## read.gvr.dual() expects VCFs annotated with BOTH VEP (CSQ INFO field)
## AND SnpEff (ANN INFO field) in the same record. The shipped example
## VCF is VEP-only, so a real dual-annotated example needs your own
## VCFs; we therefore guard the call on the directory existing, so the
## example skips cleanly on machines without the data.
dual_dir <- "/path/to/dual-annotated-vcfs/"
if (dir.exists(dual_dir)) {
    gvr <- read.gvr.dual(folder = dual_dir)

    ## Or via the auto-router in read.gvr() when the VCF header declares
    ## both VEP and SnpEff INFO fields:
    gvr <- read.gvr(dual_dir)

    ## Compare VEP vs SnpEff picks on high-impact variants:
    gvr[IMPACT == "HIGH" & snpeff_impact != "" & IMPACT != snpeff_impact,
        .(Hugo_Symbol, Consequence, IMPACT, snpeff_gene, snpeff_consequence,
            snpeff_impact)]
}
# }
```
