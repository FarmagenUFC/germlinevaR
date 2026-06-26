# Modular, individually-toggleable filtering of a read.gvr table

Applies a set of independent variant filters to the MAF-like table
produced by
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md).
Each distinct filter is its own argument; setting an argument to `NULL`
disables that filter entirely (no rows removed by it). With all
defaults, `gvr_filter(gvr)` reproduces the canonical rare /
clinically-relevant / called-genotype pipeline (AF filters + CLIN_SIG +
GT exclusion).

## Usage

``` r
gvr_filter(
  gvr,
  gnomADe_AF = 0.01,
  AF = 0.01,
  ABraOM_AF = 0.01,
  gnomADe_AF_keep_missing = TRUE,
  AF_keep_missing = TRUE,
  ABraOM_AF_keep_missing = TRUE,
  clin_sig_terms = c("likely_pathogenic", "pathogenic", "uncertain_significance"),
  clin_sig_keep_missing = TRUE,
  remove_benign = FALSE,
  biotype_keep = NULL,
  gt_exclude = c("0", "0/0"),
  vc_nonSyn = FALSE,
  missense_only = FALSE,
  genes = NULL,
  save_excel = FALSE,
  out_dir = NULL,
  file_prefix = "gvr_filter",
  verbose = TRUE
)
```

## Arguments

- gvr:

  data.table / data.frame from
  [`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
  (or compatible). Filtered in a copy; the input object is not modified.

- gnomADe_AF:

  Numeric upper threshold for the gnomAD exome AF column `gnomADe_AF`
  (keep rows with AF \< threshold). `NULL` disables this filter. Default
  0.01.

- AF:

  Numeric upper threshold for the `AF` column (gnomAD genome / VCF
  allele frequency). `NULL` disables this filter. Default 0.01.

- ABraOM_AF:

  Numeric upper threshold for the Brazilian-cohort (ABraOM SABE 609)
  column `ABraOM_AF`. `NULL` disables this filter. Default 0.01.

- gnomADe_AF_keep_missing:

  Logical; if TRUE (default), keep rows whose `gnomADe_AF` is missing
  (NA or ""); if FALSE drop them. Ignored when `gnomADe_AF` is NULL.

- AF_keep_missing:

  Logical; missing-value handling for the `AF` filter. TRUE (default)
  keeps missing. Ignored when `AF` is NULL.

- ABraOM_AF_keep_missing:

  Logical; missing-value handling for the ABraOM filter. TRUE (default)
  retains variants absent from the Brazilian cohort (where absence often
  means "not catalogued", not "common"). Ignored when `ABraOM_AF` is
  NULL.

- clin_sig_terms:

  Character vector of clinical-significance terms to keep (substring,
  case-insensitive, OR-combined). `NULL` disables the CLIN_SIG filter.
  Default: c("likely_pathogenic","pathogenic","uncertain_significance").

- clin_sig_keep_missing:

  Logical; if TRUE (default) rows with missing CLIN_SIG (NA/"") are
  kept. Only relevant when `clin_sig_terms` is non-NULL.

- remove_benign:

  Logical; if TRUE, remove rows whose `CLIN_SIG` contains "benign"
  (substring, case-insensitive). This catches `benign`, `likely_benign`,
  and compound annotations like
  `"uncertain_significance&likely_benign"`. Applied AFTER the
  `clin_sig_terms` keep-filter, so a row that matched a wanted term but
  also contains "benign" is still removed. `FALSE` (default) does not
  remove benign rows.

- biotype_keep:

  Character vector of BIOTYPE values to keep (exact match via %in%).
  `NULL` (default) disables the biotype filter — all biotypes are kept.
  Pass e.g. `c("protein_coding", "protein_coding_LoF")` to restrict to
  protein-coding transcripts.

- gt_exclude:

  Character vector of GT values to remove (exact match). `NULL` disables
  the genotype filter. Default: c("0","0/0").

- vc_nonSyn:

  Logical or character vector. Controls which `Variant_Classification`
  values are retained. `FALSE` (default) keeps all. `TRUE` keeps only
  the 9 protein-altering classes (Frame_Shift_Del, Frame_Shift_Ins,
  Splice_Site, Translation_Start_Site, Nonsense_Mutation,
  Nonstop_Mutation, In_Frame_Del, In_Frame_Ins, Missense_Mutation). A
  custom character vector keeps only those classifications. Rows with
  missing/blank `Variant_Classification` are removed when this filter is
  active.

- missense_only:

  Logical; if `TRUE`, keep only rows whose `Variant_Classification`
  equals `"Missense_Mutation"` (added in vN+5). Default `FALSE`
  preserves prior behaviour byte-for-byte. Combines non-contradictorily
  with `vc_nonSyn`: `vc_nonSyn` runs first (keeping 9 protein-altering
  classes), then `missense_only` narrows to the missense subset. Errors
  with a clear message if `Variant_Classification` is missing.

- genes:

  Character vector of `Hugo_Symbol`s to keep (exact, case-insensitive),
  or `NULL` (default) to keep all genes.

- save_excel:

  Logical; if TRUE, also write the FILTERED table to an `.xlsx` workbook
  (single `"Filtered"` sheet) at `<out_dir>/<file_prefix>.xlsx`.
  Requires the openxlsx package (a `Suggests` dependency); if it is not
  installed the export is skipped with a warning. Default FALSE. The
  write is a side effect only: the returned `data.table` is identical
  whether or not `save_excel` is TRUE.

- out_dir:

  Output directory for the Excel file. `NULL` (default) uses the current
  working directory. Created if it does not exist. Only used when
  `save_excel = TRUE`.

- file_prefix:

  Filename prefix (without extension) for the Excel file. Default
  `"gvr_filter"` -\> `gvr_filter.xlsx`. Only used when
  `save_excel = TRUE`.

- verbose:

  Logical; if TRUE (default) print a per-filter breakdown (rows in -\>
  out and rows removed by each active step).

## Value

A `data.table` of the surviving rows, with the same columns as the
input. A plain `data.frame` input is returned as a `data.table`. The
input object is not modified. With `verbose = TRUE`, a per-filter
breakdown (rows in -\> out, and rows removed by each active step) is
printed as it runs.

## Details

Filters are applied in a fixed order; each step operates on the
survivors of the previous one:

1.  gnomAD exome AF - `gnomADe_AF` (+ `gnomADe_AF_keep_missing`)

2.  gnomAD genome / VCF AF - `AF` (+ `AF_keep_missing`)

3.  ABraOM AF - `ABraOM_AF` (+ `ABraOM_AF_keep_missing`)

4.  Clinical significance - `clin_sig_terms` (+ `clin_sig_keep_missing`)

5.  Remove benign - `remove_benign`

6.  Biotype - `biotype_keep`

7.  Genotype exclusion - `gt_exclude`

8.  Variant classification - `vc_nonSyn`

9.  Gene subset - `genes`

Important data notes (true of
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
output):

- AF columns are CHARACTER (e.g. `"0.8781"`), so they are coerced with
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html) before
  comparison.

- "Missing" means EITHER `NA` OR empty string `""`
  ([`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
  uses `""` for absent values). Both are treated as missing everywhere
  in this function.

- `CLIN_SIG` matching is SUBSTRING + case-insensitive. A compound
  annotation such as `"pathogenic&benign"` or
  `"uncertain_significance&likely_benign&benign"` is KEPT because it
  CONTAINS a wanted term. This matches the dplyr `str_detect()` /
  [`grepl()`](https://rdrr.io/r/base/grep.html) convention. Use
  exact-token matching only if you split `CLIN_SIG` yourself.

- The default `gt_exclude = c("0", "0/0")` is a no-op on data whose `GT`
  column only contains called alt genotypes (e.g. `0/1`, `1/1`, `1/2`);
  it is retained for portability to data that does carry `"0"`/`"0/0"`.

The input is never modified: filtering operates on an internal
`data.table` copy.

## See also

[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
to build the table,
[`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md)
to summarise the filtered variants.

Other germlinevaR:
[`gvr_list_panels()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_list_panels.md),
[`gvr_panel_genes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_panel_genes.md),
[`gvr_plot()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_plot.md),
[`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md),
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md),
[`read.gvr.snpeff()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.snpeff.md)

## Author

germlinevaR authors

## Examples

``` r
## Load the shipped example table
gvr <- readRDS(system.file("extdata", "example_gvr.rds",
                           package = "germlinevaR"))
## Default filter (rare + clinically relevant + called genotypes)
filt <- gvr_filter(gvr, verbose = FALSE)
#> Warning: gvr_filter: column 'ABraOM_AF' not found in data; skipping ABraOM_AF filter. (Set ABraOM_AF = NULL to silence this warning.)
dim(filt)
#> [1]   7 115

if (FALSE) { # \dontrun{
gvr <- read.gvr("/path/to/vcf_folder")

## Default pipeline: rare variants + clinically relevant + called genotypes:
gvr_clean <- gvr_filter(gvr)

## Add protein-coding biotype filter:
gvr_filter(gvr, biotype_keep = c("protein_coding", "protein_coding_LoF"))

## Only the rarity filter on gnomAD exome AF, nothing else:
gvr_filter(gvr, gnomADe_AF = 0.001, AF = NULL, ABraOM_AF = NULL,
           clin_sig_terms = NULL, gt_exclude = NULL,
           vc_nonSyn = FALSE, genes = NULL)

## Pathogenic-only, protein-coding:
gvr_filter(gvr, clin_sig_terms = c("pathogenic", "likely_pathogenic"),
           biotype_keep = "protein_coding")

## Remove benign annotations (including likely_benign and compound entries):
gvr_filter(gvr, remove_benign = TRUE)

## Keep only protein-altering variants and a gene panel:
gvr_filter(gvr, vc_nonSyn = TRUE, genes = c("TP53", "BRCA1", "BRCA2"))
} # }
```
