# Multi-section summary of a germline gvr table (read.gvr / gvr_filter output)

Produces a multi-section overview of a MAF-like table - the output of
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md),
or of
[`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md) -
covering variant burden, affected genes, functional classes, clinical
significance and predicted impact. Every section is returned as a tidy
`data.table` with one column per sample plus a `Total` column.
Optionally writes a multi-sheet Excel workbook and/or a multi-page PDF
report (both into a `gvr_summary/` subfolder of `out_dir`). For a cohort
top-genes variant matrix, see
[`gvr_plot()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_plot.md).

## Usage

``` r
gvr_summary(
  gvr,
  sample_col = "Tumor_Sample_Barcode",
  top_n_genes = 20,
  top_n_variants = 20,
  save_excel = TRUE,
  save_pdf = TRUE,
  save_html = TRUE,
  out_dir = ".",
  file_prefix = "gvr_summary",
  verbose = TRUE
)
```

## Arguments

- gvr:

  An MAF-like `data.table`/`data.frame` from
  [`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
  or
  [`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md).
  Required columns: `Hugo_Symbol`, `Variant_Classification`,
  `Variant_Type`, `IMPACT`, `CLIN_SIG`.

- sample_col:

  Name of the per-sample column. Default `"Tumor_Sample_Barcode"`. If
  absent, all rows are pooled into a single sample `"All"` (with a
  warning).

- top_n_genes:

  Integer; number of genes to report in `top_genes` (by total variant
  count). Default `20`.

- top_n_variants:

  Integer; number of variants to report in `top_variants` (by `dbSNP_RS`
  frequency). Default `20`. Ignored if the `dbSNP_RS` column is absent
  from the input.

- save_excel:

  Logical; if `TRUE` (default), write a multi-sheet `.xlsx`. The
  workbook is written into the `gvr_summary/` subfolder of `out_dir`
  (see `out_dir`). Pass `FALSE` for a compute-only run.

- save_pdf:

  Logical; if `TRUE` (default), write a multi-page PDF dashboard report
  into the `gvr_summary/` subfolder of `out_dir`. Page 1 is a hero page
  (a row of KPI cards above two grouped/faceted bar charts - top genes
  and variant classification); the following pages hold the section
  tables (packed two-per-row where they fit, else full-width) plus the
  functional-impact chart. The layout adapts to cohort size (faceting
  and column pagination for many samples; see the examples). Requires
  gridExtra, ggplot2 and scales; if unavailable, the PDF is skipped with
  a warning and the sections are still returned. Pass `FALSE` for a
  compute-only run.

- save_html:

  Logical; if `TRUE` (default), write an interactive HTML dashboard
  (`<file_prefix>_report.html`) into the `gvr_summary/` subfolder of
  `out_dir`. It mirrors the PDF dashboard - a row of KPI cards, bar
  charts (top genes, variant classification, functional impact, top
  variants) as interactive plotly charts (grouped for \\\le 6\\ samples,
  faceted small-multiples for \\\> 6\\), and all section tables as
  sortable, searchable DT tables. The Clinical significance table is
  interactive: clicking a CLIN_SIG token (e.g. "pathogenic") expands a
  detail panel showing the individual variants with that annotation. By
  default a single self-contained file is produced (assets inlined via
  pandoc); if pandoc is unavailable the report is written as
  `<file_prefix>_report.html` plus a sibling
  `<file_prefix>_report_files/` asset folder (a `verbose` message notes
  this). Requires plotly, DT, htmlwidgets and htmltools; if any are
  unavailable the HTML is skipped with a warning and the sections are
  still returned. Pass `FALSE` for a compute-only run.

- out_dir:

  Parent output directory. All written outputs (Excel, PDF and/or HTML)
  are placed in a `gvr_summary/` subfolder created inside `out_dir`. The
  subfolder is created only when `save_excel`, `save_pdf` or `save_html`
  is `TRUE`. Default `"."` (current working directory), i.e. outputs go
  to `./gvr_summary/`.

- file_prefix:

  Base filename for written outputs. Default `"gvr_summary"`, giving
  `<file_prefix>.xlsx`, `<file_prefix>_report.pdf` and
  `<file_prefix>_report.html` (no timestamp). Filenames are fixed, so
  re-running into the same `out_dir` overwrites the previous files (a
  message is printed when `verbose = TRUE`).

- verbose:

  Logical; if `TRUE` (default) print a compact console digest and the
  path(s) of any file(s) written.

## Value

Invisibly, a named list of `data.table`s: `overview`, `top_genes`,
`top_genes_per_sample`, `variant_classification`, `variant_type`,
`clin_sig`, `top_variants`, `impact`. The `top_genes_per_sample` element
is itself a named list (one data.table per sample). The `top_variants`
section is absent if `dbSNP_RS` is not in the input. The return value is
identical regardless of whether the Excel/PDF/HTML files are written.

## Details

Sections returned (as a named list of `data.table`s):

- `overview` - cohort-level counts: total variants, variants per sample,
  distinct genes affected, and variants with no gene symbol.

- `top_genes` - the `top_n_genes` genes with the most variants
  (per-sample + `Total`); unknown/blank genes excluded.

- `variant_classification` - counts per `Variant_Classification`
  (functional class), per-sample + `Total`, sorted by `Total`
  descending.

- `variant_type` - counts per `Variant_Type` (SNP/DEL/INS/ONP/DNP/TNP).

- `clin_sig` - counts per `CLIN_SIG` token (ClinVar categories).
  `CLIN_SIG` is split on `&` and `/`, so a variant annotated
  `"pathogenic&benign"` increments BOTH categories; category counts
  therefore sum to \\\ge\\ the number of variants. A
  `missing/unclassified` row counts NA/"" `CLIN_SIG`.

- `top_genes_per_sample` - a named list of `data.table`s, one per
  sample, each containing the `top_n_genes` genes with the most variants
  in that sample (columns: `Hugo_Symbol`, `<sample_name>`).
  Unknown/blank genes excluded.

- `top_variants` - the `top_n_variants` most frequent variants by
  `dbSNP_RS` (rsID), with per-sample counts + `Total`. Rows with
  blank/missing `dbSNP_RS` are excluded. If the `dbSNP_RS` column is
  absent, this section is skipped with a warning.

- `impact` - counts per VEP `IMPACT` (HIGH/MODERATE/LOW/MODIFIER), in
  severity order rather than count order.

The section tables are the core return value. By default the function
also writes an Excel workbook (`save_excel = TRUE`) and a PDF report
(`save_pdf = TRUE`) into `out_dir/gvr_summary/`; set either to `FALSE`
to skip it. The cohort top-genes variant matrix lives in
[`gvr_plot()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_plot.md).

Data conventions:

- "Missing" means `NA` OR empty string `""`.

- Unknown/blank gene symbols are `Hugo_Symbol` in
  `c(".", "", "Unknown")`; these are excluded from the distinct-gene
  tally and from `top_genes`, but their variants are still counted in
  the totals (and reported as "variants with no gene symbol").

- Works on ANY MAF-like table; it makes no assumption about prior
  filtering. It is commonly run on
  [`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md)
  output to summarise the retained hits.

## Dependencies

Core summary uses data.table. The optional Excel export uses openxlsx;
the optional PDF dashboard uses gridExtra + ggplot2 + scales, rendered
via the [`grDevices::cairo_pdf`](https://rdrr.io/r/grDevices/cairo.html)
device (full Unicode, so en-dashes, multiplication signs and similar
punctuation render correctly). The optional interactive HTML dashboard
uses plotly (charts), DT (tables) and htmlwidgets + htmltools
(assembly); a single self-contained file is produced when pandoc is
available (used only for the optional asset-inlining step), otherwise a
`.html` + `_files/` asset folder is written. Each optional output
degrades gracefully: if its package(s) are unavailable, that output is
skipped with a warning and the section tables are still returned.

## See also

[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
to build the table,
[`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md)
to filter it before summarising,
[`gvr_plot()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_plot.md)
for a cohort top-genes variant matrix.

Other germlinevaR:
[`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md),
[`gvr_list_panels()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_list_panels.md),
[`gvr_panel_genes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_panel_genes.md),
[`gvr_plot()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_plot.md),
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md),
[`read.gvr.snpeff()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.snpeff.md)

## Author

germlinevaR authors

## Examples

``` r
## Load the shipped example table and run summary (no file output)
gvr <- readRDS(system.file("extdata", "example_gvr.rds",
    package = "germlinevaR"))
summ <- gvr_summary(gvr, save_excel = FALSE, save_pdf = FALSE,
    save_html = FALSE, verbose = FALSE)
names(summ)
#> [1] "overview"               "top_genes"              "top_genes_per_sample"  
#> [4] "variant_classification" "variant_type"           "clin_sig"              
#> [7] "impact"                 "top_variants"          

## Write the XLSX workbook + multi-page PDF + interactive HTML
## dashboard to a temp directory.
out_dir <- file.path(tempdir(), "gvr_summary")
s <- gvr_summary(gvr, out_dir = out_dir,
    save_excel = TRUE, save_pdf = TRUE, save_html = TRUE,
    verbose = FALSE)
## Inspect a section table
head(s$variant_classification)
#>    Variant_Classification Sample_01 Total
#>                    <char>     <int> <num>
#> 1:                 Intron        14    14
#> 2:                 Silent         7     7
#> 3:      Missense_Mutation         5     5
#> 4:                    RNA         4     4
#> 5:            Splice_Site         4     4
#> 6:                  3'UTR         3     3
```
