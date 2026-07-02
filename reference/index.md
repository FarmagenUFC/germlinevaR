# Package index

## Readers

Read VEP-, SnpEff-, or dual-annotated VCFs into an MAF-like data.table.

- [`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
  : Convert VEP-annotated germline VCF(s) to a tabular variant
  data.table
- [`read.gvr.snpeff()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.snpeff.md)
  : Convert SnpEff-annotated germline VCF(s) to a tabular variant
  data.table
- [`read.gvr.dual()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.dual.md)
  : Convert dual-annotated germline VCF(s) (VEP + SnpEff) to a tabular
  variant data.table

## Filtering

Modular filtering and novel-variant subsetting.

- [`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md)
  : Modular, individually-toggleable filtering of a read.gvr table
- [`gvr_novel()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_novel.md)
  : Filter a read.gvr table down to candidate novel variants

## Gene panels

- [`gvr_panel_genes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_panel_genes.md)
  : Genes for a Disease Panel
- [`gvr_list_panels()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_list_panels.md)
  : List Available Disease Gene Panels
- [`gvr_hpo_genes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_hpo_genes.md)
  : Resolve HPO phenotype terms to associated genes

## Cohort summary

- [`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md)
  : Multi-section summary of a germline gvr table (read.gvr / gvr_filter
  output)
- [`gvr_sum_plots()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_sum_plots.md)
  : Export gvr_summary() plots as standalone image files

## Per-gene plots

- [`gvr_plot()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_plot.md)
  : Cohort top-genes variant matrix from a germline gvr table (read.gvr
  / gvr_filter output)
- [`gvr_lollipop()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_lollipop.md)
  : Per-gene amino-acid lollipop plot for a germline gvr table
- [`gvr_genepos.plot()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_genepos.plot.md)
  : Gene-track lollipop plot on a cDNA axis

## Palette and cache utilities

- [`gvr_color_palette()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_color_palette.md)
  : Generate a Vector of Colors From a Named Palette

- [`gvr_list_palettes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_list_palettes.md)
  : List All Available Palettes

- [`gvr_domain_cache_clear()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_domain_cache_clear.md)
  :

  Clear the auto-fetched protein-domain cache used by `gvr_lollipop`
