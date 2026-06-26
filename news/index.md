# Changelog

## germlinevaR 0.99.0

### New features

- Initial Bioconductor submission.
- [`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md):
  auto-routing entry point for VEP-, SnpEff-, and dual-annotated
  germline VCFs. Inspects each VCF header and delegates to the
  appropriate sibling reader automatically.
- [`read.gvr.snpeff()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.snpeff.md):
  dedicated reader for SnpEff-annotated VCFs (INFO `ANN` field).
- [`read.gvr.dual()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.dual.md):
  reader for VCFs carrying both VEP `CSQ` and SnpEff `ANN` annotations
  on the same records.
- [`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md):
  modular filtering by population allele frequency (gnomAD exomes, 1000
  Genomes, ABraOM), ClinVar clinical significance, biotype, variant
  classification, genotype quality, and gene panel.
- [`gvr_novel()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_novel.md):
  dedicated subsetter for candidate novel variants (no dbSNP rsID and no
  allele frequency in any catalogue).
- [`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md):
  8-section cohort summary returning a list of `data.table`s, with
  optional Excel workbook, multi-page PDF, and self-contained
  interactive HTML dashboard (plotly + DT).
- [`gvr_plot()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_plot.md):
  `ComplexHeatmap`-based top-genes variant matrix across samples,
  coloured by variant classification.
- [`gvr_lollipop()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_lollipop.md):
  per-gene protein-domain lollipop plot with on-the-fly InterPro domain
  fetching (disk-cached) and a 52-palette colour system.
- [`gvr_genepos.plot()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_genepos.plot.md):
  per-gene gene-structure (cDNA) lollipop plot with exon/intron/UTR
  regions resolved via Ensembl REST or a local GTF file (disk-cached).
- [`gvr_sum_plots()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_sum_plots.md):
  standalone PNG/SVG/PDF exports of all dashboard panels produced by
  [`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md).
- [`gvr_color_palette()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_color_palette.md),
  [`gvr_list_palettes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_list_palettes.md):
  unified 52-palette colour system with colourblind-safe defaults shared
  across every figure.
- [`gvr_panel_genes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_panel_genes.md),
  [`gvr_list_panels()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_list_panels.md):
  curated disease gene-panel registry (breast cancer, Lynch syndrome,
  GIST, and others).
- [`gvr_domain_cache_clear()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_domain_cache_clear.md):
  cache management for InterPro domain data.
