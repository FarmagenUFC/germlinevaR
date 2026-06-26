# germlinevaR 0.99.0

## New features

- Initial Bioconductor submission.
- `read.gvr()`: auto-routing entry point for VEP-, SnpEff-, and dual-annotated
  germline VCFs. Inspects each VCF header and delegates to the appropriate
  sibling reader automatically.
- `read.gvr.snpeff()`: dedicated reader for SnpEff-annotated VCFs (INFO `ANN`
  field).
- `read.gvr.dual()`: reader for VCFs carrying both VEP `CSQ` and SnpEff `ANN`
  annotations on the same records.
- `gvr_filter()`: modular filtering by population allele frequency (gnomAD
  exomes, 1000 Genomes, ABraOM), ClinVar clinical significance, biotype,
  variant classification, genotype quality, and gene panel.
- `gvr_novel()`: dedicated subsetter for candidate novel variants (no dbSNP
  rsID and no allele frequency in any catalogue).
- `gvr_summary()`: 8-section cohort summary returning a list of `data.table`s,
  with optional Excel workbook, multi-page PDF, and self-contained interactive
  HTML dashboard (plotly + DT).
- `gvr_plot()`: `ComplexHeatmap`-based top-genes variant matrix across samples,
  coloured by variant classification.
- `gvr_lollipop()`: per-gene protein-domain lollipop plot with on-the-fly
  InterPro domain fetching (disk-cached) and a 52-palette colour system.
- `gvr_genepos.plot()`: per-gene gene-structure (cDNA) lollipop plot with
  exon/intron/UTR regions resolved via Ensembl REST or a local GTF file
  (disk-cached).
- `gvr_sum_plots()`: standalone PNG/SVG/PDF exports of all dashboard panels
  produced by `gvr_summary()`.
- `gvr_color_palette()`, `gvr_list_palettes()`: unified 52-palette colour
  system with colourblind-safe defaults shared across every figure.
- `gvr_panel_genes()`, `gvr_list_panels()`: curated disease gene-panel
  registry (breast cancer, Lynch syndrome, GIST, and others).
- `gvr_domain_cache_clear()`: cache management for InterPro domain data.
