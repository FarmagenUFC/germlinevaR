# germlinevaR 0.99.1

## New features

- New `gvr_hpo_genes()` helper: resolves one or more Human Phenotype Ontology
  (HPO) term id(s) to their annotated gene vectors using the HPO
  `phenotype_to_genes.txt` file. Lenient input forms are accepted and
  normalised (`"HP:0003002"`, `"hp:0003002"`, `"hp:3002"`, `"3002"`,
  `"0003002"` all resolve to canonical `"HP:0003002"`). Downloads and caches
  the annotation once per user under `tools::R_user_dir("germlinevaR", "cache")`;
  the cache is auto-refreshed after 30 days. Descendants are not expanded
  (exact-term match only). For offline / hermetic runs, point
  `options(gvr.hpo_path = )` at a local copy of the file.
- New `hpo` argument on `read.gvr()`, `read.gvr.snpeff()`, `read.gvr.dual()`,
  and `gvr_filter()`. Each accepts one or more HPO term id(s); resolved
  genes are UNION-ed with the existing `genes` and `panel` vectors before
  the downstream filtering machinery runs. Passing `hpo = NULL` (the
  default) is byte-identical to omitting the argument. Verbose mode prints
  a post-filter `hpo subset:` coverage line showing which HPO-derived genes
  survived filtering, alongside the existing `panel subset:` line.
- `gvr_filter()` gains a `panel` argument for API symmetry with the readers.
  Previously only `read.gvr()` / `read.gvr.snpeff()` / `read.gvr.dual()`
  accepted `panel`. `gvr_filter()` now unions `genes`, `panel`, and `hpo`
  identically, so a pre-loaded table can be sliced by disease panel or
  phenotype without going through the readers again.

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
