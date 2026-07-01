# germlinevaR 0.99.2

## Bug fixes

- `read.gvr.dual()` no longer fails on VCFs whose multi-ALT indel records
  share a left-anchored REF base. Previously, the reader built the SnpEff
  side-table keyed on `(chrom, pos, ref, alt)` using untrimmed coordinates,
  which collapsed real distinct ALT alleles onto the same key and threw a
  duplicate-key error at the join step. The reader now trims common
  left-anchor and right-anchor bases from each `(ref, alt)` pair before
  computing MAF-style coordinates (see below), so each ALT allele gets a
  unique key that matches the canonical form used by the VEP side.

## Breaking changes

- Coordinates emitted by `read.gvr()`, `read.gvr.snpeff()`, and
  `read.gvr.dual()` for indel ALT alleles are now MAF-normalised via a
  `bcftools norm`-equivalent trim before conversion. Downstream tables
  (`Start_Position`, `Reference_Allele`, `Tumor_Seq_Allele2`) reflect the
  trimmed representation. SNVs are unaffected.

  Example (`chr1:6095864`, `TCCCCCCCCCTGCCC`,
  `TCCCCCCCCCCTGCCC,TCCCCACCCCCTGCCC`):

  |                     | ALT1 (real 1 bp C ins) | ALT2 (compound event)  |
  |---------------------|------------------------|------------------------|
  | 0.99.1 and earlier  | `6095878 - C` (INS)    | `6095878 - C` (INS)    |
  | 0.99.2 default      | `6095874 - T` (INS)    | `6095869 - C` (INS)    |

  The pre-0.99.2 behaviour is available via `normalize_alleles = FALSE`
  on every reader, purely as an escape hatch for reproducing historical
  outputs; it re-exposes the multi-ALT collision bug and is not
  recommended for new analyses.

## New features

- New `normalize_alleles` argument (default `TRUE`) on `read.gvr()`,
  `read.gvr.snpeff()`, and `read.gvr.dual()`. Controls whether ALT-allele
  coordinates are MAF-normalised via `bcftools norm`-style left-and-right
  trimming before conversion. Setting `FALSE` restores pre-0.99.2 output
  for byte-level comparison against legacy tables.
- `read.gvr.dual()` output now carries an `snpeff_collisions_discarded`
  attribute (via `attr()`) when the SnpEff side-table has residual
  `(chrom, pos, ref, alt)` collisions after MAF-normalisation. The
  attribute is a `data.table` of the rows that were dropped by the
  Annotation_Impact-ranked deduper, so downstream code can audit which
  SnpEff annotations were superseded. Absent when no collisions occur.
  Residual collisions are expected on homopolymer runs where two ALT
  alleles collapse to the same MAF-key even after trimming
  (e.g. `chr10:68173832` in the accompanying test fixture).
- Verbose mode on `read.gvr.dual()` prints a new one-line summary of the
  Annotation_Impact-ranked dedupe pass on the SnpEff side-table when
  residual collisions are collapsed (format: `snpeff collisions: N
  MAF-key tuple(s) had M duplicate row(s); collapsed by Annotation_Impact
  rank`).

## Internal

- New unexported helper `.gvr_trim_alleles()` implements the bcftools-norm
  left-and-right trim used by `.gvr_coords()`.
- New unexported helper `.gvr_dual_dedupe_snpeff()` implements the
  Annotation_Impact-ranked collision-resolution pass on the SnpEff
  side-table, gated behind the `normalize_alleles` path.
- SnpEff-only reader's inlined `gvr_coords()` closure now delegates to
  the shared `.gvr_coords()` for consistency across all three readers.

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
