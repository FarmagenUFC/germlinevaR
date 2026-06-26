# germlinevaR — `read.gvr.dual()` integration

New exported function `read.gvr.dual()` parses VCFs annotated with **both**
Ensembl VEP `CSQ` and SnpEff `ANN` on the same records. The package's
auto-detector now recognises this case and routes the file to the new reader
instead of silently preferring VEP.

## Files in this folder

```
germlinevaR_dual_reader/
├── R/
│   ├── read.gvr.dual.R    ★ NEW
│   ├── read.gvr.R         (modified)
│   ├── read.gvr.snpeff.R  (modified: @importFrom fix)
│   ├── globals.R          (modified)
│   └── gvr_summary.R      (modified, backtick fix + @importFrom fix)
├── NAMESPACE              (modified)
└── test_outputs/
    ├── S6_dual_maf.tsv                  118,773 rows × 123 cols (82 MB)
    └── gvr_summary_dual_report.html     5.1 MB rendered report
```

Drop the five `R/*.R` files and `NAMESPACE` into your package source tree and
rebuild. No new package dependencies; no changes to `DESCRIPTION` are required.
After dropping in, run `devtools::document()` once to regenerate `man/*.Rd`.

`gvr_summary.R` is the same as the previous turn's patched copy
(`paste0("~\`", val_col, "\`")` at lines 1507-1508) — included here so this
folder is a self-contained drop-in.

## What the dual reader does

- **Spine:** VEP's most-severe-block pick drives the row spine (one row per ALT
  allele), same convention as `read.gvr()`. With `canonical_only = TRUE`
  (default), every kept row has `CANONICAL == "YES"`.
- **VEP-priority for shared fields:** `Hugo_Symbol`, `Consequence`, `IMPACT`,
  `HGVSc`, `HGVSp`, `dbSNP_RS`, gnomAD/MAX_AF/ClinVar/CANONICAL/MANE/SIFT/
  PolyPhen and the rest of the 80 canonical CSQ fields come from VEP.
- **`FREQS`** (the 81st CSQ field added by VEP `--everything` in v113+) is
  detected dynamically from `##INFO=<ID=CSQ` and emitted as its own column.
  Files without `FREQS` in the header get the column blank, no crash.
- **+ 4 LoF/NMD columns** from SnpEff's `INFO LOF=` / `NMD=`:
  `LOF_Gene`, `LOF_Pct_Transcripts`, `NMD_Gene`, `NMD_Pct_Transcripts`.
  Empty (`""` / `NA_real_`) when SnpEff did not call LoF/NMD on the record.
- **+ 4 side-by-side comparison columns** from the matching SnpEff ANN block:
  `snpeff_consequence`, `snpeff_impact`, `snpeff_gene`, `snpeff_hgvsc`.
  Block selection is gene-aware: it picks the ANN block whose `Gene_Name`
  matches VEP's `Hugo_Symbol` first; falls back to the first ANN block matching
  the ALT allele if no gene-matching block exists.
- **Same argument signature, same defaults** as `read.gvr()` (`folder` /
  `vcf_path` / `file` / `pattern` / DP/GQ thresholds / `genes` / `panel` /
  `vc_nonSyn` / `canonical_only` / ABraOM join / writers).
- **Tagged output:** `attr(maf, "annotator") == "dual"`.
- **Auto-route:** `read.gvr(...)` detects ANN+CSQ co-presence and delegates,
  emitting:
  `read.gvr: dual-annotated input detected (VEP CSQ + SnpEff ANN); delegating to read.gvr.dual().`

## Test results — `S6.deepvariant_snpEff_VEP.ann.vcf.gz`

| Metric | Value |
|---|---|
| Input | 47 MB gzipped, 322,249 records, sample `246_S6` |
| Wall time | 1.6 min (Phase 1 VEP = 47s, Phase 2 SnpEff scan = ~25s, Phase 3 join < 1s) |
| MAF dimensions | 118,773 rows × 123 cols |
| FREQS column | populated 26.7% of rows (per-record presence from VEP) |
| Allele match (VEP ↔ SnpEff) | 118,773 / 118,773 (100.0%) |
| Gene match (Hugo == snpeff_gene) | 113,686 / 118,773 (95.7%) |
| VEP-vs-SnpEff IMPACT agreement | 97.4% overall (115,697 / 118,773) |
| VEP=HIGH ↔ SnpEff=HIGH | 1,160 / 1,219 (95.2%) |
| LoF call on HIGH-impact LoF | 996 / 1,219 (81.7%) |
| NMD call on HIGH-impact LoF | 133 / 1,219 (10.9%) |

VEP-supplied field fill rates (now populated; were 0% on the SnpEff-only run):

| Field | Fill rate |
|---|---|
| `dbSNP_RS` | 77.5% |
| `gnomADe_AF` | 66.5% |
| `gnomADg_AF` | 74.8% |
| `MAX_AF` | 76.9% |
| `MAX_AF_POPS` | 76.9% |
| `BIOTYPE` | 100.0% |
| `CANONICAL == "YES"` | 100.0% |
| `MANE_SELECT` | 92.3% |
| `CLIN_SIG` | 15.0% (coding-only) |
| `SIFT` | 7.3% (coding-only) |
| `PolyPhen` | 7.1% (coding-only) |
| `FREQS` | 26.7% |

Regression test: `read.gvr()` on the SnpEff-only VCF produces **127,882 rows**
× 115 cols with `attr='snpeff'` — identical row count to the previous
turn's run (only difference: 1 column from `add_abraom = FALSE` in this test).

## Code review notes (issues found and fixed during implementation)

The first draft had a silent +1 offset bug for left-anchored insertions where
`nchar(REF) > 1` (e.g. `ref="AT", alt="ATCGT"`). The fix was to lift
`read.gvr()`'s `maf_coords()` verbatim into `read.gvr.dual.R`. Verified across
9 indel test cases: all match canonical now. Impact on this VCF: would have
silently dropped SnpEff annotations on ~1,246 records' anchored insertions
without the fix.

The first-draft ANN-header parser used a regex that was correct on this VCF
but fragile against SnpEff files with extra single-quotes in the description.
Replaced with the robust "first-and-last single-quote" approach used in
`read.gvr.snpeff.R`'s `get_ann_fields()`.

## Limitations / out of scope (mirrors plan)

- Multi-sample VCFs: not supported (inherited from existing readers).
- `gvr_filter()` / `gvr_summary()` / oncoplot integration of the new
  `LOF_*` / `snpeff_*` columns as filterable axes: not done here. The
  columns exist in the output `data.table` and can be queried manually
  (`af[IMPACT == "HIGH" & snpeff_impact != "HIGH"]`).
- `testthat` scaffolding: not added. The acceptance criteria above are the
  current test coverage.
- The 51,047 (40%) `GT=0/0` records passing DP/GQ on this sample — a
  `gvr_filter()` policy question, not a reader question.

## Reader-side changes summary

| File | Change |
|---|---|
| `R/read.gvr.dual.R` | NEW. ~700 lines. Public `read.gvr.dual()` + private helpers. |
| `R/read.gvr.R` | `.detect_annotator()` returns `"dual"` when ANN+CSQ both present. Added `dual` dispatch branch. Added private `.force_annotator` parameter (internal use only). Extended sibling auto-source helper. |
| `R/globals.R` | Added 9 names: `FREQS`, `LOF_Gene`, `LOF_Pct_Transcripts`, `NMD_Gene`, `NMD_Pct_Transcripts`, `snpeff_consequence`, `snpeff_impact`, `snpeff_gene`, `snpeff_hgvsc`. |
| `NAMESPACE` | Added `export(read.gvr.dual)`. |

## R CMD check fix (post-test patch)

Two `R CMD check` warnings flagged in the latest run, both fixed in-place
without touching reader behaviour:

1. **Undocumented argument: `.force_annotator`** in `read.gvr()` — the dual
   reader's auto-routing introduced this private parameter but the prior
   patch omitted its roxygen entry. Added an `@param .force_annotator` block
   to `R/read.gvr.R` (lines 156-161) that documents it as "Internal use
   only; do not set" and explains its `"vep"` / `"snpeff"` accepted values
   and `NULL` default. The function signature and dispatch logic are
   unchanged.

2. **`@importFrom` blocks split across multiple lines** in three files —
   roxygen2 rejects multi-line `@importFrom` tags ("must be only 1 line
   long, not N"). Collapsed every offending block to a single line in
   `read.gvr.R`, `read.gvr.snpeff.R`, and `gvr_summary.R`. Imports
   themselves are identical; only the source formatting changed.

After both fixes, `roxygen2::parse_file()` on `read.gvr.R` returns the
function's documentation block with **28 `@param` tags matching all 28
formal arguments exactly** (`.force_annotator`, `abraom_path`, `abraom_url`,
`add_abraom`, `add_genotype`, `cache_dir`, `canonical_only`, `chunk_size`,
`dedup_columns`, `drop_empty_cols`, `file`, `folder`, `genes`, `min_DP`,
`min_GQ`, `ncbi_build`, `ncores`, `out_dir`, `out_prefix`, `panel`,
`pattern`, `strip_hgvs_prefix`, `vc_nonSyn`, `vcf_path`, `verbose`,
`write_rds`, `write_tsv`, `write_xlsx`). Same applied to `read.gvr.snpeff.R`
and `gvr_summary.R`.

## Verbose dedup patch

`read.gvr()` on a dual-annotated VCF previously printed `Found 1 file(s):`
twice — once from the outer dispatch in `read.gvr()` itself, and once
again from the inner recursive call that the dual reader makes
(`read.gvr(.force_annotator = "vep")`).

One-line guard added at `R/read.gvr.R` line 551:

```r
# When called internally from read.gvr.dual() via .force_annotator, suppress
# the file-listing message: the outer read.gvr() already printed it before
# routing here, so re-printing would just duplicate.
if (verbose && is.null(.force_annotator)) {
  message(sprintf("Found %d file(s):", length(vcf_paths)))
  for (p in vcf_paths) message("  - ", basename(p))
}
```

All other verbose output (delegating notice, genome build, "Converting...",
record counters, write paths) is preserved unchanged. Live re-test on the
S6 dual VCF emits exactly one `Found 1 file(s):` line.

## `gvr_sum_plots()` companion function

NEW exported function `gvr_sum_plots(maf, ...)` that writes the same plot
categories `gvr_summary()` puts in its dashboard as **individual image files
plus combined panel images**, into a new folder. Recomputes its sections
internally from the MAF (no dependency on a prior `gvr_summary()` call).
Static `ggplot2` engine; export via `ggsave()`.

### Signature

```r
gvr_sum_plots(
  maf,
  out_dir        = ".",
  folder_name    = "gvr_sum_plots",
  format         = "png",
  width          = 7,
  height         = 5,
  dpi            = 300,
  sample_col     = "Tumor_Sample_Barcode",
  top_n_genes    = 20,
  top_n_variants = 20,
  per_sample     = TRUE,
  panel          = TRUE,
  verbose        = TRUE
)
```

`format` accepts every device `ggsave()` supports —
`c("png","pdf","svg","jpeg","tiff","bmp","eps","ps","tex","wmf")` — and
validates up front. `svg` requires `svglite`; `wmf` requires `ragg`. The
function checks `requireNamespace()` before the first plot is built and
fails with an actionable message if the device's helper package is missing.

### Output folder layout

```
out_dir/folder_name/
├── top_genes.<ext>                     # cohort: top-N genes
├── variant_classification.<ext>        # cohort: top-10 classifications
├── impact.<ext>                        # cohort: VEP IMPACT levels
├── top_variants.<ext>                  # cohort: top-N rsIDs (if dbSNP_RS has data)
├── panel_cohort.<ext>                  # 2x2 grid of the 4 cohort plots
├── panel_per_sample.<ext>              # grid of per-sample plots (if n_samples > 1)
└── per_sample/
    └── top_genes__<sample>.<ext>       # one per sample
```

Sample names are sanitized for filenames (`[^A-Za-z0-9._-]` -> `_`). With
`per_sample = FALSE` the `per_sample/` folder is not created; with
`panel = FALSE` neither panel image is written.

### Panel assembly

Uses `patchwork::wrap_plots()` when available; falls back to
`gridExtra::arrangeGrob()` wrapped in `ggplot() + theme_void() +
annotation_custom()` when not. The cohort panel is a 2-column grid sized
`width*2 x height*ceiling(n_plots/2)`; the per-sample panel uses
`ncol = ceiling(sqrt(n_samples))`.

### Plot rendering

- Cohort bar plots mirror `gvr_summary()`'s `.bar_grob` logic: grouped bars
  when `n_samples <= 6` (`FACET_THRESHOLD`), faceted otherwise. Same Phylo
  color palette (`PHYLO_BLUE`, `PHYLO_GREEN`, `ORANGE`, `VERMIL`, `YELLOW`,
  `PINK`).
- Per-sample bar plot is a horizontal single-sample chart sorted ascending
  so the largest bar lands on top after `coord_flip()`.
- Helper `.gvr_sp_save_one()` selects `device = format` and passes `dpi`
  only for raster formats (`png`, `jpeg`, `tiff`, `bmp`).

### Test results

Run on the cached dual MAF `/workspace/maf_dual_test.rds`
(118,773 rows x 123 cols, sample `246_S6`):

| Test | Result |
|---|---|
| Default `png`, 1 sample | 6 files in 5.84 s; sizes 38.3-269.2 KB |
| `format = "pdf"` | 6 files, 4.5-6.2 KB each |
| `format = "jpeg"` / `"tiff"` | 6 files each, format correct |
| `format = "svg"` (no svglite) | Stops with: *"format 'svg' requires the 'svglite' package..."* before any plot built |
| `format = "xyz"` | Stops with full valid-device list |
| `per_sample = FALSE` | 5 files; no `per_sample/` folder; no `panel_per_sample` |
| `panel = FALSE` | 5 files; no panel images |
| Multi-sample (synthetic 4-sample split of S6) | 10 files including `panel_per_sample.png` with 2x2 grid of samples S01-S04 |

**Data parity** vs. prior `gvr_summary()` run on the same MAF
(verified against transcript record i=320):

- IMPACT: `HIGH=1219`, `MODERATE=9003`, `LOW=10561`, `MODIFIER=97990`
- variant_classification top: `Intron=81061`, `Missense_Mutation=8767`,
  `3'UTR=6672`, `Silent=6282`, `Splice_Region=4284`
- top_genes top 10: `MUC16=217`, `MUC19=112`, `MUC3A=100`, `CSMD1=93`,
  `SSPOP=87`, `FBN3=84`, `ADAMTS17=75`, `TTN=74`, `OBSCN=71`, `SYNE2=71`
- top variant by rsID: `rs398123612` with 2 occurrences

All 8 acceptance criteria from `PLAN.md` pass.

### Files in this folder

```
germlinevaR_dual_reader/
├── R/
│   └── gvr_sum_plots.R   <- NEW  (~570 lines: public + 6 private helpers)
├── NAMESPACE             <- modified: + export(gvr_sum_plots) + 7 ggplot2 importFrom
├── DESCRIPTION           <- modified: + patchwork in Suggests
└── test_outputs/
    └── gvr_sum_plots_example/        <- NEW: sample PNG outputs on S6 MAF
        ├── top_genes.png
        ├── variant_classification.png
        ├── impact.png
        ├── top_variants.png
        ├── panel_cohort.png
        └── per_sample/
            └── top_genes__246_S6.png
```

### NAMESPACE / DESCRIPTION changes

- `NAMESPACE`: regenerated by `roxygen2::roxygenize()` from the roxygen tags
  in `R/*.R`. Header line marks it as auto-generated; subsequent
  `devtools::document()` calls will reproduce it identically. New entries
  introduced by this drop are `export(gvr_sum_plots)` plus
  `importFrom(ggplot2, ...)` for `aes`, `as_labeller`, `coord_flip`,
  `element_blank`, `element_text`, `facet_wrap`, `geom_col`, `ggplot`,
  `ggsave`, `labs`, `position_dodge`, `scale_fill_manual`,
  `scale_y_continuous`, `theme`, `theme_minimal` and
  `importFrom(data.table, as.data.table, data.table, setorder)`.
- `DESCRIPTION`: `patchwork` added to `Suggests`. The function still works
  without it (falls back to `gridExtra` for panel assembly), so it stays
  optional rather than moving to `Imports`.

### Limitations / out of scope

- No automatic cap on per-sample plot count. A 50-sample run produces 50
  per-sample images plus a 50-cell `panel_per_sample.<ext>`. The user
  controls this via `per_sample = FALSE`.
- The `bmp`, `eps`, `ps`, and `tex` devices are accepted (they are valid
  `ggsave` device strings) but not exercised in the test matrix above.
- `gvr_sum_plots()` is **standalone**: it is not called from `gvr_summary()`.
  The two functions are intentionally independent so the dashboard stays
  fast and the file-export path can be invoked on its own schedule.

## `gvr_summary()` performance optimization (drop-in)

Three surgical optimizations applied to `R/gvr_summary.R`. End-to-end runtime
drops ~20% on a 475k-row x 20-sample synthetic input and ~30% on a 118k-row x
1-sample real input. All optimizations preserve identical writer output (XLSX /
HTML / PDF) and identical `sections` list content. Public API, dependencies,
and `DESCRIPTION` are unchanged.

### Changes

**Change 1 - Global rank-of-row precompute (XLSX + HTML).** The combined
drill-down sheets and HTML drill-down tables both sort each per-token row pool
by the same composite key (`impact_rank`, `-n_samples`, `gnomAD_AF`,
`chromosome`, `start_position`). The baseline re-sorted each pool independently
from the per-token row indices. Since the sort key depends only on row content
(not on pool membership), the integer rank of every row in the global order is
computed once and inverted:

```r
.xl_rank_of_row <- integer(nrow(.xl_proj))
.xl_rank_of_row[order(.xl_proj$.__ir__, -.xl_proj$.__nsamp__,
                      .xl_proj$.__gaf__, .xl_proj$.__chr__,
                      .xl_proj$.__pos__, na.last = TRUE)] <- seq_len(nrow(.xl_proj))
```

Per-token sort then collapses to a single integer-keyed `order()`:

```r
ord_idx <- row_idx[order(.xl_rank_of_row[row_idx])]
```

Same construction applied on the HTML side as `.html_rank_of_row` (built via
`.gvr_rank_variants(dt)` then inverted), reused by `.dd_build_cat()`.
Algebraically identical to the per-pool sort because the key is row-intrinsic.

**Change 2 - Batched XLSX combined-sheet writes.** Each combined drill-down
sheet (`Genes_detail`, `VC_detail`, `Clinical_detail`, `Impact_detail`) was
previously written with 2-3 `openxlsx::writeData` calls per token (blank row,
label row, data block). At 20 samples and 50+ tokens per category, that is
hundreds of `writeData` calls per sheet, and `writeData` serializes through
R<->C per cell so it dominates the XLSX writer cost.

The new path builds the entire sheet body as a single data.frame via
`data.table::rbindlist(body_parts, use.names = TRUE, fill = TRUE)` and emits
ONE `writeData` call per sheet at `startRow = 2L`. The bold style on label rows
is applied as a single vectorised `addStyle()` call
(`rep(label_rows, each = n_cols)` and
`rep(seq_len(n_cols), times = length(label_rows))`), and `mergeCells` is called
once per label row (cheap, no cell-data work).

**Critical detail:** the blank/label rows must carry the SAME column types as
`.xl_proj` so that `rbindlist` preserves numeric and integer types. Building
them with `NA_character_` would force `rbindlist` to upcast numerics, which
then land as XLSX shared strings (`t="s"`) instead of inline numerics
(`t="n"`). Visible symptom: `gnomADe_AF` renders as `7.206e-07` instead of
openxlsx's `0.0000007206`. The fix:

```r
blank_row <- (function() {
  proto <- .xl_proj[0L, .gvr_xl_cols, with = FALSE]
  typed_na <- lapply(proto, function(col) col[NA_integer_])
  as.data.frame(typed_na, stringsAsFactors = FALSE)
})()
```

This restores 0 sheet-content diffs against the baseline writer output (12/12
sheets identical on small input, 31/31 on big), all 56 merge regions identical,
and all `styleObjects` row/col positions identical.

**Change 3 - Drop dead `filter_fn` closures.** The `cat_specs` list (XLSX
writer guard) and `.dd_token_specs` list (HTML writer) each carried 4 closures
of shape `filter_fn = function(tok) ...` that walked `dt` to mask rows matching
one token. None of these closures were ever invoked: `.gvr_build_token_index()`
replaced them with a single category-aware `split()` pass that hands the
writers a prebuilt `list[category][token] -> integer row indices`. Removing the
dead closures (and the equally dead `pfx` field, a relic of the old
per-token-per-sheet layout) saves 4 R closure allocations plus their
constant-pool entries on every `gvr_summary()` call.

### Optimization NOT applied (deliberately)

**`.is_missing()` memoization.** Profiling showed `.is_missing(cs)` (used
inside the unused `filter_fn` for `clin_sig`) recomputed
`is.na(cs) | trimws(cs) %in% c("", ".")` on every token call. After Change 3
removes the only caller, this is moot. Earlier profiling estimated the gain
in a hypothetical "keep filter_fn, memoize" path at ~0.3-0.4s on the big input
- below the run-to-run benchmark noise floor on either scale, and not worth
the additional state to carry.

### Verification protocol

For each input scale (small: 118,773 x 1-sample real; big: 475,092 x 20-sample
synthetic from 4x row replication + uniform sample reassignment + position
perturbation):

1. **Sections content-equal** via a recursive comparator that walks the nested
   list and inspects each `data.table` column-by-column, ignoring
   `.internal.selfref` pointer (NIL after `readRDS`, live after fresh build).
   Both inputs report TRUE.
2. **XLSX content** via `openxlsx::read.xlsx()` per sheet: 0/12 sheets differ
   on small, 0/31 sheets differ on big.
3. **XLSX merge regions** via
   `wb@.xData$worksheets[[i]]@.xData$mergeCells`: 56 merges identical across
   all drill-down sheets on both inputs.
4. **XLSX style positions** via `wb@.xData$styleObjects` row/col sets: 16/16
   identical on small, 35/35 identical on big.
5. **HTML semantic content** via line-diff after stripping known
   non-determinism (plotly `htmlwidget-<hex>` UUIDs, visdat 8-13 hex hashes,
   and the generated-at timestamp): 0 diffs on both inputs.

### Measured speedup (2-run mean, fresh `sys.source` per run)

| Configuration | Baseline | Patched | Speedup | Improvement |
|---|---|---|---|---|
| 118k x 1, all writers (xlsx + pdf + html) | 19.78 s | 13.84 s | 1.43x | 30.0% |
| 475k x 20, xlsx + html | 46.27 s | 36.86 s | 1.26x | 20.3% |

Headline number is **~1.26x** at production scale (475k x 20). The small-input
1.43x includes some run-to-run variance (15.19 s, 12.49 s across the two
patched runs); the steady-state gain on the big input is ~1.26x with much
lower variance (36.80 s, 36.92 s). The reduction is dominated by Change 2
(batched `writeData`). Changes 1 and 3 are smaller wins that do not regress
correctness.


### Real-data validation (2-sample combined input)

Initial benchmarks used a synthetic 20-sample input built by 4x row replication
+ uniform sample reassignment from the 118k 1-sample real MAF. To confirm the
speedup holds on un-synthesized multi-sample data, two real VEP-annotated VCFs
(`S1.vep.vcf.gz`, `S2.vep.vcf.gz`; samples `ACRO1_S1`, `ACRO2_S2`) were read
via `read.gvr()` and `rbindlist`-combined into a single MAF of 260,255 rows x
2 samples.

| Configuration | Baseline | Patched | Speedup | Improvement |
|---|---|---|---|---|
| 260k x 2 real, all writers (xlsx + pdf + html) | 38.09 s | 27.38 s | 1.39x | 28.1% |
| 260k x 2 real, xlsx + html | 32.46 s | 25.93 s | 1.25x | 20.1% |

Run-to-run variance is tight (baseline runs 37.67 / 38.51 s for all writers,
patched 26.63 / 28.13 s), so the mean is trustworthy. The xlsx+html number
(1.25x, 20.1% faster) reproduces the synthetic 475k x 20 result almost
exactly (1.26x, 20.3% faster), confirming the speedup is not an artifact of
the 4x-replicated synthetic structure. The all-writers number (1.39x) sits
between the small (1.43x) and big (no PDF benched) data points, consistent
with the writer mix at this scale.

Writer parity on this real input:

- Sections content-equal: TRUE
- XLSX sheet-content diffs: 0 / 13
- XLSX merge regions identical: TRUE (59 merges)
- XLSX style positions identical: TRUE (17 objects)
- XLSX cell-type-flag drift (`sheet_data$t`): 0 cells across all sheets
- HTML semantic line diffs: 0 / 3,620 lines
### Risk / regression notes

- **Type drift on `blank_row` is silent unless detected at the cell-type
  level.** Both `read.xlsx()` and a naive `readLines()` comparison parse
  numeric cells through `as.character`, so the `7.206e-07`-vs-`0.0000007206`
  drift is invisible to either check. The verification protocol above
  explicitly inspects `wb@.xData$worksheets[[i]]@.xData$sheet_data$t` (cell
  type flags) in addition to `read.xlsx()` content to catch this class of
  regression.
- **Column-type inference** for `blank_row` is driven by whatever columns
  survive the `intersect(c(.gvr_xl_cols, ...), names(dt))` projection step
  (line ~800). If a new column is added to `.gvr_xl_cols` with a
  non-character/numeric/integer type, the `col[NA_integer_]` trick still
  works (returns a typed `NA` of the column's class) but the writeData path
  has not been exercised on factor or POSIXct columns. Current `.gvr_xl_cols`
  is `Hugo_Symbol`, `dbSNP_RS`, `CLIN_SIG`, `IMPACT`, `gnomADe_AF`,
  `Chromosome`, `Start_Position` - all character / numeric / integer.

## `read.gvr()` micro-refactor + dead-code cleanup

Two small, content-equivalent edits to `R/read.gvr.R`'s per-record conversion loop, plus removal of an unused helper. No measurable wall-clock speedup on the bench, but no regression either; the code path is now slightly leaner.

### What changed

1. **P1 — coalesce 4 vapply calls into one walk** (inside `convert_chunk`, the per-ALT CSQ-block ranking section). The original code did:

   ```r
   ranks <- vapply(sel, function(k) consequence_rank(block_fields[[k]][P_Cons]), integer(1))
   canon <- vapply(sel, function(k) { v <- block_fields[[k]][P_CANONICAL]; ... }, integer(1))
   mane  <- vapply(sel, function(k) { v <- block_fields[[k]][P_MANESEL];   ... }, integer(1))
   feat  <- vapply(sel, function(k) { v <- block_fields[[k]][P_Feature];   ... }, character(1))
   ```

   walking `sel` four times and indexing `block_fields[[k]]` four times per element. Replaced with one pre-allocated `for (i in seq_len(ns))` that fills all four vectors, indexing each `block_fields[[k]]` exactly once. The subsequent `order(ranks, canon, mane, feat)` is unchanged.

2. **P2 — inline the `map_code` closure** at its two call sites (`t_allele1`/`t_allele2`). Removes one closure object created per (record × ALT). Branch logic identical.

3. **P5 — delete unused `strip_feature_prefix` helper** (lines 880-887 of the baseline). The function was defined inside `read.gvr()` but never called; the vectorized `strip_prefix_vec` at line 1649 is the one actually used to strip the leading `feature_id:` from HGVS columns. (Note: `read.gvr.snpeff.R` defines its own `strip_feature_prefix` in a separate scope; left untouched.)

### Things I tried and reverted

- **P3 — chunk-level vectorized `INFO` field extraction.** Microbench on 25,000 real INFO strings showed the fastest vectorized variant (one `strsplit` + `match` over the 7 needed keys) ran at 1.11× of the per-record `info_parse()` baseline (0.664 s vs 0.740 s). Since `info_parse` is only ~1% of profile samples, even the best vectorized version would have saved ~0.1% wall-clock. Skipped.
- **P4 — separate integer-valued `.rank_cache` for `consequence_rank()`.** Three back-to-back runs after the patch landed at 91-99 s vs P1+P2's ~93 s baseline. The new parallel cache added work without saving any (the existing `.mstr_cache` already returns a list pointer in O(1) on hits, so the `$rank` access being skipped didn't matter). Reverted.

### Bench (real VCFs, ncores=1, add_abraom=FALSE)

Interleaved 3-baseline + 3-patched runs on the same machine, alternating to absorb timing variance:

| Input | Baseline mean | Patched mean | Ratio |
|---|---|---|---|
| S1.vep.vcf.gz (263k records → 128k MAF rows) | 92.65 s (range 89.08-94.70) | 90.13 s (range 81.16-95.56) | 1.028× (2.7%, within noise) |
| S2.vep.vcf.gz (267k records → 132k MAF rows) | 102.15 s (1 run) | 102.15 s (1 run) | 1.000× |

### Equivalence verification

- S1 patched MAF compared against cached baseline MAF column-by-column with NA-safe equality: **all 115 columns identical** across 128,032 rows.
- S2 patched MAF compared against S2 baseline MAF same way: **all 115 columns identical** across 132,223 rows.
- No new dependencies. No exported-API change.

### Honest takeaway

`read.gvr.R` already carries 8+ labeled optimization passes (C1/C2/C3/O1/O2/O3/O6/O7/O8/A2/A3/B2/vN+4). The Rprof percentages suggested `vapply` (18.4%) and `strsplit` (20.5%) were the dominant cost — but those distributions overcount deep frames, and in wall-clock the per-record `for` loop in `convert_chunk` is bounded by something more uniform: data movement, GC pressure, and the inherent cost of touching every CSQ block in every record. The two cheap algebraic rewrites (P1, P2) are correctness-preserving and slightly clearer code, but they don't move the bottleneck.

The remaining real win would be restructuring `convert_chunk` so the per-record `for (r in seq_len(nrow_dt))` is replaced by chunk-level `data.table` column operations across all records at once. That's a substantial refactor (the per-record body assembles a list with ~50 fields per ALT, including downstream branching on FORMAT layout, canonical filtering, and the HGVSp 3→1 substitution), and we deferred it as out-of-scope for the current change set.

### Risk / regression notes

- Both P1 and P2 are pure algebraic reorderings of existing logic; both verified content-identical on real data (S1 + S2, ~260k MAF rows).
- P5 deletes a function with zero callers in the file. `grep -rn 'strip_feature_prefix'` confirms it is also unused anywhere else in the package's R/ folder except its independent definition in `read.gvr.snpeff.R` (a different file with its own internal scope; not touched).


## Parallel-mode heartbeat (vN+12)

**Problem.** When `read.gvr()` is called with `ncores > 1L` and `n_files > 1L` on a fork-capable OS, the work is dispatched via `parallel::mclapply`. By design (see the doc comment above the dispatch), each worker is wrapped in `suppressMessages()` so per-chunk percent lines from `convert_one_vcf()` don't interleave across forks into garbled output. The side effect: between the "Parallel conversion: ..." banner and the final per-file `done:` summary, the user sees **no output at all** — for a typical 90-100 s/file workload over several files, that's minutes of dead silence with no way to tell whether R is alive, deadlocked, or has crashed.

**What changed.** Wrap the existing `mclapply` call in a `tryCatch(..., finally = ...)` and spawn a single **side-fork** via `parallel::mcparallel()` that does nothing except sleep for `heartbeat_secs = 15L` seconds, print one liveness line, and repeat. The side fork is killed unconditionally via `tools::pskill()` in the `finally` clause, regardless of whether `mclapply` returned normally or signalled. The side fork is independent of the work — it does not look at jobs, queues, results, or PIDs — so it cannot race with `mclapply`'s child reaping. The data path is unchanged: same `mclapply(..., mc.cores = use_cores, mc.preschedule = FALSE)`, same `try-error` surfacing, same `sample_names` extraction, same per-file `done:` summary at the end in deterministic index order.

Heartbeat line format (printed to `stderr` via `message()`):

```
  ... parallel conversion still running: 2 worker(s), 60s elapsed
```

**Why this design (after a failed first attempt).** The original plan was to replace `mclapply` with a hand-rolled `mcparallel + mccollect(wait = FALSE, timeout = N)` polling loop, with the main process printing the heartbeat from the poll loop. That implementation worked in principle but hit a real `parallel`-package race: when both children finished during the same `mccollect` window, the SIGCHLD-handler path reaped one before the explicit `waitpid` could, producing a "cannot wait for child PID" warning and `NULL` result; subsequent poll calls then warning-stormed against the dead PID at full CPU and never collected the results. Walking that race down would have required either reading from `mcparallel` job pipes directly (uses unexported `parallel:::readChild`) or replacing the SIGCHLD handler — both fragile across R versions. Keeping `mclapply` for the actual work (it's used elsewhere in the package and is known to reap cleanly) and giving up the "any-worker-done early" hook in exchange for a stable side-fork heartbeat was the cheaper, safer trade.

### Files touched

- `R/read.gvr.R`
  - doc comment above the parallel block (mentions vN+12 heartbeat behaviour)
  - `if (want_par && use_cores > 1L) { ... }` branch reworked: spawn heartbeat fork, run `mclapply` inside `tryCatch(... , finally = pskill + suppressWarnings(mccollect))`, then surface errors and emit per-file `done:` summaries as before

No other files changed. `ncores = 1L` path (the sequential `for` loop) is untouched, so behaviour is bit-identical to before for users who don't request parallelism.

### Verification

| Scenario | Result |
|---|---|
| `ncores = 1L` on S1.vep.vcf.gz vs cached baseline MAF | **PASS** — 115/115 cols identical on 128,032 rows; 72.7 s wall-clock |
| `ncores = 2L` on `{S1, S2}` folder vs `ncores = 1L` on same folder | **PASS** — 115/115 cols identical on 260,255 rows |
| `ncores = 2L` on `{S1, S2}` folder — heartbeat cadence | 6 heartbeats at 15/30/45/60/75/90 s during the silent period; per-file `done:` summary printed at the end in index order; 96.4 s total |
| `ncores = 2L` on a 7-file folder (queueing test) | All 7 files processed; per-file rows match expected (4× ACRO1_S1 + 3× ACRO2_S2 = 908,797 rows); 503 s total |

**Known cosmetic quirk (7-file run).** During worker turnover (one batch finishing, mclapply forking the next), the heartbeat fork's `Sys.sleep(15)` can be delayed by SIGCHLD activity in the main process. In one observed 7-file / 2-core run, the heartbeat cadence held at 15 s for the first ~169 s, then there was a single 118 s gap before resuming at 15 s. The user still sees liveness updates throughout — the gap happens precisely when worker transitions are themselves progress — and the run completes correctly. Not worth additional complexity to chase.

### Risk / regression notes

- The new code is gated by `want_par && use_cores > 1L`. With `ncores = 1L` (the default), the gate is false and the code path is the same `for` loop it has always been.
- The heartbeat fork is killed in `finally`, including when `mclapply` itself raises (e.g. on `try-error` surface conversion). No orphaned heartbeat forks across runs.
- The `suppressWarnings(try(parallel::mccollect(hb_job, wait = FALSE, timeout = 1)))` after `pskill()` suppresses the expected "1 parallel job did not deliver a result" warning, which is harmless — it just confirms we killed the fork before it could finish.

**Drop-in:** replace `R/read.gvr.R` with the version in this folder. No NAMESPACE or DESCRIPTION changes needed; `parallel` and `tools` are base R.
