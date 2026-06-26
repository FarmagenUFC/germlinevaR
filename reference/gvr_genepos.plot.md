# Gene-track lollipop plot on a cDNA axis

Draws a per-gene track plot with exon, intron, and UTR segments for one
transcript, overlaid with lollipops placed on a cDNA-position x-axis
using the `HGVSc` field of each table row. Colours follow
`Variant_Classification` using the same palette set as
[`gvr_lollipop()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_lollipop.md).

## Usage

``` r
gvr_genepos.plot(
  gvr,
  gene,
  transcript_id = NULL,
  vc_keep = NULL,
  color_by = c("vc", "region"),
  intron_scale = c("fixed", "proportional", "log"),
  intron_visual_bp = 200L,
  utr_visual_bp = NULL,
  label_top = 5L,
  hotspot_window = 20L,
  hotspot_min_n = 4,
  stem_alpha = 0.6,
  point_size = 3,
  bar_color = "grey85",
  bar_border = "grey40",
  utr_color = "grey70",
  intron_color = "grey60",
  exon_color = "steelblue",
  variant_palette = "gvr",
  base_size = 12,
  axis_text_size = 11,
  axis_title_size = 12,
  ensembl_release = NULL,
  gtf_path = NULL,
  cache_dir = NULL,
  out_dir = ".",
  out_subdir = "gvr_genepos",
  out_prefix = NULL,
  format = c("png", "svg", "pdf", "tiff"),
  width = 10,
  height = 4,
  dpi = 300,
  verbose = TRUE
)
```

## Arguments

- gvr:

  A MAF-like data.table produced by
  [`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md).
  Required columns: `Hugo_Symbol`, `Transcript_ID`, `HGVSc`,
  `Variant_Classification`, `Tumor_Sample_Barcode`, `NCBI_Build`.
  `MANE_SELECT` and `CANONICAL` are used when `transcript_id` is
  auto-resolved.

- gene:

  Character(1). HGNC symbol (matches `Hugo_Symbol`).

- transcript_id:

  Character(1) or `NULL`. Ensembl stable id (`ENST...`). `NULL` triggers
  auto-resolution described above.

- vc_keep:

  Character vector or `NULL`. If non-`NULL`, only variants whose
  `Variant_Classification` is in `vc_keep` are kept. Otherwise the same
  default filter as
  [`gvr_lollipop()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_lollipop.md)
  is applied (non-synonymous plus splice plus UTR).

- color_by:

  Character(1). One of `"vc"` (colour by `Variant_Classification`;
  default) or `"region"` (colour by track region: cds / utr5 / utr3 /
  splice / intron). Currently only `"vc"` has a stable visual contract.

- intron_scale:

  Character(1). One of `"fixed"` (default; constant visual width per
  intron, `intron_visual_bp` pixels), `"proportional"` (true bp),
  `"log"` (log10-scaled).

- intron_visual_bp:

  Integer(1). Visual width of every intron when
  `intron_scale == "fixed"`. Default `200`.

- utr_visual_bp:

  Integer(1) or `NULL`. When non-`NULL`, total UTR visual length is
  rescaled to this many bp regardless of real UTR length (useful for
  very long UTRs). `NULL` (default) means real length.

- label_top:

  Integer(1). Number of top-counted cDNA positions to label (matches
  [`gvr_lollipop()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_lollipop.md)).

- hotspot_window:

  Numeric(1). Sliding-window width (bp on cDNA axis) for hotspot
  detection. Default `20`.

- hotspot_min_n:

  Numeric(1). Minimum distinct cDNA positions inside a window for it to
  be drawn as a hotspot band. Default `4`. Pass `Inf` to disable.

- stem_alpha:

  Numeric(1). Lollipop stem opacity. Default `0.6`.

- point_size:

  Numeric(1). Dot size. Default `3`.

- bar_color, bar_border:

  Character(1). Fill / border colour for the protein-body bar drawn over
  coding exons. Default `"grey85"` / `"grey40"`.

- utr_color, intron_color, exon_color:

  Character(1). Track region fill colours. Default `"grey70"` /
  `"grey60"` / `"steelblue"`.

- variant_palette:

  See
  [`gvr_lollipop()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_lollipop.md)
  for the full grammar (`"gvr"`, palette name, or named override
  vector).

- base_size, axis_text_size, axis_title_size:

  Numeric. Text sizes.

- ensembl_release:

  Integer or `NULL`. Future use; currently unused (REST `latest` is
  always queried).

- gtf_path:

  Character(1) or `NULL`. Offline GTF override. If non-`NULL`, REST is
  not called.

- cache_dir:

  Character(1), `NULL`, or `FALSE`. Cache directory. `FALSE` disables
  on-disk caching.

- out_dir:

  Character(1) or `NULL`. Output directory root.

- out_subdir:

  Character(1). Subfolder under `out_dir`. Default `"gvr_genepos"`.

- out_prefix:

  Character(1) or `NULL`. File prefix; default
  `paste(gene, transcript_id, sep = "_")`.

- format:

  Character(1). One of `"png"`, `"svg"`, `"pdf"`, `"tiff"`. Default
  `"png"`.

- width, height, dpi:

  Numeric. Plot size and resolution.

- verbose:

  Logical(1). If `TRUE`, emit progress messages.

## Value

A ggplot object, returned invisibly. As a side effect the plot is
written to `<out_dir>/<out_subdir>/<out_prefix>.<format>` when `out_dir`
is non-`NULL`.

## Details

Companion to
[`gvr_lollipop()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_lollipop.md)
which draws protein-domain rectangles on a protein-position axis. This
function instead draws gene structure on a cDNA axis.

## Transcript resolution

If `transcript_id` is `NULL`, the chosen transcript is, in order: the
first non-empty `MANE_SELECT` among table rows for `gene`; otherwise the
first `CANONICAL == "YES"` row; otherwise the transcript with the most
rows for that gene. The genome build is read from `NCBI_Build`.

## HGVSc parsing

The leading `c.` integer of `HGVSc` is the variant anchor: positive
integers map to CDS positions, `-N` maps to the 5' UTR (negative axis),
`*N` maps to the 3' UTR (axis position `cds_len + N`). Strings with `+K`
or `-K` immediately after the leading anchor (splice / deep intronic)
are not plotted and counted under "splice"; anything else that fails to
parse is counted under "unparsed". The caption lists both counts.

## Ensembl source

By default the function calls `<host>/lookup/id/<ENST>?expand=1;utr=1`
on `rest.ensembl.org` (GRCh38) or `grch37.rest.ensembl.org` (GRCh37) and
caches the lean parsed structure as
`<cache_dir>/genestruct_ensembl_<ENST>_<assembly>.rds`. The cache
directory follows the same resolution chain as
[`gvr_lollipop()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_lollipop.md)
(explicit arg -\> env `GVR_CACHE_DIR` -\> option `germlinevaR.cache_dir`
-\> `tools::R_user_dir("germlinevaR","cache")` -\>
[`tempdir()`](https://rdrr.io/r/base/tempfile.html)). Set `gtf_path` to
a Gencode/Ensembl GTF to skip REST entirely; this requires suggesting
`rtracklayer` or falls back to a tiny streaming parser.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ## Load the shipped example table; use gtf_path=NULL + cache_dir=FALSE
  ## to skip the Ensembl REST call in the example (requires network).
  ## A real call would be: gvr_genepos.plot(gvr, "BRCA1")
  gvr <- readRDS(system.file("extdata", "example_gvr.rds",
                             package = "germlinevaR"))
  ## Confirm the function is available
  is.function(gvr_genepos.plot)
}
#> [1] TRUE

if (FALSE) { # \dontrun{
# Auto-resolve MANE/CANONICAL transcript for BRCA1
p <- gvr_genepos.plot(gvr, "BRCA1")

# Pin transcript and use proportional intron scaling
gvr_genepos.plot(gvr, "BRCA1",
                 transcript_id = "ENST00000357654",
                 intron_scale  = "proportional")

# Fully offline using a local GTF
gvr_genepos.plot(gvr, "BRCA1",
                 gtf_path = "gencode.v44.annotation.gtf.gz")
} # }
```
