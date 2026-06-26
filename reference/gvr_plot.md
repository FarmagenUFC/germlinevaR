# Cohort top-genes variant matrix from a germline gvr table (read.gvr / gvr_filter output)

Draws a top-genes variant matrix from an MAF-like table - the output of
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md),
or of
[`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md) -
and writes it to a PNG file. Rows are the top-`top_n` genes (ranked by
number of distinct samples mutated, then by variant count); columns are
samples. Each gene x sample cell shows the single MOST-SEVERE
`Variant_Classification` observed for that gene in that sample. The plot
is drawn with ComplexHeatmap.

## Usage

``` r
gvr_plot(
  gvr,
  top_n = 20,
  sample_col = "Tumor_Sample_Barcode",
  out_dir = ".",
  file_prefix = "gvr_plot",
  sample_name_rot = 45,
  gene_name_size = 10,
  sample_name_size = 10,
  axis_tick_size = 7,
  legend_label_size = 10,
  legend_title_size = 10,
  legend_label_wrap_chars = Inf,
  impact_title_side = c("left", "right"),
  verbose = TRUE
)
```

## Arguments

- gvr:

  An MAF-like `data.table`/`data.frame` from
  [`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
  or
  [`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md).
  Required columns: `Hugo_Symbol`, `Variant_Classification`.

- top_n:

  Integer; number of genes (rows) shown, ranked by number of distinct
  samples mutated then by variant count. Default `20`.

- sample_col:

  Name of the per-sample column. Default `"Tumor_Sample_Barcode"`. If
  absent, all rows are pooled into a single sample `"All"` (with a
  warning).

- out_dir:

  Output directory for the PNG. Created if it does not exist. Default
  `"."` (current working directory).

- file_prefix:

  Base filename for the written PNG. Default `"gvr_plot"`; the file is
  written as `<file_prefix>.png` (fixed name, no timestamp), e.g.
  `gvr_plot.png`. An existing file at that path is overwritten (a
  message is emitted when `verbose = TRUE`).

- sample_name_rot:

  Numeric; rotation angle (degrees) for the sample-name labels at the
  top of the heatmap. Default `45`. Common alternatives are `0`
  (horizontal) and `90` (vertical). Must be a single finite numeric.

- gene_name_size:

  Numeric; font size in points for the gene-name row labels (left side
  of the heatmap). Default `10`. Must be a single finite positive
  numeric.

- sample_name_size:

  Numeric; font size in points for the sample-name column labels (top of
  the heatmap). Default `10`. Must be a single finite positive numeric.

- axis_tick_size:

  Numeric; font size in points for the y-axis tick labels of the top
  "Variant impact" bar and the right gene-burden bar. Default `7`. Must
  be a single finite positive numeric.

- legend_label_size:

  Numeric; font size in points for the label text of both side legends
  ("Impact" and "Most severe class"). Default `10`. Must be a single
  finite positive numeric.

- legend_title_size:

  Numeric; font size in points for the title text of both side legends
  ("Impact" and "Most severe class"). The bold face is preserved.
  Default `10`. Must be a single finite positive numeric.

- legend_label_wrap_chars:

  Numeric; if a (prettified) legend label is longer than this many
  characters, it is wrapped onto two lines at the space closest to the
  middle. Default `Inf` disables wrapping (so all labels stay on one
  line, matching the legacy behaviour). Use a finite integer (e.g. `14`)
  when long labels would clip the right edge at large
  `legend_label_size`. Labels without internal spaces (e.g. `5'UTR`,
  `RNA`) are never wrapped regardless of length.

- impact_title_side:

  One of `"left"` (default) or `"right"`; controls where the
  `"Variant impact"` annotation title is drawn relative to the top
  stacked-bar panel. `"left"` renders the title vertically (acting as a
  y-axis title for the impact panel); `"right"` renders it horizontally
  on the right of the panel (the previous default).

- verbose:

  Logical; if `TRUE` (default) print the path of the file written.

## Value

Invisibly, the path of the written PNG (character), or `NA_character_`
if no known-gene variants are present in the table.

## Details

This is the standalone top-genes variant matrix previously produced
inside
[`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md).
It needs only the `Hugo_Symbol` and `Variant_Classification` columns
(plus the per-sample column).

Cell collapse: when a gene has several variant classes in one sample,
the cell is coloured by the most severe class, using this ranking (high
to low): Translation_Start_Site, Nonsense_Mutation, Nonstop_Mutation,
Splice_Site, Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del,
In_Frame_Ins, Missense_Mutation, Splice_Region,
Protein_altering_variant, Silent, 5'UTR, 3'UTR, 5'Flank, 3'Flank, RNA,
Intron, IGR, Targeted_Region. Any class outside this list ranks last and
is coloured grey. Colours follow a colourblind-safe (Okabe-Ito) palette.

Annotations: a right-side bar shows each gene's total variant burden; a
top bar shows each sample's total variant burden (axis labelled in
thousands). Empty cells (gene not mutated in that sample) are light
grey.

Data conventions:

- "Missing" means `NA` OR empty string `""`.

- Unknown/blank gene symbols are `Hugo_Symbol` in
  `c(".", "", "Unknown")`; these are excluded from the top-genes variant
  matrix.

- Works on ANY MAF-like table; it makes no assumption about prior
  filtering. It is commonly run on
  [`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md)
  output.

## Dependencies

Uses ComplexHeatmap (a Bioconductor package, listed in `Imports`).

## See also

[`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md)
for the tabular summary,
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
to build the table,
[`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md)
to filter it before plotting.

Other germlinevaR:
[`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md),
[`gvr_list_panels()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_list_panels.md),
[`gvr_panel_genes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_panel_genes.md),
[`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md),
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md),
[`read.gvr.snpeff()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.snpeff.md)

## Author

germlinevaR authors

## Examples

``` r
## Load the shipped example table; write plot to a temp directory
gvr <- readRDS(system.file("extdata", "example_gvr.rds",
                           package = "germlinevaR"))
p <- gvr_plot(gvr, out_dir = tempdir(), verbose = FALSE)
class(p)
#> [1] "character"

# \donttest{
  ## Smaller top-genes variant matrix of filtered hits to a temp folder
  gvr <- readRDS(system.file("extdata", "example_gvr.rds",
                             package = "germlinevaR"))
  filt <- gvr_filter(gvr, ABraOM_AF = NULL, verbose = FALSE)
  if (nrow(filt) > 0L) {
    gvr_plot(filt, top_n = 15, out_dir = tempdir(), verbose = FALSE)
  }
# }
```
