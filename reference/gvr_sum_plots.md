# Export gvr_summary() plots as standalone image files

Produces the same cohort-level and per-sample plots that
[`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md)
renders in its HTML/PDF dashboard, but as individual image files and as
multi-panel composites, written into a new folder.

## Usage

``` r
gvr_sum_plots(
  gvr,
  out_dir = ".",
  folder_name = "gvr_sum_plots",
  format = "png",
  width = 7,
  height = 5,
  dpi = 300,
  sample_col = "Tumor_Sample_Barcode",
  top_n_genes = 20,
  top_n_variants = 20,
  per_sample = TRUE,
  panel = TRUE,
  verbose = TRUE
)
```

## Arguments

- gvr:

  data.table or data.frame produced by
  [`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
  /
  [`read.gvr.snpeff()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.snpeff.md)
  /
  [`read.gvr.dual()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.dual.md).
  Required columns: `Hugo_Symbol`, `Variant_Classification`,
  `Variant_Type`, `IMPACT`, `CLIN_SIG`; `dbSNP_RS` is optional.

- out_dir:

  parent directory under which the output folder is created. Default
  `"."` (current working directory).

- folder_name:

  name of the new folder created under `out_dir`. Default
  `"gvr_sum_plots"`.

- format:

  one of `"png"` (default), `"pdf"`, `"svg"`, `"jpeg"`, `"tiff"`,
  `"bmp"`, `"eps"`, `"ps"`, `"tex"`, `"wmf"`. Mirrors the device list of
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

- width:

  plot width in inches (passed to
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)).
  Default 7.

- height:

  plot height in inches. Default 5.

- dpi:

  resolution for raster formats (ignored for vector). Default 300.

- sample_col:

  column holding the sample identifier. Default
  `"Tumor_Sample_Barcode"`.

- top_n_genes:

  number of top genes to keep in `top_genes` and per-sample plots.
  Default 20.

- top_n_variants:

  number of top variants to keep in `top_variants`. Default 20.

- per_sample:

  logical; if `FALSE`, skip per-sample plots and the per-sample panel.
  Default `TRUE`.

- panel:

  logical; if `FALSE`, skip both panel images. Default `TRUE`.

- verbose:

  logical; print progress messages. Default `TRUE`.

## Value

The output folder path, invisibly.

## Details

The function recomputes the summary sections from the table (it does not
depend on a prior call to
[`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md)),
then writes:

Cohort-level plots (always produced):

- `top_genes.<ext>`: top genes by variant count

- `variant_classification.<ext>`: top variant classifications

- `impact.<ext>`: VEP IMPACT severity

- `top_variants.<ext>`: top variants by `dbSNP_RS` recurrence (omitted
  if the column is missing or all empty)

Per-sample plots (when `per_sample = TRUE`, one file per sample):

- `per_sample/top_genes__<sample>.<ext>`

Panels (when `panel = TRUE`):

- `panel_cohort.<ext>`: 2x2 grid of the four cohort plots

- `panel_per_sample.<ext>`: grid of every per-sample plot (omitted if
  `n_samples == 1` since it would duplicate the standalone per-sample
  image)

All files are written under `out_dir/folder_name/`. Existing files of
the same name are overwritten; other files in the folder are left alone.

## See also

[`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md)
for the dashboard view,
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
for the reader.

## Author

Thiago Loreto Matos

## Examples

``` r
## Load the shipped example table and write plots to a temp directory
gvr <- readRDS(system.file("extdata", "example_gvr.rds",
                           package = "germlinevaR"))
out_dir <- gvr_sum_plots(gvr, out_dir = tempdir(), verbose = FALSE)
dir.exists(out_dir)
#> [1] TRUE

## Write PDF outputs into a named subfolder under tempdir()
gvr_sum_plots(gvr, out_dir = tempdir(), folder_name = "S6_plots",
              format = "pdf", verbose = FALSE)
```
