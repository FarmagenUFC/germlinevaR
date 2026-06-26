# List Available Disease Gene Panels

Returns the sorted character vector of every disease name accepted by
[`gvr_panel_genes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_panel_genes.md)
and by
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)'s
`panel` argument.

## Usage

``` r
gvr_list_panels()
```

## Value

Character vector of canonical panel names (lowercase, spaces).

## Details

Each canonical name corresponds to a curated Hugo_Symbol gene list.
Panel names are matched case-insensitively and underscores are treated
as spaces. A small alias table is also recognised: e.g.
`"gastrointestinal stromal tumor"` is accepted as a synonym for
`"gist"`, and the common typo `"pheocromocytoma"` resolves to
`"pheochromocytoma"`.

Panels currently shipped:

- `"men1"`:

  13 genes.

- `"acromegaly"`:

  23 genes.

- `"pheochromocytoma"`:

  24 genes. (alias: `"pheocromocytoma"`)

- `"hereditary cancer"`:

  87 genes.

- `"gist"`:

  11 genes. (alias: `"gastrointestinal stromal tumor"`)

- `"lynch syndrome"`:

  5 genes.

- `"li-fraumeni syndrome"`:

  4 genes.

- `"hereditary gastric cancer"`:

  28 genes.

- `"hereditary colorectal cancer"`:

  23 genes.

- `"familial adenomatous polyposis"`:

  19 genes.

- `"hereditary melanoma cancer"`:

  20 genes.

- `"hereditary prostate cancer"`:

  19 genes.

- `"hereditary breast and ovarian cancer"`:

  25 genes.

- `"breast cancer"`:

  24 genes.

- `"breast cancer somatic"`:

  6 genes. (alias: `"breast cancer somatic panel"`)

Use
[`gvr_panel_genes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_panel_genes.md)
to retrieve the gene vector for a specific disease, or pass the panel
name (or vector of names) to
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
via its `panel` argument to filter a gvr table.

## See also

[`gvr_panel_genes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_panel_genes.md),
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md).

Other germlinevaR:
[`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md),
[`gvr_panel_genes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_panel_genes.md),
[`gvr_plot()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_plot.md),
[`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md),
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md),
[`read.gvr.snpeff()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.snpeff.md)

## Author

germlinevaR authors

## Examples

``` r
gvr_list_panels()
#>  [1] "acromegaly"                          
#>  [2] "breast cancer"                       
#>  [3] "breast cancer somatic"               
#>  [4] "familial adenomatous polyposis"      
#>  [5] "gist"                                
#>  [6] "hereditary breast and ovarian cancer"
#>  [7] "hereditary cancer"                   
#>  [8] "hereditary colorectal cancer"        
#>  [9] "hereditary gastric cancer"           
#> [10] "hereditary melanoma cancer"          
#> [11] "hereditary prostate cancer"          
#> [12] "li-fraumeni syndrome"                
#> [13] "lynch syndrome"                      
#> [14] "men1"                                
#> [15] "pheochromocytoma"                    

## Use in read.gvr() to keep only breast-cancer genes
if (FALSE) { # \dontrun{
gvr <- read.gvr("/path/to/folder", panel = "breast cancer")
} # }
```
