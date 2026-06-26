# List All Available Palettes

Returns the sorted character vector of every palette name accepted by
[`gvr_color_palette()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_color_palette.md).
Use this for documentation, autocomplete, or programmatic iteration over
palettes.

## Usage

``` r
gvr_list_palettes()
```

## Value

Character vector of palette names.

## Details

The 56 palette names returned by this function are organized into the
following 11 families. The same taxonomy is used by
[`gvr_color_palette()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_color_palette.md)
and the sibling `plot_multiple_compounds()` plotting function so the two
packages share the same color vocabulary.

*Base:* `"hue"`, `"ggplot2"`.

*ColorBrewer qualitative (good for distinct categories):* `"set1"`,
`"set2"`, `"set3"`, `"dark2"`, `"paired"`, `"accent"`, `"pastel1"`,
`"pastel2"`.

*ColorBrewer sequential (ordered data):* `"blues"`, `"reds"`,
`"greens"`, `"purples"`, `"oranges"`, `"greys"`.

*ColorBrewer diverging:* `"spectral"`, `"rdylbu"`, `"rdylgn"`, `"piyg"`,
`"prgn"`.

*Viridis (perceptually uniform):* `"viridis"`, `"magma"`, `"inferno"`,
`"plasma"`.

*Journal palettes:* `"nature"`, `"science"`, `"cell"`, `"plos"`,
`"elife"`.

*Publisher palettes:* `"bmc"`, `"frontiers"`, `"wiley"`, `"elsevier"`,
`"oxford"`, `"springer"`, `"acs"`, `"rsc"`.

*Corporate brand palettes:* `"ibm"`, `"google"`, `"microsoft"`,
`"twitter"`.

*Colorblind-friendly:* `"okabe_ito"`, `"colorblind"`, `"cud"`, `"tol"`.

*Gradient palettes:* `"blue_red"`, `"green_red"`, `"purple_orange"`,
`"cool_warm"`, `"blue_yellow"`.

*Classic R:* `"rainbow"`, `"heat"`, `"terrain"`, `"topo"`, `"cm"`.

See
[`gvr_color_palette()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_color_palette.md)
for example hex codes from each family.

## See also

[`gvr_color_palette()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_color_palette.md),
[`gvr_lollipop()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_lollipop.md).

## Examples

``` r
gvr_list_palettes()
#>  [1] "hue"           "ggplot2"       "set1"          "set2"         
#>  [5] "set3"          "dark2"         "paired"        "accent"       
#>  [9] "pastel1"       "pastel2"       "blues"         "reds"         
#> [13] "greens"        "purples"       "oranges"       "greys"        
#> [17] "spectral"      "rdylbu"        "rdylgn"        "piyg"         
#> [21] "prgn"          "viridis"       "magma"         "inferno"      
#> [25] "plasma"        "nature"        "science"       "cell"         
#> [29] "plos"          "elife"         "bmc"           "frontiers"    
#> [33] "wiley"         "elsevier"      "oxford"        "springer"     
#> [37] "acs"           "rsc"           "ibm"           "google"       
#> [41] "microsoft"     "twitter"       "okabe_ito"     "colorblind"   
#> [45] "cud"           "tol"           "blue_red"      "green_red"    
#> [49] "purple_orange" "cool_warm"     "blue_yellow"   "rainbow"      
#> [53] "heat"          "terrain"       "topo"          "cm"           
# Generate a 5-color sample of every palette
if (FALSE) { # \dontrun{
sapply(gvr_list_palettes(), gvr_color_palette, n = 5)
} # }
```
