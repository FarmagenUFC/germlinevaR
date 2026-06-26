# Generate a Vector of Colors From a Named Palette

`gvr_color_palette()` returns a character vector of `n` hex codes drawn
from one of 52 named palettes covering ColorBrewer, viridis, scientific
journals, publishers, corporate brands, colorblind-friendly schemes,
gradients, and classic R palettes. This is the helper that
[`gvr_lollipop()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_lollipop.md)
uses internally to resolve the `variant_palette` and `domain_palette`
arguments, but it is also exported so users can build their own color
vectors for any plot.

## Usage

``` r
gvr_color_palette(palette = "okabe_ito", n = 8L)
```

## Arguments

- palette:

  Character string. Name of the palette to use. See
  [`gvr_list_palettes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_list_palettes.md)
  for the full catalog of 52 names. Defaults to `"okabe_ito"`, the
  colorblind-friendly palette recommended by Okabe & Ito (2008).

- n:

  Integer. Number of colors to return. When `n` exceeds the base length
  of a discrete palette, intermediate colors are interpolated via
  [`grDevices::colorRampPalette()`](https://rdrr.io/r/grDevices/colorRamp.html).

## Value

A character vector of length `n` containing hex color codes.

## Details

The 52 palettes are grouped into the following families. Names are
identical to those exposed by the sibling `plot_multiple_compounds()`
function so the two packages share the same color vocabulary.

*Base:* `"hue"`, `"ggplot2"`, `"default"` (the ggplot2 default).

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

**Optional dependencies:** ColorBrewer palettes use RColorBrewer when
available and fall back to interpolated base colors otherwise. Viridis
palettes use viridisLite when available and fall back to manual hex
sequences otherwise. Both packages are in `Suggests:`.

## See also

[`gvr_list_palettes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_list_palettes.md),
[`gvr_lollipop()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_lollipop.md).

## Examples

``` r
# Get 8 Nature-journal colors
gvr_color_palette("nature", 8)
#> [1] "#0072B2" "#986332" "#797931" "#1D987A" "#AE7E9F" "#E0B66D" "#C4D671"
#> [8] "#56B4E9"

# Default Okabe-Ito palette for 6 categories
gvr_color_palette(n = 6)
#> [1] "#E69F00" "#56B4E9" "#009E73" "#F0E442" "#0072B2" "#D55E00"

# List all available palette names
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

# \donttest{
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    # Use directly in a ggplot2 plot
    ggplot2::ggplot(iris,
        ggplot2::aes(Sepal.Length, Sepal.Width, color = Species)) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(
        values = gvr_color_palette("viridis", 3))
  }

# }
```
