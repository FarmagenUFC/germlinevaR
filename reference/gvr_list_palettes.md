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
lapply(gvr_list_palettes(), gvr_color_palette, n = 5)
#> [[1]]
#> [1] "#F8766D" "#A3A500" "#00BF7D" "#00B0F6" "#E76BF3"
#> 
#> [[2]]
#> [1] "#F8766D" "#A3A500" "#00BF7D" "#00B0F6" "#E76BF3"
#> 
#> [[3]]
#> [1] "#E41A1C" "#4DAF4A" "#FF7F00" "#A65628" "#999999"
#> 
#> [[4]]
#> [1] "#66C2A5" "#A89BB0" "#C6B18B" "#F8D348" "#B3B3B3"
#> 
#> [[5]]
#> [1] "#8DD3C7" "#EB8E8B" "#D8C965" "#D1C2D2" "#FFED6F"
#> 
#> [[6]]
#> [1] "#1B9E77" "#8D6B86" "#A66753" "#D59D08" "#666666"
#> 
#> [[7]]
#> [1] "#A6CEE3" "#52AF43" "#F06C45" "#B294C7" "#B15928"
#> 
#> [[8]]
#> [1] "#7FC97F" "#EDBB99" "#9BB5A4" "#E31864" "#666666"
#> 
#> [[9]]
#> [1] "#FBB4AE" "#CCEBC5" "#FED9A6" "#E5D8BD" "#F2F2F2"
#> 
#> [[10]]
#> [1] "#B3E2CD" "#D7D3D9" "#EDDFD6" "#FBEDB5" "#CCCCCC"
#> 
#> [[11]]
#> [1] "#F7FBFF" "#C6DBEF" "#6BAED6" "#2171B5" "#08306B"
#> 
#> [[12]]
#> [1] "#FFF5F0" "#FCBBA1" "#FB6A4A" "#CB181D" "#67000D"
#> 
#> [[13]]
#> [1] "#F7FCF5" "#C7E9C0" "#74C476" "#238B45" "#00441B"
#> 
#> [[14]]
#> [1] "#FCFBFD" "#DADAEB" "#9E9AC8" "#6A51A3" "#3F007D"
#> 
#> [[15]]
#> [1] "#FFF5EB" "#FDD0A2" "#FD8D3C" "#D94801" "#7F2704"
#> 
#> [[16]]
#> [1] "#FFFFFF" "#D9D9D9" "#969696" "#525252" "#000000"
#> 
#> [[17]]
#> [1] "#9E0142" "#F88D51" "#FFFFBF" "#88CFA4" "#5E4FA2"
#> 
#> [[18]]
#> [1] "#A50026" "#F88D51" "#FFFFBF" "#8FC3DD" "#313695"
#> 
#> [[19]]
#> [1] "#A50026" "#F88D51" "#FFFFBF" "#86CB66" "#006837"
#> 
#> [[20]]
#> [1] "#8E0152" "#E796C3" "#F7F7F7" "#9BCE63" "#276419"
#> 
#> [[21]]
#> [1] "#40004B" "#AD8ABC" "#F7F7F7" "#80C480" "#00441B"
#> 
#> [[22]]
#> [1] "#440154FF" "#3B528BFF" "#21908CFF" "#5DC863FF" "#FDE725FF"
#> 
#> [[23]]
#> [1] "#000004FF" "#51127CFF" "#B63679FF" "#FB8861FF" "#FCFDBFFF"
#> 
#> [[24]]
#> [1] "#000004FF" "#56106EFF" "#BB3754FF" "#F98C0AFF" "#FCFFA4FF"
#> 
#> [[25]]
#> [1] "#0D0887FF" "#7E03A8FF" "#CC4678FF" "#F89441FF" "#F0F921FF"
#> 
#> [[26]]
#> [1] "#0072B2" "#D55E00" "#009E73" "#CC79A7" "#F0E442"
#> 
#> [[27]]
#> [1] "#3070B0" "#B03070" "#30B070" "#B07030" "#7030B0"
#> 
#> [[28]]
#> [1] "#DC143C" "#4682B4" "#2E8B57" "#FF8C00" "#9370DB"
#> 
#> [[29]]
#> [1] "#3498DB" "#E74C3C" "#2ECC71" "#F39C12" "#9B59B6"
#> 
#> [[30]]
#> [1] "#F04E4E" "#4EA5F0" "#4EF0A5" "#F0A54E" "#A54EF0"
#> 
#> [[31]]
#> [1] "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99"
#> 
#> [[32]]
#> [1] "#4DBBD5" "#E64B35" "#00A087" "#3C5488" "#F39B7F"
#> 
#> [[33]]
#> [1] "#5A9BD5" "#ED7D31" "#A5A5A5" "#FFC000" "#4472C4"
#> 
#> [[34]]
#> [1] "#F39800" "#DC143C" "#004080" "#009944" "#8B4513"
#> 
#> [[35]]
#> [1] "#002147" "#8B0000" "#006A4E" "#FF6B35" "#5E2A84"
#> 
#> [[36]]
#> [1] "#B22222" "#0066CC" "#228B22" "#FF8C00" "#9400D3"
#> 
#> [[37]]
#> [1] "#0066CC" "#CC0000" "#009966" "#FF9900" "#660099"
#> 
#> [[38]]
#> [1] "#B31B1B" "#005F8C" "#2E8B57" "#FF7F0E" "#9467BD"
#> 
#> [[39]]
#> [1] "#648FFF" "#785EF0" "#DC267F" "#FE6100" "#FFB000"
#> 
#> [[40]]
#> [1] "#4285F4" "#EA4335" "#FBBC05" "#34A853" "#FF6D00"
#> 
#> [[41]]
#> [1] "#F65314" "#7CBB00" "#00A1F1" "#FFBB00" "#A0A0A0"
#> 
#> [[42]]
#> [1] "#1DA1F2" "#14171A" "#657786" "#AAB8C2" "#E1E8ED"
#> 
#> [[43]]
#> [1] "#E69F00" "#56B4E9" "#009E73" "#F0E442" "#0072B2"
#> 
#> [[44]]
#> [1] "#0072B2" "#D55E00" "#009E73" "#CC79A7" "#F0E442"
#> 
#> [[45]]
#> [1] "#0072B2" "#E69F00" "#009E73" "#F0E442" "#56B4E9"
#> 
#> [[46]]
#> [1] "#4477AA" "#EE6677" "#228833" "#CCBB44" "#66CCEE"
#> 
#> [[47]]
#> [1] "#313695" "#8FC3DC" "#FFFFBF" "#F88D52" "#A50026"
#> 
#> [[48]]
#> [1] "#1A9850" "#A3D76A" "#EBE78B" "#FCA165" "#D73027"
#> 
#> [[49]]
#> [1] "#542788" "#A8A1CC" "#EBDCD0" "#F4B25D" "#B35806"
#> 
#> [[50]]
#> [1] "#3A5F8F" "#8EAECA" "#DBDDDB" "#EBB891" "#CC673B"
#> 
#> [[51]]
#> [1] "#053061" "#6AABD0" "#F7F7F7" "#E58267" "#67001F"
#> 
#> [[52]]
#> [1] "#FF0000" "#CCFF00" "#00FF66" "#0066FF" "#CC00FF"
#> 
#> [[53]]
#> [1] "#FF0000" "#FF5500" "#FFAA00" "#FFFF00" "#FFFF80"
#> 
#> [[54]]
#> [1] "#00A600" "#E6E600" "#EAB64E" "#EEB99F" "#F2F2F2"
#> 
#> [[55]]
#> [1] "#4C00FF" "#004CFF" "#00E5FF" "#00FF4D" "#FFFF00"
#> 
#> [[56]]
#> [1] "#80FFFF" "#BFFFFF" "#FFFFFF" "#FFBFFF" "#FF80FF"
#> 
```
