# Genes for a Disease Panel

Returns the Hugo_Symbol vector associated with a given disease panel
name. The same registry is used by
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
when the user passes a `panel` argument.

## Usage

``` r
gvr_panel_genes(panel)
```

## Arguments

- panel:

  Single character string naming a disease panel (case-insensitive;
  underscores are treated as spaces, e.g. `"breast_cancer"` is
  equivalent to `"breast cancer"`; selected aliases like
  `"gastrointestinal stromal tumor"` -\> `"gist"` are recognised). Must
  resolve to exactly one entry in the registry.

## Value

Character vector of Hugo_Symbols (uppercase, sorted, unique).

## Details

If `panel` does not resolve to a known disease, an error is thrown
listing all available panels with their gene counts, so the caller sees
the catalog at the point of failure.

## See also

[`gvr_list_panels()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_list_panels.md),
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md).

Other germlinevaR:
[`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md),
[`gvr_hpo_genes()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_hpo_genes.md),
[`gvr_list_panels()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_list_panels.md),
[`gvr_plot()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_plot.md),
[`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md),
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md),
[`read.gvr.snpeff()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.snpeff.md)

## Author

germlinevaR authors

## Examples

``` r
gvr_panel_genes("breast cancer")
#>  [1] "AKT1"   "ATM"    "BARD1"  "BRCA1"  "BRCA2"  "BRIP1"  "CDH1"   "CHEK2" 
#>  [9] "EPCAM"  "ERBB2"  "ESR1"   "MLH1"   "MSH2"   "MSH6"   "NBN"    "NF1"   
#> [17] "PALB2"  "PIK3CA" "PMS2"   "PTEN"   "RAD51C" "RAD51D" "STK11"  "TP53"  
gvr_panel_genes("Breast_Cancer")            # case + underscore alias
#>  [1] "AKT1"   "ATM"    "BARD1"  "BRCA1"  "BRCA2"  "BRIP1"  "CDH1"   "CHEK2" 
#>  [9] "EPCAM"  "ERBB2"  "ESR1"   "MLH1"   "MSH2"   "MSH6"   "NBN"    "NF1"   
#> [17] "PALB2"  "PIK3CA" "PMS2"   "PTEN"   "RAD51C" "RAD51D" "STK11"  "TP53"  
gvr_panel_genes("gastrointestinal stromal tumor")  # alias of "gist"
#>  [1] "BRAF"   "KIT"    "NF1"    "NTRK1"  "NTRK2"  "NTRK3"  "PDGFRA" "SDHA"  
#>  [9] "SDHB"   "SDHC"   "SDHD"  
gvr_panel_genes("hereditary prostate cancer")
#>  [1] "ATM"    "ATR"    "BRCA1"  "BRCA2"  "BRIP1"  "CHEK2"  "EPCAM"  "HOXB13"
#>  [9] "MLH1"   "MRE11"  "MSH2"   "MSH6"   "NBN"    "PALB2"  "PMS2"   "PTEN"  
#> [17] "RAD51C" "RAD51D" "TP53"  

## Combine multiple panels and post-hoc filter the pre-parsed example
## table (equivalent to read.gvr(..., panel = c(...)) but instantaneous):
gvr <- readRDS(system.file("extdata", "example_gvr.rds",
    package = "germlinevaR"))
multi_panel <- unique(c(gvr_panel_genes("breast cancer"),
    gvr_panel_genes("hereditary prostate cancer")))
nrow(gvr[gvr$Hugo_Symbol %in% multi_panel, ])
#> [1] 2
```
