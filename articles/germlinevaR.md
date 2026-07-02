# germlinevaR walkthrough

## Introduction

**germlinevaR** is a self-contained R toolchain for single-sample
germline VCFs annotated with Ensembl VEP, SnpEff, or both. It turns
on-disk VCFs into a tabular variant `data.table`, filters them with a
modular set of toggles, isolates candidate novel variants, summarises a
cohort, and renders publication-ready per-gene plots — all inside one R
session, without an external pipeline.

### Motivation for Bioconductor

Germline variant analysis sits at the intersection of clinical genetics
and bioinformatics: researchers need to go from annotated VCFs to
filtered, interpretable variant tables and publication-ready figures
without stitching together multiple unrelated tools. germlinevaR fills
this gap by providing a single, cohesive workflow that handles the full
journey from VCF to cohort summary and per-gene visualisation, with
built-in support for the two most widely used annotators (Ensembl VEP
and SnpEff) and optional Brazilian population frequency annotation
(ABraOM SABE-609). Inclusion in Bioconductor makes the package
discoverable to the clinical genomics community and ensures long-term
maintenance through the Bioconductor release cycle.

### Comparison with related packages

Several Bioconductor and CRAN packages address parts of the germline
variant analysis workflow:

- **`VariantAnnotation`** (Bioconductor): provides low-level VCF parsing
  and annotation infrastructure. germlinevaR builds on top of this
  conceptual layer but targets a higher-level workflow: it reads
  VEP/SnpEff-annotated VCFs directly, applies opinionated filtering, and
  produces cohort summaries and plots without requiring users to
  manipulate `VCF` or `GRanges` objects.
- **`maftools`** (Bioconductor): a comprehensive somatic mutation
  analysis toolkit centred on the MAF format. germlinevaR is
  complementary: it targets the **germline** single-sample case,
  produces a tabular variant table with MAF-compatible core columns for
  interoperability with maftools, and adds VEP/SnpEff multi-annotator
  support and population-frequency filtering that maftools does not
  provide.
- **`vcfR`** (CRAN): focuses on VCF import and manipulation. germlinevaR
  goes further downstream: after reading the VCF it applies
  clinical-significance and population-frequency filters, generates a
  cohort summary, and renders per-gene lollipop plots.

germlinevaR’s unique combination is the auto-routing multi-annotator
reader + modular clinical filter + cohort summary + per-gene
visualisation, all in one package with no external pipeline dependency.

### Vignette scope

This vignette walks through the package on the **tiny example VCF
shipped with the package**: a single anonymised sample with 62 variants
covering all 17 VEP `Variant_Classification` levels, drawn from a real
exome to keep the example faithful. The whole walkthrough runs in a few
seconds.

What this vignette does **not** cover: variant calling, VCF annotation
(do it upstream with VEP or SnpEff), and somatic tumour/normal
pipelines. germlinevaR is built for germline single-sample analyses.

## Installation

``` r

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("germlinevaR")
```

For the development version from GitHub:

``` r

# install.packages("remotes")
remotes::install_github("FarmagenUFC/germlinevaR")
```

## Quick start

``` r

library(germlinevaR)
library(data.table)
#> 
#> Attaching package: 'data.table'
#> The following object is masked from 'package:base':
#> 
#>     %notin%
```

## Reading a VCF into a tabular variant table

[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
is the auto-routing entry point: it inspects the VCF header, detects
whether annotations come from VEP, SnpEff, or both, and dispatches to
the correct sibling reader (`read.gvr`, `read.gvr.snpeff`, or
`read.gvr.dual`). When you read your own data, point
[`read.gvr()`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
at the folder containing your VCF file(s) — the function discovers every
`*.vcf.gz` in that directory, auto-detects the annotator on each file,
and returns one tabular variant row per ALT allele. All files in a batch
must share the same annotator; mixed VEP/SnpEff folders are rejected.
For this vignette we stage the shipped example into a temporary
directory and read it from there:

``` r

vcf_dir <- tempfile("gvr_vig_")
dir.create(vcf_dir)
file.copy(
    system.file("extdata", "example.vep.vcf.gz", package = "germlinevaR"),
    file.path(vcf_dir, "example.vep.vcf.gz")
)
#> [1] TRUE

gvr <- read.gvr(vcf_dir, verbose = FALSE)
#> Warning: read.gvr: ABraOM reference unreadable; 'ABraOM_AF' left blank.
dim(gvr)
#> [1]  62 116
```

The returned `data.table` carries one row per ALT allele on the
most-severe transcript, with a canonical 116-column schema:

``` r

head(gvr[, .(Hugo_Symbol, Variant_Classification, IMPACT,
    Tumor_Sample_Barcode)])
#>    Hugo_Symbol Variant_Classification   IMPACT Tumor_Sample_Barcode
#>         <char>                 <char>   <char>               <char>
#> 1:       OR4F5                 Silent      LOW            Sample_01
#> 2:       OR4F5      Missense_Mutation MODERATE            Sample_01
#> 3:       OR4F5                 Silent      LOW            Sample_01
#> 4:     Unknown                 Intron MODIFIER            Sample_01
#> 5:      WASH9P                    RNA MODIFIER            Sample_01
#> 6:      SAMD11                  5'UTR MODIFIER            Sample_01
```

## Filtering: `gvr_filter()`

[`gvr_filter()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md)
exposes population-frequency, clinical-significance, biotype,
variant-classification, genotype, and panel filters as independent
toggles. Defaults drop variants common in gnomAD exomes, 1000 Genomes,
or ABraOM (above 1 %), and require either a pathogenic ClinVar
annotation or no annotation at all:

``` r

filt <- gvr_filter(gvr, verbose = FALSE)
dim(filt)
#> [1]   7 116
table(filt$Variant_Classification)
#> 
#>        Frame_Shift_Del        Frame_Shift_Ins                 Intron 
#>                      1                      1                      1 
#>      Missense_Mutation       Nonstop_Mutation            Splice_Site 
#>                      1                      1                      1 
#> Translation_Start_Site 
#>                      1
```

Defaults take the 62-row example down to 7 rows, one per surviving
variant classification. To restrict to a curated gene panel, pass
`genes =`:

``` r

gvr_list_panels()                    # registered panels
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
panel <- gvr_panel_genes("breast cancer")
panel
#>  [1] "AKT1"   "ATM"    "BARD1"  "BRCA1"  "BRCA2"  "BRIP1"  "CDH1"   "CHEK2" 
#>  [9] "EPCAM"  "ERBB2"  "ESR1"   "MLH1"   "MSH2"   "MSH6"   "NBN"    "NF1"   
#> [17] "PALB2"  "PIK3CA" "PMS2"   "PTEN"   "RAD51C" "RAD51D" "STK11"  "TP53"

filt_bc <- gvr_filter(gvr, genes = panel, verbose = FALSE)
filt_bc[, .(Hugo_Symbol, HGVSp_Short, Variant_Classification, IMPACT)]
#>    Hugo_Symbol HGVSp_Short Variant_Classification   IMPACT
#>         <char>      <char>                 <char>   <char>
#> 1:       BRCA2    p.G1376A      Missense_Mutation MODERATE
```

The one BC-panel survivor — BRCA2 p.G1376A — is a missense with no rsID
and no allele frequency in any catalogue, which is exactly the profile
that motivates the next section.

## Candidate novel variants: `gvr_novel()`

[`gvr_novel()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_novel.md)
is a dedicated subsetter that returns rows with no dbSNP rsID **and** no
allele frequency in gnomAD, 1000 Genomes, or ABraOM. These are the
variants most worth following up:

``` r

novel <- gvr_novel(gvr, verbose = FALSE)
dim(novel)
#> [1]   3 116
novel[, .(Hugo_Symbol, HGVSp_Short, Variant_Classification, IMPACT)]
#>    Hugo_Symbol HGVSp_Short Variant_Classification   IMPACT
#>         <char>      <char>                 <char>   <char>
#> 1:      ATAD3B                             Intron MODIFIER
#> 2:    C11orf21                        Splice_Site     HIGH
#> 3:       BRCA2    p.G1376A      Missense_Mutation MODERATE
```

Note that the BRCA2 p.G1376A row from the panel filter above also
appears here — a single variant can satisfy several criteria at once.

## Cohort summary: `gvr_summary()`

[`gvr_summary()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md)
returns an 8-section list-of-`data.table`s ready for direct inspection
or downstream tabulation. It can also write the same content as an Excel
workbook, a multi-page PDF, and a self-contained interactive HTML
dashboard with plotly drill-downs and DT tables; we disable those side
effects here to keep the vignette fast:

``` r

summ <- gvr_summary(
    gvr,
    save_excel = FALSE, save_pdf = FALSE, save_html = FALSE,
    verbose    = FALSE
)
names(summ)
#> [1] "overview"               "top_genes"              "top_genes_per_sample"  
#> [4] "variant_classification" "variant_type"           "clin_sig"              
#> [7] "impact"                 "top_variants"
summ$top_genes
#>     Hugo_Symbol Sample_01 Total
#>          <char>     <int> <num>
#>  1:       OR4F5         3     3
#>  2:      ACTRT2         1     1
#>  3:     ARHGEF7         1     1
#>  4:      ATAD3A         1     1
#>  5:      ATAD3B         1     1
#>  6:    ATP6V1B1         1     1
#>  7:       ATRIP         1     1
#>  8:       BRCA1         1     1
#>  9:       BRCA2         1     1
#> 10:        BRD7         1     1
#> 11:    C11orf21         1     1
#> 12:    C1QTNF12         1     1
#> 13:         CA6         1     1
#> 14:      CACHD1         1     1
#> 15:       CCNL2         1     1
#> 16:        CD1C         1     1
#> 17:      COL5A1         1     1
#> 18:    CSNKA2IP         1     1
#> 19:        EFL1         1     1
#> 20:        ESPN         1     1
#>     Hugo_Symbol Sample_01 Total
#>          <char>     <int> <num>
```

For the full HTML dashboard call
`gvr_summary(..., save_html = TRUE, out_dir = ".")`.

## Static plots

germlinevaR produces several families of plots: per-gene
**protein-domain lollipops** (`gvr_lollipop`), per-gene **gene-structure
lollipops** (`gvr_genepos.plot`), a `ComplexHeatmap`-based **top-genes
variant matrix** (`gvr_plot`), and the **HTML-dashboard panels** that
`gvr_summary(save_html = TRUE)` embeds in its interactive output (also
writable as standalone PNGs via
[`gvr_sum_plots()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_sum_plots.md)).
The figures below were rendered from a real cohort — the tiny example
tabular variant table shipped with the package is too sparse for a
meaningful top-genes variant matrix or dashboard, so those two families
are described in the text only; call the functions yourself against your
own tabular variant table to reproduce the layouts.

**Protein-domain lollipop** — `gvr_lollipop(gvr, gene = "TP53")` draws
variants along the protein, with InterPro domains fetched and cached on
the fly:

![Lollipop for TP53 from gvr_lollipop() — InterPro domains
shaded](figures/lollipop_TP53.png)

Lollipop for TP53 from gvr_lollipop() — InterPro domains shaded

**Gene-structure lollipop** — `gvr_genepos.plot(gvr, gene = "BRCA1")`
draws variants along the gene’s cDNA, with exon/intron/UTR regions
resolved via Ensembl REST or a local GTF:

![Gene-structure track for BRCA1 from gvr_genepos.plot() — exons/UTRs
from Ensembl REST](figures/genepos_BRCA1.png)

Gene-structure track for BRCA1 from gvr_genepos.plot() — exons/UTRs from
Ensembl REST

**Top-genes variant matrix** —
`gvr_plot(gvr, top_n = 20, out_dir = ".")` writes a `ComplexHeatmap`
top-genes variant matrix PNG of the top mutated genes across samples,
coloured by the most severe variant class per cell. The top bar shows
per-sample variant impact (HIGH / MODERATE / LOW / MODIFIER); the right
bar shows per-gene total burden. The figure below uses illustrative
multi-sample data to demonstrate the layout:

![Top-genes variant matrix from gvr_plot() — illustrative multi-sample
data (8 samples, 20 cancer genes). In a real cohort, replace with your
own gvr table.](figures/gvr_plot_demo.png)

Top-genes variant matrix from gvr_plot() — illustrative multi-sample
data (8 samples, 20 cancer genes). In a real cohort, replace with your
own gvr table.

**HTML cohort dashboard** —
`gvr_summary(gvr, save_html = TRUE, out_dir = ".")` writes a
self-contained interactive HTML dashboard with plotly drill-downs and DT
tables. It opens with four KPI cards (total variants, samples, distinct
genes, HIGH-impact count) followed by three bar charts (top mutated
genes, variant classification, IMPACT) and a top-variants table. The
figure below is a screenshot of a dashboard rendered from a **synthetic
8-sample cohort** fabricated purely to demonstrate the layout; do not
read biological meaning into any specific value.

![Screenshot of gvr_summary(save_html = TRUE) output — synthetic
8-sample cohort rendered for layout demonstration only, not real
data.](figures/gvr_summary_dashboard.png)

Screenshot of gvr_summary(save_html = TRUE) output — synthetic 8-sample
cohort rendered for layout demonstration only, not real data.

[`gvr_sum_plots()`](https://farmagenufc.github.io/germlinevaR/reference/gvr_sum_plots.md)
writes the same panels that `gvr_summary(save_html = TRUE)` embeds in
its interactive dashboard (top genes, variant classification, IMPACT,
top variants) as standalone PNGs.

The demo cohort shown above is fabricated by
`inst/scripts/build_synthetic_dashboard.R` (also reachable via
`system.file("scripts", "build_synthetic_dashboard.R", package = "germlinevaR")`).
Source that script in a fresh R session to reproduce the exact
106-variant / 19-gene / 8-sample `gvr` table used here.

## Where to next

- [`?read.gvr`](https://farmagenufc.github.io/germlinevaR/reference/read.gvr.md)
  — every reader option (parallel cohort folders, ABraOM joining, custom
  output paths)
- [`?gvr_filter`](https://farmagenufc.github.io/germlinevaR/reference/gvr_filter.md)
  — full toggle reference for each filter axis
- [`?gvr_summary`](https://farmagenufc.github.io/germlinevaR/reference/gvr_summary.md)
  — multi-section cohort summary with optional Excel, PDF, and HTML
  outputs
- [`?gvr_lollipop`](https://farmagenufc.github.io/germlinevaR/reference/gvr_lollipop.md)
  and
  [`?gvr_genepos.plot`](https://farmagenufc.github.io/germlinevaR/reference/gvr_genepos.plot.md)
  — per-gene plotting
- GitHub: <https://github.com/FarmagenUFC/germlinevaR>

## Session info

    #> R version 4.6.1 (2026-06-24)
    #> Platform: x86_64-pc-linux-gnu
    #> Running under: Ubuntu 24.04.4 LTS
    #> 
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    #> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    #> 
    #> locale:
    #>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    #>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    #>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    #> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    #> 
    #> time zone: UTC
    #> tzcode source: system (glibc)
    #> 
    #> attached base packages:
    #> [1] stats     graphics  grDevices utils     datasets  methods   base     
    #> 
    #> other attached packages:
    #> [1] data.table_1.18.4  germlinevaR_0.99.3 BiocStyle_2.40.0  
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] sass_0.4.10           generics_0.1.4        shape_1.4.6.1        
    #>  [4] stringi_1.8.7         magrittr_2.0.5        digest_0.6.39        
    #>  [7] evaluate_1.0.5        grid_4.6.1            RColorBrewer_1.1-3   
    #> [10] bookdown_0.47         iterators_1.0.14      circlize_0.4.18      
    #> [13] fastmap_1.2.0         foreach_1.5.2         doParallel_1.0.17    
    #> [16] jsonlite_2.0.0        zip_3.0.0             GlobalOptions_0.1.4  
    #> [19] BiocManager_1.30.27   ComplexHeatmap_2.28.0 scales_1.4.0         
    #> [22] codetools_0.2-20      textshaping_1.0.5     jquerylib_0.1.4      
    #> [25] cli_3.6.6             rlang_1.2.0           crayon_1.5.3         
    #> [28] cachem_1.1.0          yaml_2.3.12           otel_0.2.0           
    #> [31] tools_4.6.1           parallel_4.6.1        dplyr_1.2.1          
    #> [34] colorspace_2.1-2      ggplot2_4.0.3         GetoptLong_1.1.1     
    #> [37] BiocGenerics_0.58.1   vctrs_0.7.3           R6_2.6.1             
    #> [40] png_0.1-9             matrixStats_1.5.0     stats4_4.6.1         
    #> [43] lifecycle_1.0.5       S4Vectors_0.50.1      fs_2.1.0             
    #> [46] htmlwidgets_1.6.4     IRanges_2.46.0        clue_0.3-68          
    #> [49] ragg_1.5.2            cluster_2.1.8.2       pkgconfig_2.0.3      
    #> [52] desc_1.4.3            openxlsx_4.2.8.1      pillar_1.11.1        
    #> [55] pkgdown_2.2.0         bslib_0.11.0          gtable_0.3.6         
    #> [58] Rcpp_1.1.1-1.1        glue_1.8.1            systemfonts_1.3.2    
    #> [61] tidyselect_1.2.1      tibble_3.3.1          xfun_0.59            
    #> [64] knitr_1.51            farver_2.1.2          rjson_0.2.23         
    #> [67] htmltools_0.5.9       rmarkdown_2.31        compiler_4.6.1       
    #> [70] S7_0.2.2
