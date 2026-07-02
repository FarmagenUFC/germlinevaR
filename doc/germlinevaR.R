## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment  = "#>",
    fig.align = "center"
)

## ----install, eval=FALSE------------------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("germlinevaR")

## ----install-dev, eval=FALSE--------------------------------------------------
# # install.packages("remotes")
# remotes::install_github("FarmagenUFC/germlinevaR")

## ----setup--------------------------------------------------------------------
library(germlinevaR)
library(data.table)

## ----read-gvr-----------------------------------------------------------------
vcf_dir <- tempfile("gvr_vig_")
dir.create(vcf_dir)
file.copy(
    system.file("extdata", "example.vep.vcf.gz", package = "germlinevaR"),
    file.path(vcf_dir, "example.vep.vcf.gz")
)

gvr <- read.gvr(vcf_dir, verbose = FALSE)
dim(gvr)

## ----show-gvr-----------------------------------------------------------------
head(gvr[, .(Hugo_Symbol, Variant_Classification, IMPACT,
    Tumor_Sample_Barcode)])

## ----filter-default-----------------------------------------------------------
filt <- gvr_filter(gvr, verbose = FALSE)
dim(filt)
table(filt$Variant_Classification)

## ----filter-panel-------------------------------------------------------------
gvr_list_panels()                    # registered panels
panel <- gvr_panel_genes("breast cancer")
panel

filt_bc <- gvr_filter(gvr, genes = panel, verbose = FALSE)
filt_bc[, .(Hugo_Symbol, HGVSp_Short, Variant_Classification, IMPACT)]

## ----novel--------------------------------------------------------------------
novel <- gvr_novel(gvr, verbose = FALSE)
dim(novel)
novel[, .(Hugo_Symbol, HGVSp_Short, Variant_Classification, IMPACT)]

## ----summary------------------------------------------------------------------
summ <- gvr_summary(
    gvr,
    save_excel = FALSE, save_pdf = FALSE, save_html = FALSE,
    verbose    = FALSE
)
names(summ)
summ$top_genes

## ----lollipop, echo=FALSE, out.width="80%", fig.cap="Lollipop for TP53 from gvr_lollipop() — InterPro domains shaded"----
knitr::include_graphics("figures/lollipop_TP53.png")

## ----genepos, echo=FALSE, out.width="100%", fig.cap="Gene-structure track for BRCA1 from gvr_genepos.plot() — exons/UTRs from Ensembl REST"----
knitr::include_graphics("figures/genepos_BRCA1.png")

## ----gvr-plot-demo, echo=FALSE, out.width="100%", fig.cap="Top-genes variant matrix from gvr_plot() — illustrative multi-sample data (8 samples, 20 cancer genes). In a real cohort, replace with your own gvr table."----
knitr::include_graphics("figures/gvr_plot_demo.png")

## ----gvr-dashboard-demo, echo=FALSE, out.width="100%", fig.cap="Screenshot of gvr_summary(save_html = TRUE) output — synthetic 8-sample cohort rendered for layout demonstration only, not real data."----
knitr::include_graphics("figures/gvr_summary_dashboard.png")

## ----session-info, echo=FALSE-------------------------------------------------
sessionInfo()

