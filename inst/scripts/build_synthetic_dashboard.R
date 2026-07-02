## ============================================================================
##  build_synthetic_dashboard.R
##
##  Reproducible generator for the synthetic 8-sample / 19-gene / 106-variant
##  cohort used to render the HTML dashboard screenshot shown in
##  README.md and vignette("germlinevaR").
##
##  IMPORTANT — this is NOT real patient data. Every field (positions,
##  variant classes, IMPACT, CLIN_SIG, samples) is fabricated purely to
##  populate `gvr_summary()`'s KPI cards and bar charts so the layout can
##  be demonstrated at a glance. Do NOT read biological meaning into any
##  specific value.
##
##  Layout locked by design:
##    - 19 germline / hereditary-cancer genes (fixed order below)
##    - Per-gene total variant burden: 1..15 (fixed vector `gene_burden`)
##    - 8 samples (Sample_01 .. Sample_08)
##    - Per-gene x sample cell counts drawn by Iterative Proportional
##      Fitting (IPF) toward target sample totals, then rounded via
##      multinomial draws seeded by `set.seed(20260702)`
##    - Variant classes drawn from a fixed 12-level distribution
##    - IMPACT determined by the class (VEP severity)
##
##  Usage:
##      source(system.file("scripts", "build_synthetic_dashboard.R",
##                         package = "germlinevaR"))
##      # `gvr_synth` is now a data.table in the calling environment
##      gvr_summary(gvr_synth, save_html = TRUE, out_dir = ".")
##
##  Or from source:
##      Rscript inst/scripts/build_synthetic_dashboard.R
## ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(20260702)

## ---- Design parameters ------------------------------------------------------

## Gene order (matches gvr_plot_demo.png right-hand burden bar, low -> high)
genes <- c("CDH1", "APC", "MLH1", "PMS2", "VHL", "RB1", "PTEN", "SMAD4",
           "MSH2", "BRAF", "PALB2", "MSH6", "KRAS", "STK11", "CHEK2", "ATM",
           "BRCA2", "TP53", "BRCA1")
stopifnot(length(genes) == 19)

## Per-gene total burden (matches demo right-hand bar); sum = 106
gene_burden <- c(
  CDH1  = 1L, APC   = 1L, MLH1  = 2L, PMS2  = 2L, VHL   = 2L, RB1   = 2L,
  PTEN  = 3L, SMAD4 = 3L, MSH2  = 4L, BRAF  = 4L, PALB2 = 5L, MSH6  = 5L,
  KRAS  = 5L, STK11 = 7L, CHEK2 = 9L, ATM   = 10L,
  BRCA2 = 12L, TP53  = 14L, BRCA1 = 15L
)
stopifnot(sum(gene_burden) == 106L)

## 8 samples with target per-sample totals (~11..17 from a visual read of the
## demo). The IPF step below rescales these so only the RATIO between samples
## matters — the absolute sum can differ from sum(gene_burden) without harm.
samples <- sprintf("Sample_%02d", 1:8)
sample_targets <- c(Sample_01 = 12, Sample_02 = 11, Sample_03 = 15,
                    Sample_04 = 13, Sample_05 = 15, Sample_06 = 13,
                    Sample_07 = 16, Sample_08 = 17)

## Real GRCh38 gene loci (start coordinate for offset base)
gene_loc <- list(
  CDH1  = list(chr = 16, start = 68737225),  APC   = list(chr = 5,  start = 112707500),
  MLH1  = list(chr = 3,  start = 36993350),  PMS2  = list(chr = 7,  start = 5970925),
  VHL   = list(chr = 3,  start = 10141635),  RB1   = list(chr = 13, start = 48303256),
  PTEN  = list(chr = 10, start = 87863625),  SMAD4 = list(chr = 18, start = 51028394),
  MSH2  = list(chr = 2,  start = 47403067),  BRAF  = list(chr = 7,  start = 140719327),
  PALB2 = list(chr = 16, start = 23603160),  MSH6  = list(chr = 2,  start = 47783145),
  KRAS  = list(chr = 12, start = 25205246),  STK11 = list(chr = 19, start = 1205798),
  CHEK2 = list(chr = 22, start = 28687820),  ATM   = list(chr = 11, start = 108222484),
  BRCA2 = list(chr = 13, start = 32315474),  TP53  = list(chr = 17, start = 7668402),
  BRCA1 = list(chr = 17, start = 43044295)
)
stopifnot(setequal(names(gene_loc), genes))

## Variant_Classification distribution (12 canonical levels)
vc_levels <- c("Missense_Mutation", "Silent", "Nonsense_Mutation",
               "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site",
               "In_Frame_Del", "In_Frame_Ins", "Translation_Start_Site",
               "Intron", "5'UTR", "3'UTR")
vc_probs  <- c(0.38, 0.15, 0.10, 0.08, 0.06, 0.07, 0.04, 0.03, 0.03,
               0.03, 0.02, 0.01)
stopifnot(length(vc_levels) == length(vc_probs), abs(sum(vc_probs) - 1) < 1e-6)

## VEP IMPACT map (severity that gvr_summary bins into KPI + bar chart)
impact_map <- c(Missense_Mutation      = "MODERATE",
                Silent                 = "LOW",
                Nonsense_Mutation      = "HIGH",
                Frame_Shift_Del        = "HIGH",
                Frame_Shift_Ins        = "HIGH",
                Splice_Site            = "HIGH",
                In_Frame_Del           = "MODERATE",
                In_Frame_Ins           = "MODERATE",
                Translation_Start_Site = "HIGH",
                Intron                 = "MODIFIER",
                `5'UTR`                = "MODIFIER",
                `3'UTR`                = "MODIFIER")

## MAF Variant_Type map (SNP / INS / DEL)
vt_map <- c(Missense_Mutation      = "SNP",
            Silent                 = "SNP",
            Nonsense_Mutation      = "SNP",
            Frame_Shift_Del        = "DEL",
            Frame_Shift_Ins        = "INS",
            Splice_Site            = "SNP",
            In_Frame_Del           = "DEL",
            In_Frame_Ins           = "INS",
            Translation_Start_Site = "SNP",
            Intron                 = "SNP",
            `5'UTR`                = "SNP",
            `3'UTR`                = "SNP")

## ---- Iterative Proportional Fitting (IPF) for per-gene x sample cell counts -

## Start with the outer product of row and column marginals, then scale rows
## and columns alternately until both marginal sums match. Finally round each
## row via a multinomial draw so integer cell counts sum to `gene_burden` per
## gene (row totals are exact; column totals drift due to multinomial variance).
cell <- outer(gene_burden, sample_targets) /
        sum(sample_targets)                    # rows sum to gene_burden

for (iter in seq_len(30)) {
  ## Scale columns to hit sample_targets
  col_sums <- colSums(cell)
  cell     <- sweep(cell, 2, sample_targets / pmax(col_sums, 1e-9), `*`)
  ## Scale rows back to hit gene_burden
  row_sums <- rowSums(cell)
  cell     <- sweep(cell, 1, gene_burden  / pmax(row_sums, 1e-9), `*`)
}

## Multinomial rounding per gene (row totals become exactly gene_burden)
cell_int <- matrix(0L, nrow = nrow(cell), ncol = ncol(cell),
                   dimnames = dimnames(cell))
for (g in rownames(cell)) {
  p_g <- cell[g, ] / sum(cell[g, ])
  cell_int[g, ] <- as.integer(
    rmultinom(1, size = as.integer(gene_burden[g]), prob = p_g))
}
stopifnot(all(rowSums(cell_int) == gene_burden))   # exact per-gene totals

## ---- Build 106-row variant table --------------------------------------------

n_variants   <- sum(cell_int)
rows_out     <- vector("list", n_variants)
row_idx      <- 0L
allele_pairs <- list(c("A", "G"), c("C", "T"), c("G", "A"), c("T", "C"),
                     c("A", "C"), c("G", "T"))

for (g in genes) {
  loc <- gene_loc[[g]]
  for (s in samples) {
    n_gs <- cell_int[g, s]
    if (n_gs == 0L) next
    for (k in seq_len(n_gs)) {
      row_idx <- row_idx + 1L
      pos     <- loc$start + sample(-30000:30000, 1)
      vc      <- sample(vc_levels, 1L, prob = vc_probs)
      vt      <- vt_map[[vc]]
      ap      <- allele_pairs[[sample.int(length(allele_pairs), 1L)]]
      ref     <- ap[1]
      alt     <- if (vt == "INS")      paste0(ref, sample(c("A","C","G","T"), 1))
                 else if (vt == "DEL") "-"
                 else                   ap[2]
      if (vt == "DEL") ref <- paste0(ref, sample(c("A","C","G","T"), 1))
      rows_out[[row_idx]] <- data.table(
        Hugo_Symbol            = g,
        Chromosome             = as.character(loc$chr),
        Start_Position         = pos,
        End_Position           = if (vt == "DEL") pos + 1L else pos,
        Reference_Allele       = ref,
        Tumor_Seq_Allele2      = alt,
        Variant_Classification = vc,
        Variant_Type           = vt,
        IMPACT                 = impact_map[[vc]],
        HGVSp_Short            = paste0("p.", sample(c("R","K","D","E","Q","P","L"), 1),
                                        sample(10:800, 1),
                                        sample(c("*","fs","del","H","W","C","G"), 1)),
        CLIN_SIG               = NA_character_,
        Tumor_Sample_Barcode   = s
      )
    }
  }
}

gvr_synth <- rbindlist(rows_out, use.names = TRUE, fill = FALSE)

## CLIN_SIG populated on ~50% of HIGH-impact rows (mimics ClinVar coverage)
high_idx <- which(gvr_synth$IMPACT == "HIGH")
n_clin   <- as.integer(round(length(high_idx) * 0.5))
if (n_clin > 0L) {
  clin_pick <- sample(high_idx, n_clin)
  gvr_synth[clin_pick,
            CLIN_SIG := sample(
              c("pathogenic", "likely_pathogenic", "uncertain_significance",
                "benign", "likely_benign"),
              n_clin, replace = TRUE,
              prob = c(0.35, 0.30, 0.20, 0.10, 0.05))]
}

## ---- Sanity ----------------------------------------------------------------

stopifnot(
  nrow(gvr_synth) == 106L,
  length(unique(gvr_synth$Hugo_Symbol))          == 19L,
  length(unique(gvr_synth$Tumor_Sample_Barcode)) == 8L,
  identical(sort(unique(gvr_synth$Hugo_Symbol)), sort(genes))
)

## Per-gene burden matches design exactly
burden_check <- gvr_synth[, .N, by = Hugo_Symbol][
                          , setNames(N, Hugo_Symbol)][names(gene_burden)]
stopifnot(all(burden_check == gene_burden))

message(sprintf(
  "[build_synthetic_dashboard] built %d variants, %d genes, %d samples",
  nrow(gvr_synth),
  length(unique(gvr_synth$Hugo_Symbol)),
  length(unique(gvr_synth$Tumor_Sample_Barcode))))

## Available in the calling environment as `gvr_synth`.
invisible(gvr_synth)
