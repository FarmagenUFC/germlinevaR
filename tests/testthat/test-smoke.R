# Smoke tests for germlinevaR
#
# Scope: end-to-end sanity checks against the shipped 62-variant example VCF
# (inst/extdata/example.vep.vcf.gz, single sample Sample_01). These tests do
# NOT exercise the full filter / plot / summary machinery — they verify that
# the three core entry points run and return the documented shape on the
# fixture that ships with the package. The fixture is offline, deterministic,
# and small (~67 KB), so the whole file runs in well under a second on Bioc CI.
#
# Each test stages the fixture into a temp directory because read.gvr() takes
# a *directory* (it scans for `*.vcf.gz`), not a file path.

stage_example_vcf <- function() {
  vcf_dir <- tempfile("gvr_smoke_")
  dir.create(vcf_dir)
  file.copy(
    system.file("extdata", "example.vep.vcf.gz", package = "germlinevaR"),
    file.path(vcf_dir, "example.vep.vcf.gz")
  )
  vcf_dir
}

test_that("read.gvr() returns the documented shape on the example VCF", {
  vcf_dir <- stage_example_vcf()
  gvr <- read.gvr(vcf_dir, verbose = FALSE)

  expect_s3_class(gvr, "data.table")
  # Shape contract (matches the vignette and README hero example)
  expect_equal(nrow(gvr), 62L)
  expect_equal(ncol(gvr), 116L)
  # Sanity check on essential columns
  expect_true(all(c("Hugo_Symbol", "Variant_Classification", "IMPACT",
                    "Tumor_Sample_Barcode", "dbSNP_RS")
                  %in% names(gvr)))
  expect_equal(unique(gvr$Tumor_Sample_Barcode), "Sample_01")
})

test_that("gvr_filter() default thresholds reduce variants to a known count", {
  vcf_dir <- stage_example_vcf()
  gvr  <- read.gvr(vcf_dir, verbose = FALSE)
  filt <- gvr_filter(gvr, verbose = FALSE)

  expect_s3_class(filt, "data.table")
  # Defaults remove rows (rare / clinically relevant / called-genotype pipeline)
  expect_lt(nrow(filt), nrow(gvr))
  # Exact value depends on the documented default thresholds. If defaults are
  # tightened or loosened on purpose, update this expectation in lockstep.
  expect_equal(nrow(filt), 7L)
})

test_that("gvr_novel() returns rows with no rsID", {
  vcf_dir <- stage_example_vcf()
  gvr <- read.gvr(vcf_dir, verbose = FALSE)
  nv  <- gvr_novel(gvr, verbose = FALSE)

  expect_s3_class(nv, "data.table")
  # The invariant of gvr_novel(): every returned row must have no dbSNP rsID
  # (either NA or empty string — read.gvr() uses "" for missing).
  expect_true(all(is.na(nv$dbSNP_RS) | nv$dbSNP_RS == ""))
  # Documented count on this fixture (3 candidates with no rsID / no AF)
  expect_equal(nrow(nv), 3L)
})
