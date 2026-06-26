## Smoke tests for germlinevaR
## These tests use only the example data shipped with the package
## (inst/extdata/example_gvr.rds and inst/extdata/example.vep.vcf.gz)
## and require no network access or external dependencies.

# ── helpers ──────────────────────────────────────────────────────────────────

example_rds <- function() {
  system.file("extdata", "example_gvr.rds", package = "germlinevaR")
}

example_vcf_dir <- function() {
  vcf_src <- system.file("extdata", "example.vep.vcf.gz",
                         package = "germlinevaR")
  d <- tempfile("gvr_test_")
  dir.create(d)
  file.copy(vcf_src, file.path(d, "example.vep.vcf.gz"))
  d
}

# ── example_gvr.rds ──────────────────────────────────────────────────────────

test_that("example_gvr.rds loads and has expected shape", {
  path <- example_rds()
  expect_true(nzchar(path), label = "example_gvr.rds not found in extdata")

  gvr <- readRDS(path)
  expect_s3_class(gvr, "data.frame")
  expect_equal(nrow(gvr), 62L)
  expect_gte(ncol(gvr), 20L)
})

test_that("example_gvr.rds has required MAF-style columns", {
  gvr <- readRDS(example_rds())
  required_cols <- c(
    "Hugo_Symbol", "Variant_Classification", "Variant_Type",
    "Chromosome", "Start_Position", "End_Position",
    "Reference_Allele", "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode", "HGVSc", "HGVSp"
  )
  missing <- setdiff(required_cols, names(gvr))
  expect_equal(missing, character(0))
})

test_that("example_gvr.rds sample barcode is Sample_01", {
  gvr <- readRDS(example_rds())
  expect_equal(unique(gvr$Tumor_Sample_Barcode), "Sample_01")
})

test_that("example_gvr.rds contains expected Variant_Classification levels", {
  gvr <- readRDS(example_rds())
  expected_vc <- c(
    "Missense_Mutation", "Silent", "Nonsense_Mutation",
    "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site"
  )
  present <- intersect(expected_vc, unique(gvr$Variant_Classification))
  expect_gte(length(present), 4L)
})

# ── read.gvr() ───────────────────────────────────────────────────────────────

test_that("read.gvr() parses the example VCF and returns a data.frame", {
  skip_if_not_installed("data.table")
  d <- example_vcf_dir()
  on.exit(unlink(d, recursive = TRUE))

  gvr <- read.gvr(d, verbose = FALSE)
  expect_s3_class(gvr, "data.frame")
  expect_gte(nrow(gvr), 1L)
  expect_true("Hugo_Symbol" %in% names(gvr))
  expect_true("Variant_Classification" %in% names(gvr))
})

test_that("read.gvr() result has Tumor_Sample_Barcode column", {
  skip_if_not_installed("data.table")
  d <- example_vcf_dir()
  on.exit(unlink(d, recursive = TRUE))

  gvr <- read.gvr(d, verbose = FALSE)
  expect_true("Tumor_Sample_Barcode" %in% names(gvr))
  expect_true(all(nzchar(gvr$Tumor_Sample_Barcode)))
})

test_that("read.gvr() result matches example_gvr.rds dimensions", {
  skip_if_not_installed("data.table")
  d <- example_vcf_dir()
  on.exit(unlink(d, recursive = TRUE))

  gvr_live <- read.gvr(d, verbose = FALSE)
  gvr_ref  <- readRDS(example_rds())
  expect_equal(nrow(gvr_live), nrow(gvr_ref))
})

# ── gvr_filter() ─────────────────────────────────────────────────────────────

test_that("gvr_filter() subsets rows correctly", {
  gvr <- readRDS(example_rds())
  # Suppress expected ABraOM_AF warning (column absent in example data)
  filtered <- suppressWarnings(gvr_filter(gvr, vc = "Missense_Mutation"))
  expect_s3_class(filtered, "data.frame")
  expect_true(nrow(filtered) < nrow(gvr))
  expect_true(all(filtered$Variant_Classification == "Missense_Mutation"))
})

test_that("gvr_filter() with no matching vc returns zero rows", {
  gvr <- readRDS(example_rds())
  filtered <- suppressWarnings(gvr_filter(gvr, vc = "THIS_DOES_NOT_EXIST"))
  expect_equal(nrow(filtered), 0L)
})

# ── gvr_novel() ──────────────────────────────────────────────────────────────

test_that("gvr_novel() returns a data.frame subset of input", {
  gvr <- readRDS(example_rds())
  novel <- gvr_novel(gvr)
  expect_s3_class(novel, "data.frame")
  expect_lte(nrow(novel), nrow(gvr))
  expect_true(all(names(gvr) %in% names(novel)))
})

# ── gvr_summary() ────────────────────────────────────────────────────────────

test_that("gvr_summary() returns a named list with expected elements", {
  gvr <- readRDS(example_rds())
  s <- gvr_summary(gvr,
                   save_excel = FALSE,
                   save_pdf   = FALSE,
                   save_html  = FALSE,
                   verbose    = FALSE)
  expect_type(s, "list")
  expect_gte(length(s), 1L)
})

# ── gvr_list_panels() ────────────────────────────────────────────────────────

test_that("gvr_list_panels() returns a character vector", {
  panels <- gvr_list_panels()
  expect_type(panels, "character")
  expect_gte(length(panels), 1L)
})
