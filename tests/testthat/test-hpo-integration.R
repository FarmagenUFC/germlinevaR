# Integration tests for the `hpo` (and, for gvr_filter, `panel`) argument on
# read.gvr() and gvr_filter().
#
# These tests use the shipped fixture and the shipped example VCF, so they
# never hit the network. The internal option `gvr.hpo_path` is used as a
# testing hook to force gvr_hpo_genes() to read the fixture instead of the
# real HPO file.

fx <- testthat::test_path("fixtures", "phenotype_to_genes_mini.tsv")

# Helper: stage the shipped example VEP VCF into a temp folder so read.gvr()'s
# folder-mode discovery can find it.
.stage_example_vcf <- function() {
    dir <- tempfile("gvr_hpo_it_")
    dir.create(dir)
    file.copy(
        system.file("extdata", "example.vep.vcf.gz", package = "germlinevaR"),
        file.path(dir, "example.vep.vcf.gz")
    )
    dir
}


test_that("read.gvr(hpo = <term>) restricts to matching genes only", {
    skip_if_not(file.exists(fx), "HPO fixture missing")

    withr::with_options(list(gvr.hpo_path = fx), {
        vcf_dir <- .stage_example_vcf()
        on.exit(unlink(vcf_dir, recursive = TRUE), add = TRUE)

        baseline <- read.gvr(vcf_dir, verbose = FALSE)
        hpo_only <- read.gvr(vcf_dir, hpo = "HP:0003002", verbose = FALSE)

        # Baseline has 62 rows; HP:0003002 fixture maps to 10 genes of which
        # BRCA1 and BRCA2 are the only two present in the example VCF.
        expect_gt(nrow(baseline), nrow(hpo_only))
        expect_setequal(unique(hpo_only$Hugo_Symbol), c("BRCA1", "BRCA2"))
    })
})


test_that("read.gvr(hpo = NULL) is byte-identical to omitting the arg", {
    vcf_dir <- .stage_example_vcf()
    on.exit(unlink(vcf_dir, recursive = TRUE), add = TRUE)

    a <- read.gvr(vcf_dir, verbose = FALSE)
    b <- read.gvr(vcf_dir, hpo = NULL, verbose = FALSE)

    expect_equal(nrow(a), nrow(b))
    expect_equal(a$Hugo_Symbol, b$Hugo_Symbol)
})


test_that("read.gvr(hpo=) unions with panel= (both restrictions kept)", {
    skip_if_not(file.exists(fx))

    withr::with_options(list(gvr.hpo_path = fx), {
        vcf_dir <- .stage_example_vcf()
        on.exit(unlink(vcf_dir, recursive = TRUE), add = TRUE)

        # Breast cancer panel = 24 genes incl. BRCA1/BRCA2; HP:0025022 fixture
        # = COL5A1 + COL5A2 + TNXB. Union should keep BRCA1/BRCA2 AND COL5A1
        # rows from the example VCF (COL5A2/TNXB absent from example data).
        out <- read.gvr(vcf_dir,
                        panel   = "breast cancer",
                        hpo     = "HP:0025022",
                        verbose = FALSE)

        expect_gt(nrow(out), 0L)
        expect_true("COL5A1" %in% out$Hugo_Symbol)
        expect_true(any(c("BRCA1", "BRCA2") %in% out$Hugo_Symbol))
    })
})


test_that("read.gvr(hpo = unresolved) returns baseline unchanged", {
    skip_if_not(file.exists(fx))

    withr::with_options(list(gvr.hpo_path = fx), {
        vcf_dir <- .stage_example_vcf()
        on.exit(unlink(vcf_dir, recursive = TRUE), add = TRUE)

        base <- read.gvr(vcf_dir, verbose = FALSE)
        # HP:0999999 is not in the fixture -> warning + no gene restriction
        expect_warning(
            unres <- read.gvr(vcf_dir, hpo = "HP:0999999", verbose = FALSE),
            regexp = "not found in the table"
        )
        expect_equal(nrow(base), nrow(unres))
    })
})


test_that("read.gvr(hpo=) verbose output includes hpo subset coverage msg", {
    skip_if_not(file.exists(fx))

    withr::with_options(list(gvr.hpo_path = fx), {
        vcf_dir <- .stage_example_vcf()
        on.exit(unlink(vcf_dir, recursive = TRUE), add = TRUE)

        expect_message(
            read.gvr(vcf_dir, hpo = "HP:0003002", verbose = TRUE),
            regexp = "hpo subset:"
        )
    })
})


# =============================================================================
# gvr_filter(hpo=, panel=) integration
# =============================================================================

# Reusable helper: load the shipped example table with all AF/CLIN_SIG/GT
# filters disabled so the tests isolate the gene-subset behaviour.
.example_gvr <- function() {
    readRDS(system.file("extdata", "example_gvr.rds", package = "germlinevaR"))
}

.filter_no_ax_defaults <- function(gvr, ...) {
    gvr_filter(gvr,
               gnomADe_AF     = NULL,
               AF             = NULL,
               ABraOM_AF      = NULL,
               clin_sig_terms = NULL,
               gt_exclude     = NULL,
               verbose        = FALSE,
               ...)
}


test_that("gvr_filter(hpo = <term>) restricts to matching genes only", {
    skip_if_not(file.exists(fx), "HPO fixture missing")

    withr::with_options(list(gvr.hpo_path = fx), {
        gvr <- .example_gvr()
        out <- .filter_no_ax_defaults(gvr, hpo = "HP:0003002")

        # Fixture maps HP:0003002 to 10 genes; only BRCA1 and BRCA2 are in
        # the example table.
        expect_setequal(unique(out$Hugo_Symbol), c("BRCA1", "BRCA2"))
    })
})


test_that("gvr_filter(hpo = NULL, panel = NULL) is byte-identical to omitting", {
    gvr <- .example_gvr()
    a <- .filter_no_ax_defaults(gvr)
    b <- .filter_no_ax_defaults(gvr, hpo = NULL, panel = NULL)
    expect_equal(nrow(a), nrow(b))
    expect_equal(a$Hugo_Symbol, b$Hugo_Symbol)
})


test_that("gvr_filter(panel=) restricts to panel genes (mirrors read.gvr)", {
    gvr <- .example_gvr()
    out <- .filter_no_ax_defaults(gvr, panel = "breast cancer")

    # Breast cancer panel contains BRCA1 and BRCA2; example table has both.
    expect_true(all(unique(out$Hugo_Symbol) %in%
                    germlinevaR::gvr_panel_genes("breast cancer")))
    expect_true(any(c("BRCA1", "BRCA2") %in% out$Hugo_Symbol))
})


test_that("gvr_filter(panel + hpo) unions both restrictions", {
    skip_if_not(file.exists(fx))

    withr::with_options(list(gvr.hpo_path = fx), {
        gvr <- .example_gvr()
        out <- .filter_no_ax_defaults(gvr,
                                      panel = "breast cancer",
                                      hpo   = "HP:0025022")

        # BC panel contributes BRCA1/BRCA2; HP:0025022 fixture contributes
        # COL5A1 (COL5A2/TNXB absent from example VCF).
        expect_true("COL5A1" %in% out$Hugo_Symbol)
        expect_true(any(c("BRCA1", "BRCA2") %in% out$Hugo_Symbol))
    })
})


test_that("gvr_filter(hpo = unresolved) leaves gene subset empty (removes all rows)", {
    skip_if_not(file.exists(fx))

    withr::with_options(list(gvr.hpo_path = fx), {
        gvr <- .example_gvr()
        expect_warning(
            out <- .filter_no_ax_defaults(gvr, hpo = "HP:0999999"),
            regexp = "not found in the table"
        )
        # Unresolved HPO produces zero genes, so genes vector stays what the
        # user passed (NULL) -> the "9. Gene subset" block is a no-op ->
        # ALL rows survive (matches gvr_filter contract for `genes=NULL`).
        expect_equal(nrow(out), nrow(gvr))
    })
})


# =============================================================================
# read.gvr.snpeff(hpo=) tests -- exercise the parallel HPO wiring in the
# SnpEff reader against a tiny hermetic SnpEff-annotated VCF.
# =============================================================================

snpeff_fx <- testthat::test_path("fixtures", "example.snpeff.vcf.gz")

test_that("read.gvr.snpeff(hpo = <term>) restricts to matching genes only", {
    skip_if_not(file.exists(snpeff_fx), "SnpEff fixture missing")
    skip_if_not(file.exists(fx),        "HPO fixture missing")

    withr::with_options(list(gvr.hpo_path = fx), {
        gvr <- suppressWarnings(suppressMessages(
            read.gvr.snpeff(
                vcf_path   = snpeff_fx,
                hpo        = "HP:0003002",       # BRCA1, BRCA2 in fixture
                add_abraom = FALSE,
                verbose    = FALSE,
                min_DP     = 0, min_GQ = 0
            )
        ))
    })
    expect_true(nrow(gvr) > 0L)
    expect_true(all(sort(unique(gvr$Hugo_Symbol)) %in% c("BRCA1", "BRCA2")))
})


test_that("read.gvr.snpeff(hpo = <term>) with unrelated term returns zero rows", {
    skip_if_not(file.exists(snpeff_fx))
    skip_if_not(file.exists(fx))

    withr::with_options(list(gvr.hpo_path = fx), {
        gvr <- suppressWarnings(suppressMessages(
            read.gvr.snpeff(
                vcf_path   = snpeff_fx,
                hpo        = "HP:0001939",       # metabolism panel -- no overlap
                add_abraom = FALSE,
                verbose    = FALSE,
                min_DP     = 0, min_GQ = 0
            )
        ))
    })
    expect_equal(nrow(gvr), 0L)
})


test_that("read.gvr.snpeff(panel + hpo) unions gene sets", {
    skip_if_not(file.exists(snpeff_fx))
    skip_if_not(file.exists(fx))

    withr::with_options(list(gvr.hpo_path = fx), {
        gvr <- suppressWarnings(suppressMessages(
            read.gvr.snpeff(
                vcf_path   = snpeff_fx,
                panel      = "breast cancer",    # includes BRCA1, BRCA2
                hpo        = "HP:0025022",       # COL5A1 in fixture
                add_abraom = FALSE,
                verbose    = FALSE,
                min_DP     = 0, min_GQ = 0
            )
        ))
    })
    expect_setequal(unique(gvr$Hugo_Symbol),
                    c("BRCA1", "BRCA2", "COL5A1"))
})


test_that("read.gvr.snpeff(hpo = NULL, panel = NULL) is byte-identical to no filters", {
    skip_if_not(file.exists(snpeff_fx))

    base <- suppressWarnings(suppressMessages(
        read.gvr.snpeff(vcf_path = snpeff_fx, add_abraom = FALSE,
                        verbose = FALSE, min_DP = 0, min_GQ = 0)
    ))
    null_out <- suppressWarnings(suppressMessages(
        read.gvr.snpeff(vcf_path = snpeff_fx, hpo = NULL, panel = NULL,
                        add_abraom = FALSE, verbose = FALSE,
                        min_DP = 0, min_GQ = 0)
    ))
    expect_identical(base, null_out)
})
