# Offline unit tests for the HPO helper (gvr_hpo_genes).
#
# All tests use the shipped fixture:
#   tests/testthat/fixtures/phenotype_to_genes_mini.tsv
# so they never hit the network. Integration with read.gvr() and gvr_filter()
# is tested in test-hpo-integration.R.

# Locate the fixture next to this test file so it works both from
# devtools::test() (cwd = tests/testthat) and R CMD check.
fx <- testthat::test_path("fixtures", "phenotype_to_genes_mini.tsv")

test_that("gvr_hpo_genes: canonical HP:NNNNNNN resolves to expected symbols", {
    skip_if_not(file.exists(fx), "HPO fixture missing")

    g <- gvr_hpo_genes("HP:0003002", hpo_path = fx, verbose = FALSE)

    # Fixture defines 10 genes for HP:0003002 (Breast carcinoma):
    #   ATM BARD1 BRCA1 BRCA2 BRIP1 CHEK2 NF1 PALB2 PMS2 TP53
    expect_type(g, "character")
    expect_length(g, 10L)
    expect_true(all(c("BRCA1", "BRCA2", "TP53", "PALB2") %in% g))

    # Output must be uppercased, unique, and sorted.
    expect_identical(g, sort(unique(toupper(g))))
})


test_that("gvr_hpo_genes: multiple terms union without duplication", {
    skip_if_not(file.exists(fx))

    # HP:0003002 has BRCA1 + BRCA2; HP:0025022 has COL5A1 + COL5A2 + TNXB.
    # Union should be 13 distinct genes (10 + 3, no overlap in fixture).
    u <- gvr_hpo_genes(
        c("HP:0003002", "HP:0025022"),
        hpo_path = fx,
        verbose  = FALSE
    )
    expect_length(u, 13L)
    expect_true(all(c("BRCA1", "COL5A1") %in% u))
    # No duplicates
    expect_identical(u, unique(u))
})


test_that("gvr_hpo_genes: lenient input forms all normalise to canonical", {
    skip_if_not(file.exists(fx))

    # Canonical form ("HP:0003002") is the reference.
    ref <- gvr_hpo_genes("HP:0003002", hpo_path = fx, verbose = FALSE)

    # Four lenient forms must produce the same output.
    for (input in list("hp:0003002", "hp:3002", "3002", "0003002")) {
        out <- suppressMessages(
            gvr_hpo_genes(input, hpo_path = fx, verbose = FALSE)
        )
        expect_identical(out, ref,
                         info = paste("Input:", input))
    }
})


test_that("gvr_hpo_genes: unresolved term warns and returns empty vector", {
    skip_if_not(file.exists(fx))

    # HP:0999999 is not present in the fixture.
    expect_warning(
        w <- gvr_hpo_genes("HP:0999999", hpo_path = fx, verbose = FALSE),
        regexp = "not found in the table"
    )
    expect_identical(w, character(0))
})


test_that("gvr_hpo_genes: garbled input errors with actionable message", {
    skip_if_not(file.exists(fx))

    expect_error(
        gvr_hpo_genes("SPAM", hpo_path = fx, verbose = FALSE),
        regexp = "unrecognisable HPO identifier"
    )
})


test_that("gvr_hpo_genes: NULL/empty input returns character(0) silently", {
    skip_if_not(file.exists(fx))

    expect_identical(gvr_hpo_genes(NULL),             character(0))
    expect_identical(gvr_hpo_genes(character(0)),     character(0))
    expect_identical(gvr_hpo_genes(""),               character(0))
    expect_identical(gvr_hpo_genes(c("", NA_character_)), character(0))
})


test_that("gvr_hpo_genes: verbose emits per-code resolution message", {
    skip_if_not(file.exists(fx))

    expect_message(
        gvr_hpo_genes(c("HP:0003002", "HP:0025022"),
                      hpo_path = fx,
                      verbose  = TRUE),
        regexp = "per-code:.*HP:0003002=10.*HP:0025022=3"
    )
})


test_that("gvr_hpo_genes: missing local file gives clear error", {
    expect_error(
        gvr_hpo_genes("HP:0003002", hpo_path = "/tmp/no_such_file_here.tsv"),
        regexp = "local HPO file not found"
    )
})


test_that("gvr_hpo_genes: fread skip logic ignores leading # comment lines", {
    skip_if_not(file.exists(fx))

    # Fixture already ships with a leading `#format:` comment line, so the
    # canonical resolve already exercises the skip path. Also verify that a
    # fixture WITHOUT any comment (header on line 1) still parses.
    tmp <- tempfile(fileext = ".tsv")
    on.exit(unlink(tmp), add = TRUE)
    writeLines(c(
        "hpo_id\thpo_name\tncbi_gene_id\tgene_symbol\tdisease_id",
        "HP:0003002\tBreast carcinoma\t672\tBRCA1\tOMIM:604370",
        "HP:0003002\tBreast carcinoma\t675\tBRCA2\tOMIM:612555"
    ), tmp)
    out <- gvr_hpo_genes("HP:0003002", hpo_path = tmp, verbose = FALSE)
    expect_identical(out, c("BRCA1", "BRCA2"))
})
