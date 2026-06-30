test_that("gvr_hpo_genes resolves HPO IDs from a local table", {
    hpo_file <- tempfile(fileext = ".txt")

    writeLines(c(
        "HPO-ID\tHPO-NAME\tGene-ID\tGene-Name",
        "HP:0003002\tBreast carcinoma\t672\tBRCA1",
        "HP:0003002\tBreast carcinoma\t675\tBRCA2",
        "HP:0002664\tNeoplasm\t7157\tTP53"
    ), hpo_file)

    genes <- gvr_hpo_genes(
        hpo = "HP:0003002",
        hpo_path = hpo_file,
        verbose = FALSE
    )

    expect_equal(genes, c("BRCA1", "BRCA2"))
})


test_that("gvr_hpo_genes validates HPO ID format", {
    hpo_file <- tempfile(fileext = ".txt")

    writeLines(c(
        "HPO-ID\tHPO-NAME\tGene-ID\tGene-Name",
        "HP:0003002\tBreast carcinoma\t672\tBRCA1"
    ), hpo_file)

    expect_error(
        gvr_hpo_genes("HP0003002", hpo_path = hpo_file, verbose = FALSE),
        "invalid HPO identifier"
    )
})


test_that("gvr_hpo_genes returns empty vector for unknown valid HPO ID", {
    hpo_file <- tempfile(fileext = ".txt")

    writeLines(c(
        "HPO-ID\tHPO-NAME\tGene-ID\tGene-Name",
        "HP:0003002\tBreast carcinoma\t672\tBRCA1"
    ), hpo_file)

    expect_warning(
        genes <- gvr_hpo_genes(
            hpo = "HP:9999999",
            hpo_path = hpo_file,
            verbose = FALSE
        ),
        "no genes found"
    )

    expect_equal(genes, character(0))
})
