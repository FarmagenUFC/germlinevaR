test_that("gvr_filter can filter by HPO-derived genes", {
    hpo_file <- tempfile(fileext = ".txt")

    writeLines(c(
        "HPO-ID\tHPO-NAME\tGene-ID\tGene-Name",
        "HP:0003002\tBreast carcinoma\t672\tBRCA1",
        "HP:0003002\tBreast carcinoma\t675\tBRCA2"
    ), hpo_file)

    gvr <- data.table::data.table(
        Hugo_Symbol = c("BRCA1", "BRCA2", "TP53"),
        gnomADe_AF = c("", "", ""),
        AF = c("", "", ""),
        ABraOM_AF = c("", "", ""),
        CLIN_SIG = c("pathogenic", "pathogenic", "pathogenic"),
        BIOTYPE = c("protein_coding", "protein_coding", "protein_coding"),
        GT = c("0/1", "0/1", "0/1"),
        Variant_Classification = c(
            "Missense_Mutation",
            "Missense_Mutation",
            "Missense_Mutation"
        )
    )

    out <- gvr_filter(
        gvr,
        hpo = "HP:0003002",
        hpo_path = hpo_file,
        verbose = FALSE
    )

    expect_equal(sort(unique(out$Hugo_Symbol)), c("BRCA1", "BRCA2"))
    expect_false("TP53" %in% out$Hugo_Symbol)
})


test_that("gvr_filter combines manual genes and HPO-derived genes by union", {
    hpo_file <- tempfile(fileext = ".txt")

    writeLines(c(
        "HPO-ID\tHPO-NAME\tGene-ID\tGene-Name",
        "HP:0003002\tBreast carcinoma\t672\tBRCA1"
    ), hpo_file)

    gvr <- data.table::data.table(
        Hugo_Symbol = c("BRCA1", "TP53", "MYC"),
        gnomADe_AF = c("", "", ""),
        AF = c("", "", ""),
        ABraOM_AF = c("", "", ""),
        CLIN_SIG = c("pathogenic", "pathogenic", "pathogenic"),
        BIOTYPE = c("protein_coding", "protein_coding", "protein_coding"),
        GT = c("0/1", "0/1", "0/1"),
        Variant_Classification = c(
            "Missense_Mutation",
            "Missense_Mutation",
            "Missense_Mutation"
        )
    )

    out <- gvr_filter(
        gvr,
        genes = "TP53",
        hpo = "HP:0003002",
        hpo_path = hpo_file,
        verbose = FALSE
    )

    expect_equal(sort(unique(out$Hugo_Symbol)), c("BRCA1", "TP53"))
    expect_false("MYC" %in% out$Hugo_Symbol)
})
