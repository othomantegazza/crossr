context("px_2_tx returns a list of orthogroups with transcript ids")

test_that("px_2_tx returns list of orthogroups", {
    # load(file = "crossr/tests/testthat/sample_100_og_pxid.Rdata")
    og <- crossr::parse_orthogroups("sample_100_OrthologousGroup.csv")

    expect_is(crossr::px_2_tx(ogroups = og,
                                    anno_file = "~/Google Drive/Cross_species_comparison/genomes/t_hassleriana_from_ncbi/GCF_000463585.1_ASM46358v1_feature_table.txt",
                                    tx_id = "product_accession",
                                    px_id = "related_accession",
                                    mc.cores = 4),
                    "list")
})
