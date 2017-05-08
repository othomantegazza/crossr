context("switch_ids returns a list of orthogroups with changed ids")

test_that("switch_ids returns list of orthogroups", {

    path2og <- system.file("exdata", "ogroups_sample.csv", package = "crossr")
    og <- parse_orthogroups(path2og)

    path2annos_at <- system.file("exdata", "thaliana_features_sample.txt", package = "crossr")
    thaliana_ids <- read.table(file = path2annos_at,
                               sep = "\t", header = FALSE,
                               stringsAsFactors = FALSE,
                               quote = "")

    ## Not sure why the header is commented with "#"
    colnames(thaliana_ids) <-  scan(file = path2annos_at,
                                    what = "",
                                    nlines = 1)[-1]

    expect_is(switch_ids(ogroups = og,
                          ids_table = thaliana_ids,
                          px_id = "product_accession",
                          tx_id = "related_accession",
                          mc.cores = 1),
                    "list")
})
