# require(devtools)
# devtools::load_all("~/Google Drive/Cross_species_comparison/crossr")

context("orthogroup parser")

test_that("output of parse_orthogroups is list of character vectors", {
    ogroups <- parse_orthogroups("sample_100_OrthologousGroup.csv")

    expect_is(ogroups, "list")
    expect_identical(unique(sapply(ogroups, class)), "character")
})

test_that("no tab signs in orthogroup output", {
    ogroups <- parse_orthogroups("sample_100_OrthologousGroup.csv")

    expect_false(any(grepl("\t", ogroups)))
})

test_that("orthogroup name is correct", {
    ogroups <- parse_orthogroups("sample_100_OrthologousGroup.csv")

    expect_named(ogroups)
    expect_true(all(grepl("^OG", names(ogroups))))
})

test_that("no empty element in orthogroups", {
    ogroups <- parse_orthogroups("sample_100_OrthologousGroup.csv")
    empty_slots <- sum(unlist(lapply(ogroups, function(i) sum(i == ""))))

    expect_equal(empty_slots, 0)
})
