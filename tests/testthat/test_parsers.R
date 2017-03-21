require(devtools)
require(testthat)

devtools::load_all("~/Google Drive/Cross_species_comparison/crossr")

context("parse quant.sf files")

test_that("output of get_TPM is a 1 numeric column data.frame", {
    expect_equal(class(get_TPM("data/quant.sf")), "data.frame")
    expect_equal(ncol(get_TPM("data/quant.sf")), 1)
    expect_equal(class(get_TPM("data/quant.sf")$TPM), "numeric")
})

