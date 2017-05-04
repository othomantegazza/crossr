context("parse quant.sf files")

test_that("output of get_TPM is a 1 numeric column data.frame", {
    expect_equal(class(get_TPM("../../inst/exdata/quant.sf")), "data.frame")
    expect_equal(ncol(get_TPM("../../inst/exdata/quant.sf")), 1)
    expect_equal(class(get_TPM("../../inst/exdata/quant.sf")$TPM), "numeric")
})

