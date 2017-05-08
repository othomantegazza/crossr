context("parse quant.sf files")

test_that("output of get_TPM is a 1 numeric column data.frame", {
    expect_equal(class(get_TPM("sample_quant.sf")), "data.frame")
    expect_equal(ncol(get_TPM("sample_quant.sf")), 1)
    expect_equal(class(get_TPM("sample_quant.sf")$TPM), "numeric")
})

