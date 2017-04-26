# require(testthat)

context("linear model and ANOVA work correctly")

test_that("make_fit takes any design within coldata", {
    load(file = "sample_oges.Rdata")
    coldata <- data.frame(spc = substr(colnames(dat), 1, 2),
                          stg = substr(colnames(dat), 4, 4),
                          row.names = colnames(dat))

    expect_error(make_fit(dat = dat, coldata = coldata, design = ciao))
    expect_s3_class(make_fit(dat = dat,
                            coldata = coldata,
                            design = ~ 0 + spc + stg + spc:stg), "data.frame")
})


test_that("add_fit deals with log_scale", {
    FALSE
})
