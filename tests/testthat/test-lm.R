# require(testthat)

context("linear model and ANOVA work correctly")

test_that("make_fit takes any design within coldata", {

    dat <- read.table("sample_oges.txt")

    coldata <- data.frame(spc = substr(colnames(dat), 1, 2),
                          stg = substr(colnames(dat), 4, 4),
                          row.names = colnames(dat))

    expect_error(make_fit(dat = dat, coldata = coldata, design = ciao, mc.cores = 1))
    expect_s3_class(make_fit(dat = dat,
                            coldata = coldata,
                            design = ~ 0 + spc + stg + spc:stg,
                            mc.cores = 1), "data.frame")
})

