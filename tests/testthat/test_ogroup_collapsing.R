context("collapse orthogroup return dataframe of collapsed orthogroups")

test_that("collapse_orthogroups returns dataframe",{
    th_dat <- read.csv2(file = "sample_expr_th.csv",
              row.names = 1)


    expect_s3_class(collapse_orthologs(th_dat, ogroups = ogroups), "data.frame")
})
