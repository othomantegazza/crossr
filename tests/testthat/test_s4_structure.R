context("s4 class ogset can take empty elements and tests for unwanted faults")

test_that("can make empty og_set", {
    expect_s4_class(make_ogset(), "ogset")
})

test_that("check ogset tests for coldata and rowdata matching", {
    load(file = "tos4.Rdata")

    rowdata <- data.frame(seq_along(og_eset$thal_0h_1))
    rownames(rowdata) <- rownames(og_eset)

    expect_s4_class(make_ogset(og = ogroups,
                               og_exp = og_eset,
                               colData = coldata,
                               rowData = rowdata,
                               design = ~ spc + treat + spc:treat,
                               stats = og_fit,
                               spec1_exp = thaliana,
                               spec2_exp = lyrata),
                    "ogset")

    og_eset_WRONG <- rbind(og_eset,
                           not_an_orthogroup = seq_along(og_eset[1,]))
    expect_error(make_ogset(og = ogroups,
                            og_exp = og_eset_WRONG))

    coldata_WRONG <- coldata[sample(rownames(coldata)), ]
    expect_error(make_ogset(og = ogroups,
                            og_exp = og_eset,
                            colData = coldata_WRONG))

    rowdata_WRONG <- rowdata[sample(rownames(rowdata)), , drop = FALSE]
    expect_error(make_ogset(og = ogroups,
                            og_exp = og_eset,
                            rowData = rowdata_WRONG))
})
