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

test_that("check ogset for variable in design and coldata match", {
    load(file = "tos4.Rdata")
    coldata <- data.frame(spc = rep(1, 12),
                          treat = rep(1, 12),
                          stringsAsFactors = FALSE)
    design <- ~ spc + treat + spc*treat
    design_wrong <- ~ spc + that + spc*that

    expect_s4_class(make_ogset(colData = coldata,
                               design = design),
                    "ogset")
    expect_error(make_ogset(coldata = coldata,
                            design = design_wrong))

})

test_that("check consistency of species colData", {
    load(file = "tos4.Rdata")

    coldata_ly <- data.frame(treat = rep(c("a", "b", "c"),
                                         each = 2),
                             stringsAsFactors = FALSE)
    rownames(coldata_ly) <- colnames(lyrata)

    coldata_th <- data.frame(treat = rep(c("a", "b", "c"),
                                         each = 2),
                             stringsAsFactors = FALSE)
    rownames(coldata_th) <- colnames(thaliana)

    expect_s4_class(make_ogset(spec1_exp = lyrata,
                               spec2_exp = thaliana,
                               spec1_colData = coldata_ly,
                               spec2_colData = coldata_th),
                    "ogset")

    ## test: check for spec1 consistency
    coldata_ly_WRONG <- rbind(coldata_ly, bye = "bye")
    expect_error(make_ogset(spec1_exp = lyrata,
                            spec2_exp = thaliana,
                            spec1_colData = coldata_ly_WRONG,
                            spec2_colData = coldata_th))

    coldata_ly_WRONG2 <- rbind(coldata_ly[1:5, , drop = FALSE], bye = "bye")
    expect_error(make_ogset(spec1_exp = lyrata,
                            spec2_exp = thaliana,
                            spec1_colData = coldata_ly_WRONG2,
                            spec2_colData = coldata_th))

    ## test: check for spec2 consistency
    coldata_th_WRONG <- rbind(coldata_th, bye = "bye")
    expect_error(make_ogset(spec1_exp = lyrata,
                            spec2_exp = thaliana,
                            spec1_colData = coldata_ly,
                            spec2_colData = coldata_th_WRONG))

    coldata_th_WRONG2 <- rbind(coldata_th[1:5, , drop = FALSE], bye = "bye")
    expect_error(make_ogset(spec1_exp = lyrata,
                            spec2_exp = thaliana,
                            spec1_colData = coldata_ly,
                            spec2_colData = coldata_th_WRONG2))

    ## test: check for exp_cond length and consistency
    exp_cond <- "treat"
    expect_s4_class(make_ogset(spec1_exp = lyrata,
                               spec2_exp = thaliana,
                               spec1_colData = coldata_ly,
                               spec2_colData = coldata_th,
                               exp_cond = exp_cond), "ogset")

    exp_cond_WRONG <- c("treat", "spc")
    expect_error(make_ogset(exp_cond = exp_cond_WRONG))

    coldata_th_WRONG <- rbind(coldata_th, bye = "bye")
    expect_error(make_ogset(spec1_colData = coldata_ly,
                            spec2_colData = coldata_th_WRONG,
                            exp_cond = exp_cond))
})
