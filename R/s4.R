# http://adv-r.had.co.nz/S4.html
# http://stackoverflow.com/search?tab=votes&q=user%3a547331%20%5bs4%5d%20is%3aanswe
# library(methods)

#' Check Validity of an \code{ogset} class element
#'
#' @param object an \code{ogset} class element
#'
#' @export

check_ogset <- function(object) {
    errors <- character()

    # Rownames of og_exp in names orthogroups
    names_og <- names(object@og)
    rownames_exp <- rownames(object@og_exp)
    if (!all(rownames_exp %in% names_og)) {
        msg <- "Rownames of og_exp must be contained in the names of the orthogroup list"
        errors <- c(errors, msg)
    }

    # column names of og_exp in rownames colData
    names_cols <- colnames(object@og_exp)
    names_coldata <- rownames(object@colData)
    if (!all(names_cols == names_coldata)) {
        msg <- "The column names of og_exp must match the row names of colData"
        errors <- c(errors, msg)
    }

    # rownames og og_exp in rownames rowData
    names_rows <- rownames(object@og_exp)
    names_rowndata <- rownames(object@rowData)
    if (!all(dim(object@rowData) == c(0,0))) {
        if(!all(names_rows == names_rowndata)) {
            msg <- "The row names of og_exp must match the row names of rowData"
            errors <- c(errors, msg)
        }}

    # variables in design must be in columns names of colData
    cols_coldata <- colnames(object@colData)
    vars <- all.vars(object@design)
    if (length(vars) > 0 && !is.null(cols_coldata)) {
        if (!all(vars %in% cols_coldata)) {
            msg <- "The variables in design must be contained in the colnames of colData"
            errors <- c(errors, msg)
        }}

    # rownames spec1_colData are equals to colnames spec1_exp
    names_coldata <- rownames(object@spec1_colData)
    names_cols <- colnames(object@spec1_exp)
    if (length(names_coldata) > 0 && length(names_cols) > 0) {
        if (!length(names_coldata) == length(names_cols)) {
            msg <- "The number of rows of spec1_colData must be the same the n of columns of spec1_exp"
            errors <- c(errors, msg)
        } else if (!all(names_cols == names_coldata)) {
            msg <- "Rownames of spec1_colData must be the same as the colnames of spec1_exp"
            errors <- c(errors, msg)
        }}

    # rownames spec2_colData are equals to colnames spec2_exp
    names_coldata <- rownames(object@spec2_colData)
    names_cols <- colnames(object@spec2_exp)
    if (length(names_coldata) > 0 && length(names_cols) > 0) {
        if (!length(names_coldata) == length(names_cols)) {
            msg <- "The number of rows of spec2_colData must be the same the n of columns of spec2_exp"
            errors <- c(errors, msg)
        } else if (!all(names_cols == names_coldata)) {
            msg <- "Rownames of spec2_colData must be the same as the colnames of spec2_exp"
            errors <- c(errors, msg)
        }}

    # experimental condition is character string
    exp_cond <- object@exp_cond
    if(length(exp_cond) > 0) {
        if(length(exp_cond) > 1) {
            msg <- "exp_cond must be a character string"
            errors <- c(errors, msg)
        }
    }

    # exp_cond matches in species datasets
    exp_cond <- object@exp_cond
    spc1_cd <- object@spec1_colData
    spc2_cd <- object@spec2_colData
    if(length(exp_cond) == 1) {
        if(!{exp_cond %in% colnames(spc1_cd) &&
                exp_cond %in% colnames(spc2_cd)}) {
            msg <- "exp_cond must be contained in columns of spec1_colData and spec2_colData"
            errors <- c(errors, msg)
        } else if(!{{all(unique(spc1_cd[[exp_cond]]) %in%
                       unique(spc2_cd[[exp_cond]]))} &&
                       {all(unique(spc2_cd[[exp_cond]]) %in%
                            unique(spc1_cd[[exp_cond]]))}}) {
            msg <- "experimental conditions in spec1_colData and spec2_colData must be the same"
            errors <- c(errors, msg)
        }}

    if (length(errors) == 0) TRUE else errors
}


#' S4 Class, Container for Ortholog Expression Data
#'
#' @export make_ogset

make_ogset <- setClass("ogset",
                       representation(og = "list",
                                      og_exp = "data.frame",
                                      rowData = "data.frame",
                                      colData = "data.frame",
                                      design = "formula",
                                      metadata = "list",
                                      stats = "data.frame",
                                      spec1_exp = "data.frame",
                                      spec2_exp = "data.frame",
                                      spec1_colData = "data.frame",
                                      spec2_colData = "data.frame",
                                      exp_cond = "character",
                                      og_nomatch = "data.frame"),
                       validity = check_ogset)

# og_S4_set <- new("ogset",
#                  og = ogroups,
#                  og_exp = og_eset,
#                  colData = coldata,
#                  design = ~ spc + treat + spc:treat)


# # tst og_set ----------------------------------------------------------------
# og_eset_WRONG <- rbind(og_eset, spiripucci = seq_along(og_eset[1,]))
# tail(og_eset_WRONG)
#
# coldata_WRONG <- coldata[sample(rownames(coldata)), ]
# rowdata_WRONG <- data.frame(a = "A")
# rownames(rowdata_WRONG)
#
# tst <- new("ogset",
#            og = ogroups,
#            og_exp = og_eset_WRONG,
#            colData = coldata_WRONG,
#            rowData = rowdata_WRONG,
#            design = ~ spc + treat + spc:treat)
#
# og_S4_set[1:2, ]


"a" %in%
    "a"
