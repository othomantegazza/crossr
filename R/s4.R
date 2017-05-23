# http://adv-r.had.co.nz/S4.html
# http://stackoverflow.com/search?tab=votes&q=user%3a547331%20%5bs4%5d%20is%3aanswe
# library(methods)

#' @import methods
NULL

#' Check Validity of an \code{ogset} class element
#'
#' The function \code{check_ogset} checks that in the provided object:
#' - the names of the orthogroups (`og`) match the rownames of the orthogroup expression set `og_eset` (if any is given),
#' - the row names of og_annos are contained in the names of the orthogroup list
#' - the row names of `colData` (if any) match the column names of the orthogroup expression set `og_eset`,
#' - the row names of `rowData` (if any) match the row names of the orthogroup expression set `og_eset`,
#' - the variable in `design` (if any) are colnames of `colData`,
#' - the rownames of `spec1_colData` match the colnames `spec1_exp` (if any),
#' - the rownames of `spec2_colData` match the colnames `spec2_exp` (if any),
#' - `exp_cond` is a character string
#' - `exp_cond` is contained in columns of `spec1_colData` and `spec2_colData` and that those columns match.
#' @md
#'
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

    # Rownames of og_annos in names orthogroups
    names_og <- names(object@og)
    rownames_annos <- rownames(object@og_annos)
    if (!all(rownames_annos %in% names_og)) {
        msg <- "Rownames of og_annos must be contained in the names of the orthogroup list"
        errors <- c(errors, msg)
    }

    # column names of og_exp in rownames colData
    names_cols <- colnames(object@og_exp)
    names_coldata <- rownames(object@colData)
    if (!all(names_cols == names_coldata)) {
        msg <- "The column names of og_exp must match the row names of colData"
        errors <- c(errors, msg)
    }

    # rownames of og_exp in rownames rowData
    names_rows <- rownames(object@og_exp)
    names_rowdata <- rownames(object@rowData)
    if (!all(dim(object@rowData) == c(0,0))) {
        if(length(names_rows) != length(names_rowdata)) {
            msg <- "The row names of og_exp must match the row names of rowData"
            errors <- c(errors, msg)
        } else if(!all(names_rows == names_rowdata)) {
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
#' None of the arguments is required to initialize a ogset class element
#'
#' @param og A \code{list} of orthogroups
#' @param og_annos a \code{data.frame} with functional annotations for the orthogroups
#' @param og_exp A orthogroup-wise expression \code{data.frame}
#' @param rowData A \code{data.frame} containing the row metadata for the orthogroup expression set
#' @param colData A \code{data.frame} containing the sample info for `og_exp`
#' @param desing A \code{formula} describing the experimental design
#' @param metadata A \code{list} of metadata
#' @param stats A \code{data.frame} with output statistics for `og_exp`
#' @param spec1_exp A \code{data.frame} with expression data for species 1
#' @param spec2_exp A \code{data.frame} with expression data for species 2
#' @param spec1_colData A \code{data.frame} with sample info for `spec1_exp`
#' @param spec2_colData A \code{data.frame} with sample info for `spec2_exp`
#' @param exp_cond A \code{character} string indicating which variable of `design`
#'  represents the experimental condition under which the two species are compared
#' @param og_nomatch orthogroup represented only in one species, filled by \code{collapse_orthologs}
#'
#'
#' @export make_ogset

make_ogset <- setClass("ogset",
                       representation(og = "list",
                                      og_annos = "data.frame",
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

