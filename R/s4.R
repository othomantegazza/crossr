# http://adv-r.had.co.nz/S4.html
# http://stackoverflow.com/search?tab=votes&q=user%3a547331%20%5bs4%5d%20is%3aanswe
# library(methods)

check_ogset <- function(object) {
    errors <- character()

    names_og <- names(object@og)
    rownames_exp <- rownames(object@og_exp)
    if (!all(rownames_exp %in% names_og)) {
        msg <- "Rownames of og_exp must be contained in the names of the orthogroup list"
        errors <- c(errors, msg)
    }

    names_cols <- colnames(object@og_exp)
    names_coldata <- rownames(object@colData)
    if (!all(names_cols == names_coldata)) {
        msg <- "The column names of og_exp must match the row names of colData"
        errors <- c(errors, msg)
    }

    names_rows <- rownames(object@og_exp)
    names_rowndata<- rownames(object@rowData)
    if (!all(dim(object@rowData) == c(0,0))) {
        if(!all(names_rows == names_rowndata)) {
            msg <- "The row names of og_exp must match the row names of rowData"
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
