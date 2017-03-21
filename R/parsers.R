#' Extracts TPM from a Salmon Output File
#'
#' \code{get_TPM} is a wrapper for \code{\link{read.csv}} that loads and extracts
#' TPM from a salmon output file (quant.sf) provided in \code{path}
#'
#' @param path a character string, the path to the (quant.sf) Salmon output

get_TPM <- function(path) {
    stopifnot(grepl("quant.sf$", path))
    TPM <- read.csv(file = path,
                    sep = "\t", row.names = "Name")[, "TPM", drop = FALSE]
    return(TPM)
}


#' Make Expression Matrix from All Salmon Output Folders
#'
#' \code{Make_TPM_df} is a wrapper for \code{\link{get_TPM}}, it takes as input the path
#' for the folder that contain the Salmon output folders for all your sample and returns
#' the TPM expression matrix as a \code{data.frame} object
#'
#' @param path character string, the path for the folder that contains the Salmon results for all the samples

make_TPM_df <- function(path) {
    stopifnot(is.character(path) & length(path) == 1)
    paths <- list.files(path = path,
                        full.names = TRUE,
                        recursive = TRUE)
    paths <- grep("quant.sf$", paths, value = TRUE)
    if(length(paths) == 0) stop("No quant.sf file in the folder provided or in its subfolder")
    dat <- lapply(paths, get_TPM)

    ## from here http://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-vector
    zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
        if (length(x) == 1) return(TRUE)
        x <- range(x) / mean(x)
        isTRUE(all.equal(x[1], x[2], tolerance = tol))
    }
    stopifnot(zero_range(vapply(dat, nrow, numeric(1))))

    names(dat) <- sub(paste0(normalizePath(path),"/"), "", paths)
    names(dat) <- sub("/quant.sf", "", names(dat))

    dat <- lapply(names(dat), function(i) {
        colnames(dat[[i]]) <- i
        return(dat[[i]])
    })
    dat <- data.frame(dat, check.rows = TRUE)

    return(dat)
}



#' Parse and Load an Orthogroup File Provided by Orthofinder
#'
#' \code{parse_orthogroups} parses an orthogroup file provided by Orthofinder and
#' outputs the result as a list of orthogroups
#'
#' @param path_2_ogroups a character string that provides the path to the orthogroup
#' file

parse_orthogroups <- function(path_2_ogroups)
{
    orthogroups <- scan(path_2_ogroups, what = "",
                        skip = 1, sep = "\n")

    ### Don't know why some values are separated by tab in file
    ### This fixes it
    orthogroups <- gsub("\t", ", ",
                        orthogroups)

    ### why are GG identifiers different in mRNA and proteins
    ### convert from CgyXXXXX to CgXXXXX  ##This might be prone
    orthogroups <- gsub("y", "",
                        orthogroups)
    orthogroups <- strsplit(orthogroups, ", ")

    ### names of orthogroup are read in first list element
    names(orthogroups) <- sapply(orthogroups,
                                 function(i) i[1])
    orthogroups <- lapply(orthogroups,
                          function(i) i[-1])

    ## Dont know why some elements in orthogroups are null ( "" )
    orthogroups <- lapply(orthogroups,
                          function(i) i[i != ""])


    return(orthogroups)
}
