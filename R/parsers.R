#' Extracts TPM for a Salmon Output File
#'
#' \code{get_TPM} is a wrapper for \code{\link{read.csv}} that loads and extracts
#' TPM from a salmon output file provided in \code{path}
#'
#' @param path a character string containing the path to the SALMON output file

get_TPM <- function(path) {
    TPM <- read.csv(file = paste(path, "quant.sf", sep = "/"),
                    sep = "\t", row.names = "Name")[, "TPM", drop = FALSE]
    return(TPM)
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
    ### convert from CgyXXXXX to CgXXXXX  ##This is error prone
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
