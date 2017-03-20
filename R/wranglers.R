#' Sum Expression of Orthologs and Make one Expression Matrix for two species
#'
#' \code{collapse_orthologs_3} takes two expression matrices and a list of orthogroups.
#' It sums the expression of the genes in the orthogroups for the two species separately
#' and outputs a single dataset with orthogroup-wise expression
#'
#' The function assumes that the info on the species are contained in the first two
#' characters of the gene names, and that the gene names are the rownames of the datasets.
#' The samples for each dataset should be the same The numbers of samples per dataset is the
#' input n_cols. Some step are considerably slow and are parallelized with the function
#' \code{mclapply} from the package \code{parallel} the number of threads defaults to 1
#'
#' @param eset1 a \code{data.frame} with the expression matrix for species 1
#' @param eset2 a \code{data.frame} with the expression matrix for species 2
#' @param ogroups a \code{list} of orthogroups, the gene names in the list must match the \code{row.names} of the datasets
#' @param n_cols numeric, the number of sample for each species
#' @param mc.cores numeric, the number of threads

collapse_orthologs_3 <- function(eset1, eset2, ogroups, n_cols = 18, mc.cores = 1)
{
    ## COLLAPSE ORTHOLOGS MIGHT BE FAULTY BECAUSE IT DOES NOT
    ## DEAL WITH DISCREPANCIES BETWEEN TRANSCRIPT NAMES IN
    ## TH ANNOS SHEETH AND IN TX DATA SHEETG

    # require(parallel)
    # require(stringr)

    ## colnames must be the same for collapsing
    # stopifnot(all(colnames(eset1) %in% colnames(eset2)))

    sort_cols <- function(df) df[ , sort(colnames(df))]
    eset1 <- sort_cols(eset1)
    eset2 <- sort_cols(eset2)

    es <- as.data.frame(rbind(eset1, eset2))

    ## This produce orthologs expression df
    reduce_XM <- function(name)
    {
        ifelse(grepl("XM", name),
               strsplit(name, "\\.")[[1]][1],
               name)
    }
    ogroups <- lapply(ogroups, function(i) sapply(i, reduce_XM))

    ogroups_es <- parallel::mclapply(ogroups, function(i) es[ grep( paste( i, collapse = "|" ), rownames(es)), ],
                                     mc.cores = mc.cores)

    ## Some ID in the ogroup is unmatched in the eset, why?
    unmatched <- ogroups_es[!sapply(ogroups_es, function(i) all(complete.cases(i)))]
    ogroups_es <- lapply(ogroups_es, function(i) i[complete.cases(i), ])

    ## Assumption that homologs in single species
    ## have same function
    sum_hom <- function(df)
    {
        spec <- substr(rownames(df), 1, 2)
        sum_on_rows <- function(i)
        {
            tapply(i, spec, sum)
        }
        lapply(df, sum_on_rows)
    }
    ogroups_es <- parallel::mclapply(ogroups_es, sum_hom,
                                     mc.cores = mc.cores)
    ogroups_es <- lapply(ogroups_es, unlist)

    ## remove orthogroups that have expression data for only one species
    ogroups_es <- ogroups_es[vapply(ogroups_es, length, numeric(1)) == n_cols*2]

    ## collapse to df
    ogroups_es <- do.call(rbind, ogroups_es)

    ## Make colnames clearer spec_stage_rep
    colnames(ogroups_es) <- paste(stringr::str_sub(colnames(ogroups_es), -2, -1),
                                  stringr::str_sub(colnames(ogroups_es), -6, -4), sep = "_")
    ogroups_es <- ogroups_es[, order(colnames(ogroups_es))]
    colnames(ogroups_es) <- sub("Cg", "gg", colnames(ogroups_es))
    colnames(ogroups_es) <- sub("XM", "th", colnames(ogroups_es))

    return(list(og_es = ogroups_es, unm = unmatched))
    #return(list(eset1, eset2, es))
}


#' Converts Protein ID into Transcript ID in Orthogroup List
#'
#' Orthogroups are generated on protein sequence, \code{px_2_tx} converts protein ID to transcript ID in an orthogroup list
#'
#' @param groups a \code{list} of orthogroups
#' @param anno_file a character string containing the path to the annotation file in tabulated format
#' @param tx_id a \code{character} string indicating the name of the column of transcript ids in the annotation file
#' @param px_id a \code{character} string indicating the name of the column of protein ids in the annotation file
#' @param mc.cores \code{numeric}, the number of threads used for \code{mclapply}

px_2_tx <- function(ogroups,
                    anno_file,
                    tx_id,
                    px_id,
                    mc.cores = 1)
{
    # th_id_file <- "~/Google Drive/Cross_species_comparison/genomes/t_hassleriana_from_ncbi/GCF_000463585.1_ASM46358v1_feature_table.txt"
    ids <- read.table(annofile,
                         sep = "\t", header = TRUE,
                         stringsAsFactors = FALSE,
                         quote = "")
    # head(th_ids)
    # tst <- th_dat[rownames(th_dat2) %in% th_ids$product_accession, ]; nrow(tst)
    # tst <- th_dat[!rownames(th_dat2) %in% th_ids$product_accession, ]; nrow(tst)
    switch_ids <- function(id) {
        if(id %in% as.character(ids[, px_id])) {
            return(ids[which(as.character(ids[, px_id])), tx_id])
        } else return(id)
    }
    ogroups <- parallel::mclapply(ogroups, function(i) {
        unname(sapply(i, switch_ids))},
        mc.cores = mc.cores)
    return(ogroups)
}
