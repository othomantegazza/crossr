#' Plot the Expression of One Orthogroup in Both Species
#'
#' \code{plot_all_stages} is a wrapper for \code{stripchart} and plots a stripchart
#' and a line plot of the expression of an orthogroup across all samples
#'
#' This function takes an orthogroup ID and the expression matrix.
#' It plots a stripchart of the expression of all the stages; The title of the
#' stripchart contains details on: how many genes are in the orthogroup, the cluster
#' that contain the ogroup (if any), the F-value for the interaction term.
#'
#' @param orthogroup a \code{character} string
#' @param dset a numeric \code{data.frame} with orthogroups id as row names
#' @param species a \code{character} vector containing the species of the
#'  samples in the columns of the dset
#' @param condition a \code{character} vector containing the experimental
#' conditions of the samples in the columns of the dset
#' @param main the title of the plot


plot_all_stages <- function(orthogroup,
                            dset,
                            species,
                            condition,
                            main = "",
                            ...)
{
    old_mar <- par()$mar
    on.exit(par(mar = old_mar))
    stopifnot(ncol(dset) == length(species) && ncol(dset) == length(condition))

    # par(mar = c(5, 4, 8, 2))

    spc_levels <- unique(species)
    a_points <- split(dset[orthogroup, species == spc_levels[1], drop = TRUE],
                       as.factor(condition[species == spc_levels[1]]))
    a_points <- lapply(a_points, unlist)

    b_points <- split(dset[orthogroup, species == spc_levels[2], drop = TRUE],
                      as.factor(condition[species == spc_levels[2]]))
    b_points <- lapply(b_points, unlist)

    stripchart(a_points,
               frame = FALSE, vertical = TRUE,
               ylim = range(unlist(c(a_points, b_points))),
               pch = 16, method = "jitter",
               main = ifelse(main == "", orthogroup, main),
               ...)
    stripchart(b_points,
               frame = FALSE, vertical = TRUE, add = TRUE,
               pch = 17, method = "jitter")
    lines(1:length(a_points), vapply(a_points, median, numeric(1)), lty = 2)
    lines(1:length(b_points), vapply(b_points, median, numeric(1)), lty = 3)
    legend("topright",
           legend = c(spc_levels[1], spc_levels[2]),
           lty = c(2, 3), bty = "n")
}

#' plot the expression of all genes within one orthogroup
#'
#' \code{plot_og_gene} is a wrapper for \code{stripchart} and plots a stripchart
#' of the expression of all the gene within a orthogroup for each species in which
#' expression data are available
#'
#' @param ogroup a character string with the orthogroup ID for the genes to plot
#' @param og_list a \code{list} containing the mapping from orthogroups to genes
#' @param eset_spec1 the expression set for species 1 as \code{data.frame},
#'  it must contain the gene ID in the rownames
#' @param eset_spec2 the expression set for species 2 as \code{data.frame},
#' it must contain the gene ID in the rownames
#' @param coldata_spec1 a \code{character} vector or a \code{factor} with the
#'  experimental condition for the samples in species 1
#' @param coldata_spec2 a \code{character} vector or a \code{factor} with the
#'  experimental condition for the samples in species 2


plot_og_genes <- function(ogroup,
                          og_list,
                          eset_spec1,
                          eset_spec2,
                          coldata_spec1,
                          coldata_spec2,
                          ylab = "TPM")
{
    ### There is an dicrepancy between ids in ogroups and expression matrix
    genes <- og_list[[ogroup]]
    print(genes)
    plot_gene <- function(gene, dset, cdata, ylab)
    {
        to_plot <- grep(gene, rownames(dset))
        if(length(to_plot) > 0) {
        stripchart(as.numeric(dset[to_plot, ]) ~ cdata,
                   vertical = TRUE, method = "jitter",
                   pch = 16, col = "blue",
                   main = gene, ylab = ylab)
        grid()
        }
    }
    genes_spec1 <- genes[genes %in% rownames(eset_spec1)]
    print(genes_spec1)
    genes_spec2 <- genes[genes %in% rownames(eset_spec2)]
    print(genes_spec2)
    tmp <- sapply(genes_spec1, plot_gene,
                  dset = eset_spec1,
                  cdata = coldata_spec1,
                  ylab = ylab)
    tmp <- sapply(genes_spec2, plot_gene,
                  dset = eset_spec2,
                  cdata = coldata_spec2,
                  ylab = ylab)
}

#' Plot the Expression of Every Orthogroup that Contains a Keyword in their Annos
#'
#' \code{plot_keyword} is a wrapper for \code{\link{plot_all_stages}} that applies this
#' function to all the orthogroups in the dataset that contain a user specified keyword
#' in their functional annotation
#'
#' The function uses grep to search for the \code{keyword} in the functional annotation
#' of the orthogroups
#'
#' @param keyword a character string
#' @inheritParams plot_all_stages

plot_keyword <- function(keyword, dset, clusters)
{

    to_plot <- rownames(dset[grep(keyword, dat_fit_log$annos_th), ])
    sapply(to_plot, plot_all_stages, dset = dset, clusters = clusters)
}

#' Plot an Histogram og Orthogroup Dimension
#'
#' \code{explore_ogroups} is a wrapper for the \code{plot} method for \code{table}
#'
#' @param groups a \code{list} of orthogroups
#' @param main an overall title for the plot

explore_ogroups <- function(groups, main = "dimension of orthogroups")
{
    stopifnot(is.list(groups))
    plot(table(sapply(groups, length)),
         frame = F, ylab = "number of groups",
         xlab = "genes per group", main = main)
    grid()
}
