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
#' @param ogroup a character string

plot_og_genes <- function(ogroup)
{
    ### There is an dicrepancy between ids in ogroups and expression matrix
    genes <- og[[ogroup]]
    print(genes)
    gg_stages <- as.factor(substr(colnames(gg_dat), 1, 4))
    plot_gene <- function(gene, dset)
    {
        to_plot <- grep(gene, rownames(dset))
        if(length(to_plot) > 0) {
        stripchart(as.numeric(dset[to_plot, ]) ~ gg_stages,
                   vertical = TRUE, method = "jitter",
                   pch = 16, col = "blue",
                   ylab = "TMPs", main = gene)
        grid()
        }
    }
    in_gg <- grep("^Cg", x = genes)
    in_th <- genes[-in_gg]
    in_th <- sapply(in_th, function(i) strsplit(i, "\\.")[[1]][1])
    in_gg <- genes[in_gg]
    # print(in_gg); print(in_th)
    # in_gg <- in_gg[in_gg %in% rownames(gg_dat)]
    # in_th <- in_th[in_th %in% rownames(th_dat)]
    print(in_gg); print(in_th)
    sapply(in_gg, plot_gene, dset = gg_dat)
    sapply(in_th, plot_gene, dset = th_dat)
    # plot_gene(genes[1])
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
