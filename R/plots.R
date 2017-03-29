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
#' @param ogroup a character string
#' @param dset a numeric data.frame with ogroups as rownames
#' @param ylab a character string
#' @param clusters a named numeric vector


plot_all_stages <- function(ogroup,
                            dset,
                            ylab = "normalized expression",
                            clusters = "cls_here")
{
    # old <- par()
    # on.exit(par(old), add = TRUE)
    par(mar = c(5, 4, 8, 2))

    coldata <- substr(colnames(dset[, 1:36]), 1, 4)
    th_points <- split(dset[ogroup, grep("^th", colnames(dset)), drop = TRUE],
                       as.factor(grep("^th", coldata, value = TRUE)))
    gg_points <- split(dset[ogroup, grep("^gg", colnames(dset)), drop = TRUE],
                       as.factor(grep("^gg", coldata, value = TRUE)))
    stripchart(th_points,
               frame = FALSE, vertical = TRUE,
               ylim = range(unlist(c(th_points, gg_points))),
               pch = 16, method = "jitter",
               ylab = ylab,
               xlab = "Stage",
               xaxt = "n",
               main = paste(ogroup, #"\n",
                            # paste(ogroup_es_both_annos[[ogroup]]$annos, collapse = "\n"), "\n",
                            # paste(og[[ogroup]], collapse = " "),
                            # gsub(" - ", "\n", dset[ogroup, "annos_th"]),
                            dset[ogroup, "annos_th"],
                            paste0("F-value = ", round( dset[ogroup, "spc:stg"], 2 )  ),
                            paste0("genes in ogroup: ", length(og[[ogroup]]),
                                   "; of which ", length( grep("Cg", og[[ogroup]]) ), " from Gg, ",
                                   "and ", length( grep("XM", og[[ogroup]]) ), " from Th"),
                            cls <- ifelse( any(clusters == "cls_here"),
                                           "",
                                           ifelse(ogroup %in% names( clusters ),
                                                  paste0("in cluster ", unname(clusters[ogroup]) ),
                                                  "not in clusters") ),
                            sep = " \n "))
    stages <- length(unique(coldata))/2
    axis(side = 1, at = 1:stages, labels = paste("Stage", 0:(stages - 1), sep = ""))
    stripchart(gg_points,
               frame = FALSE, vertical = TRUE, add = TRUE,
               pch = 17, method = "jitter",
               ylab = "size adjusted CPM",
               xlab = "species - stage")
    lines(1:stages, sapply(th_points, function(i) median(unlist(i))), lty = 2)
    lines(1:stages, sapply(gg_points, function(i) median(unlist(i))), lty = 3)
    legend("topright", legend = c("C3", "C4"), lty = c(2, 3), bty = "n")
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
