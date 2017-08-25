#' Plot the Expression of One Orthogroup in Both Species
#'
#' \code{plot_all_stages} is a wrapper for \code{stripchart} and plots a stripchart
#' and a line plot of the expression of an orthogroup across all samples
#'
#' This function takes an orthogroup ID and an ogset class element as arguments.
#' It plots a stripchart of the expression of all the stages; The title of the
#' stripchart contains details on: how many genes are in the orthogroup, the cluster
#' that contain the ogroup (if any), the F-value for the interaction term.
#'
#' @param orthogroup a \code{character} string
#' @param ogset an ogset class element
#' @param species a string used to encode the species info in
#' the design formula
#' @param condition a string used to encode the experimental
#' conditions of the samples in the columns of the dset
#' @param main the title of the plot
#' @param use_annos \code{logical}; if TRUE, functional annotation from og_annos will be displayed in
#'  the plot subtitle. Defaults to FALSE
#' @param annos_col a \code{character} string with the name of the column of og_annos that contains
#'  the functional annotation that should be displayed. Defaults to NULL. It must be provided if \code{use_annos}
#'  is set to TRUE.
#' @param ... arguments to be passed to \code{stripchart}
#'
#' @export


plot_all_stages <- function(orthogroup,
                            ogset,
                            species,
                            condition,
                            main = "",
                            use_annos = FALSE,
                            annos_col = NULL,
                            ...)
{
    old_mar <- graphics::par()$mar
    on.exit(graphics::par(mar = old_mar))
    if(use_annos) stopifnot(annos_col %in% colnames(ogset@og_annos))
    # stopifnot(ncol(dset) == length(species) && ncol(dset) == length(condition))

    # if(use_annos) par(mar = c(3, 4, 7, 2))

    spc_levels <- unique(ogset@colData[[species]])
    spec_fac <- ogset@colData[[species]]
    cond_fac <- ogset@colData[[condition]]
    dset <- ogset@og_exp
    a_points <- split(dset[orthogroup, spec_fac == spc_levels[1], drop = TRUE],
                       as.factor(cond_fac[spec_fac == spc_levels[1]]))
    a_points <- lapply(a_points, unlist)

    b_points <- split(dset[orthogroup, spec_fac == spc_levels[2], drop = TRUE],
                      as.factor(cond_fac[spec_fac == spc_levels[2]]))
    b_points <- lapply(b_points, unlist)

    f_annos <- ifelse(use_annos,
                      yes = gsub(pattern = "--",
                                 replacement = "\n",
                                 x = ogset@og_annos[orthogroup, annos_col, drop = TRUE]),
                      no = "")

    main <- ifelse(main == "", orthogroup, main)
    main <- ifelse(use_annos, paste(main, f_annos, sep = "\n"), main)


    graphics::stripchart(a_points,
                         frame = FALSE, vertical = TRUE,
                         ylim = range(unlist(c(a_points, b_points))),
                         pch = 16, method = "jitter",
                         main = main,
                         ...)
    graphics::stripchart(b_points,
                         frame = FALSE, vertical = TRUE, add = TRUE,
                         pch = 17, method = "jitter")
    graphics::lines(1:length(a_points), vapply(a_points, stats::median, numeric(1)), lty = 2)
    graphics::lines(1:length(b_points), vapply(b_points, stats::median, numeric(1)), lty = 3)
    graphics::legend("topright",
                     legend = c(spc_levels[1], spc_levels[2]),
                     lty = c(2, 3), bty = "n")
}

#' Plot the Expression of One Orthogroup in Both Species with ggplot2
#'
#' \code{ggplot_all_stages} is a wrapper for ggplot2 that plots a stripchart
#' and a line plot of the expression of an orthogroup in both species.
#'
#' This function takes one (ot many) orthogroup IDs and an ogset class element as arguments.
#' It plots a stripchart of the expression of all the stages.
#'
#' You can also provided a vector of multiple orthogroup IDs, \code{ggplot_all_stages} will plot
#' all of them using faceting.
#'
#' The function returns a ggplot object that can be further customized with all basic
#' ggplot2 functionalities.
#'
#' @param orthogroups a \code{character} vector with the orthogroups for the plot
#' @param ogset an ogset class element
#' @param species_var a \code{character} string with the name of column of coldata that encodes
#' for the species from which the samples originates.
#' @param condition_var a \code{character} string with the name of column of coldata that encodes the experimental
#' conditions under which you are comparing the samples
#' @param colours a \code{character} vector with the colours for the plot
#' @param use_annos \code{logical}; if \code{TRUE} the functional annotation from og_annos will
#'  be displayed in the plot subtitle. Defaults to \code{FALSE}
#'
#' @export

ggplot_all_stages <- function(orthogroups,
                               ogset,
                               species_var,
                               condition_var,
                               colours = c("darkblue", "#56B4E9"),
                               use_annos = FALSE)
{
    if(use_annos) {
        if(nrow(ogset@og_annos) == 0) {
            stop("the og_annos slot is empty,
        please provide functional annotations in og_annos slot of ogset or set use_annos = FALSE")
        }
    }
    dset <- ogset@og_exp
    dset <- lapply(orthogroups, function(i) {
        dat <- cbind(t(dset[i, ]), ogset@colData)
        colnames(dat)[1] <- "TPM"
        dat$ogroup <- i
        return(dat)
    })
    dset <- do.call(rbind, dset)
    if(use_annos) {
        for(i in 1:length(dset$ogroup)) {
            dset$ogroup[i] <- paste(dset$ogroup[i],
                                    ogset@og_annos[dset$ogroup[i], 1],
                                    sep = "\n")
        }
    }
    dset$ogroup <- as.factor(dset$ogroup)
    p_out <- ggplot2::ggplot(data = dset,
                             ggplot2::aes_string(x = condition_var,
                                                 y = "TPM",
                                                 colour = species_var,
                                                 group = species_var,
                                                 pch = species_var)) +
        ggplot2::geom_jitter(width = .1,
                             height = 0) +
        ggplot2::stat_summary(fun.y = stats::median,
                              geom= "line",
                              lty = 2) +
        ggplot2::scale_colour_manual(values = colours,
                                     name = "Species") +
        ggplot2::scale_shape(name = "Species") +
        ggplot2::xlab(label = "Condition") +
        ggplot2::facet_wrap(~ ogroup,
                            scales = "free_y")
    return(p_out)

}

#' Plot the Expression of all Genes within one Orthogroup
#'
#' \code{plot_og_gene} is a wrapper for \code{stripchart} and plots a stripchart
#' of the expression of all the gene within a orthogroup for each species in which
#' expression data are available.
#'
#' @param ogroup a character string with the ID of the orthogroups that contains
#' the genes to be plotted
#' @param ogset an ogset class element
#' @param ylab a character string, the name of the y-axis in the plots
#'
#' @export

plot_og_genes <- function(ogroup,
                          ogset,
                          ylab = "TPM")
{
    genes <- ogset@og[[ogroup]]
    print(paste("the ", ogroup, " orthogroup contains ", length(genes), " genes: ",
                paste(genes, collapse = ", "),
                sep = ""))
    plot_gene <- function(gene, dset, cdata, ylab)
    {
        to_plot <- grep(gene, rownames(dset))
        if(length(to_plot) > 0) {
            graphics::stripchart(as.numeric(dset[to_plot, ]) ~ cdata,
                                 vertical = TRUE, method = "jitter",
                                 pch = 16, col = "blue",
                                 main = gene, ylab = ylab)
            graphics::grid()
        }
    }
    genes_spec1 <- genes[genes %in% rownames(ogset@spec1_exp)]
    genes_spec2 <- genes[genes %in% rownames(ogset@spec2_exp)]
    tmp <- sapply(genes_spec1,
                  plot_gene,
                  dset = ogset@spec1_exp,
                  cdata = ogset@spec1_colData[[ogset@exp_cond]],
                  ylab = ylab)
    tmp <- sapply(genes_spec2,
                  plot_gene,
                  dset = ogset@spec2_exp,
                  cdata = ogset@spec2_colData[[ogset@exp_cond]],
                  ylab = ylab)
}



#' Plot an Histogram of Orthogroup Dimension
#'
#' \code{explore_ogroups} is a wrapper for the \code{plot} method for \code{table}
#'
#' @param groups a \code{list} of orthogroups
#' @param main an overall title for the plot
#'
#' @export

explore_ogroups <- function(groups, main = "dimension of orthogroups")
{
    stopifnot(is.list(groups))
    graphics::plot(table(sapply(groups, length)),
         frame = F, ylab = "number of groups",
         xlab = "genes per group", main = main)
    graphics::grid()
}
