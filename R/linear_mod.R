#' Calculate F-values from a Linear Model on an Orthogroup Expression Set
#'
#' @param dat \code{data.frame} contains data for any gene
#' @param coldata \code{data.frame} categorical variables with info for every samples
#' @param design \code{formula}, encodes the assumption on how the variables in coldata
#' explain the observed gene expression

make_fit <- function(dat, coldata, design)
{
    vars <- all.vars(design)
    stopifnot(all(vars %in% colnames(coldata)))
    design <- stats::as.formula(paste("exp", deparse(design)))
    fit_anova_int <- function(ogroup, dset, coldata)
    {
        dat_temp <- t(dset)
        to_fit <- cbind(exp = dat_temp[, ogroup], coldata)
        fit <- stats::lm(formula = design, data = to_fit)
        anova_of_fit <- stats::anova(fit)
        out <- anova_of_fit$`F value`
        names(out) <- attributes(anova_of_fit)$row.names
        return(out)
    }

    fits <- parallel::mclapply(rownames(dat),
                     fit_anova_int, dset = dat, coldata = coldata,
                     mc.cores = 4)


    names(fits) <- rownames(dat)
    fits <- as.data.frame(do.call(rbind, fits))
    fits$Residuals <- NULL

    return(fits)

}

#' Estimates and Add F-values from a Linear Model to an og_set class element
#'
#' @param ogset a \code{ogset} class element
#' @param log_scale logical, should the data be natural log trasnsformed before fitting the model? (actually \code{log(data + 1)} transformed)
#'
#' The info on the \code{log_scaled} parameters get stored in the metadata slot
#'
#' @export

add_fit <- function(ogset, log_scale = FALSE)
{
    stopifnot(isS4(ogset))

    if (log_scale) {
        to_fit <- log(ogset@og_exp + 1)
        ogset@metadata <- c(ogset@metadata, "fit performed on log scaled data")
    } else {
        to_fit <- ogset@og_exp
        ogset@metadata <- c(ogset@metadata, "fit performed on not scaled data")
    }

    ogset@stats <- make_fit(dat = to_fit,
                            coldata = ogset@colData,
                            design = ogset@design)

    return(ogset)
}

#' Extract Top Tags from Ogset Class Element
#'
#' @param ogset an ogset class element
#' @param stat \code{character}, the name of the column of ogset@stat used as ranking feature
#' @param n \code{numeric} how many tags should the function extract
#'
#' \code{get_top_tags} extract the top tags from an ogset class element and returns
#' a named numeric vector with the top statistics and the ID of the associated elements
#' (orthogroups ID)
#'
#' if \code{n} is \code{Inf}, the function return all the tags ranked by the ranking statistics
#'
#' @export

get_top_tags <- function(ogset, rank_stat, n = 100)
{
    if(nrow(ogset@stats) < 1) {stop("The stats slot is empty, please fill the stats slot
                                   by running add_fit on the ogset class element")}
    if(n < 1) stop("n must be a positive integer")
    if(n > nrow(ogset@stats) && is.finite(n)) stop(paste0("n is too high, the ranking statistic is available only for ",
                                           nrow(ogset@stats),
                                           " genes"))
    rankd <- ogset@stats[order(ogset@stats[[rank_stat]], decreasing = TRUE), ]
    out <- rankd[[rank_stat]]
    out <- stats::setNames(out, rownames(rankd))
    if(is.infinite(n)) {
        return(out)
        } else {return(out[1:n])}
}
