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
    design <- as.formula(paste("exp", deparse(design)))
    fit_anova_int <- function(ogroup, dset, coldata)
    {
        dat_temp <- t(dset)
        to_fit <- cbind(exp = dat_temp[, ogroup], coldata)
        fit <- lm(formula = design, data = to_fit)
        anova_of_fit <- anova(fit)
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
