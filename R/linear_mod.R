#' Add F-values From Linear Model to Orthologroup Expression Set
#'
#' @param dat \code{data.frame} contains data for any gene
#' @param coldata \code{data.frame} categorical variables with info for every samples
#' @param design \code{formula}, encodes the assumption on how the variables in coldata
#' explain the observed gene expression

add_fit <- function(dat, coldata, design)
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

    dat_fit <- merge(dat, fits, by = "row.names"); rownames(dat_fit) <- dat_fit$Row.names; dat_fit$Row.names <- NULL
    dat_fit <- dat_fit[order(dat_fit$`spc:stg`, decreasing = T), ]
    return(dat_fit)
}

if(FALSE) {
    load("tests/testthat/sample_oges.Rdata")
    require(devtools)
    load_all()
    coldata <- data.frame(spc = substr(colnames(dat), 1, 2),
                          stg = substr(colnames(dat), 4, 4),
                          row.names = colnames(dat))
    tst1 <- add_fit(dat = dat, coldata = coldata, design = exp ~ 0 + spc + stg + spc:stg)
    tst2 <- add_fit(dat = dat, coldata = coldata, design = . ~ 0 + spc + stg + spc:stg)
    tst3 <- add_fit(dat = dat, coldata = coldata, design = ~ 0 + spc + stg + spc:stg)


    dt_tst <- cbind(exp = dat[1, ], coldata)
    fit1 <- lm(exp ~ 0 + spc + stg + spc:stg, data = dt_tst)
    fit2 <- lm( ~ 0 + spc + stg + spc:stg, data = dt_tst)
}
