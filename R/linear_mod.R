#' Add F-value from Linear Model to Orthologroup Expression Set

add_fit <- function(dat, coldata)
{
    fit_anova_int <- function(ogroup, dset, coldata)
    {
        dat_temp <- t(dset)
        to_fit <- cbind(exp = dat_temp[, ogroup], coldata)
        fit <- lm(exp ~ 0 + spc + stg + spc:stg, data = to_fit)
        anova_of_fit <- anova(fit)
        out <- anova_of_fit$`F value`
        names(out) <- attributes(anova_of_fit)$row.names
        return(out)
    }

    fits <- mclapply(rownames(dat),
                     fit_anova_int, dset = dat, coldata = coldata,
                     mc.cores = 4)

    names(fits) <- rownames(dat)
    fits <- as.data.frame(do.call(rbind, fits))

    dat_fit <- merge(dat, fits, by = "row.names"); rownames(dat_fit) <- dat_fit$Row.names; dat_fit$Row.names <- NULL
    dat_fit <- dat_fit[order(dat_fit$`spc:stg`, decreasing = T), ]
    return(dat_fit)
}
