#' Add F-value from Linear Model to Orthologroup Expression Set

add_fit <- function(dat, coldata, formula)
{
    vars <- all.vars(formula)
    vars <- vars[!vars == "exp"]
    stopifnot(all(vars %in% colnames(coldata)))
    fit_anova_int <- function(ogroup, dset, coldata)
    {
        dat_temp <- t(dset)
        to_fit <- cbind(exp = dat_temp[, ogroup], coldata)
        fit <- lm(formula = formula, data = to_fit)
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
model.matrix(coldata, 0 + spc + stg + spc:stg)
add_fit(dat = dat, coldata = coldata, formula = exp ~ 0 + spc + stg + spc:stg)
}
