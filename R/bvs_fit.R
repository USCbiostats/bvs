bvs_fit <- function(z,
                    num_active,
                    y,
                    x,
                    n,
                    p,
                    rare,
                    hap,
                    region_ind,
                    num_regions,
                    forced,
                    p_forced,
                    family_func,
                    control,
                    weights,
                    offset,
                    mustart,
                    m)
{

    # If rare variants make risk index and new design matrix
    if (rare) {
        if (num_active > 0) {
            region_z <- region_ind * z
            nonnull_col <- colSums(region_z) > 0
            data.i <- x %*% region_z[, nonnull_col]
            num_vars <- sum(nonnull_col)
        } else {
            data.i <- NULL
            nonnull_col <- rep(FALSE, num_regions)
            num_vars <- 0
        }
    } else if (hap && num_active > 1) {
        data.i <- hapBVS(data[, colnames(x)[z]], min.Hap.freq = 0.02)
    } else {
        data.i <- x[, z, drop = FALSE]
        num_vars <- num_active
    }

    # Fit glm
    fit <- suppressWarnings(glm_fit_custom(x = cbind(rep(1, n), data.i, forced),
                                           y = y,
                                           nobs = n,
                                           nvars = 1 + num_vars + p_forced,
                                           weights = weights,
                                           mustart = mustart,
                                           m = m,
                                           offset = offset,
                                           family = family_func,
                                           control = control))

    fit$num_vars <- num_vars
    if (rare && !all(nonnull_col)) {
        coef <- rep(0.0, num_regions)
        coef[nonnull_col] <- fit$coef
        fit$coef <- coef
    }
    return(fit)
}





