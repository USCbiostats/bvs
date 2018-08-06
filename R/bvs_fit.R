bvs_fit <- function(z,
                    num_active,
                    y,
                    x,
                    n,
                    p,
                    intercept,
                    rare,
                    hap,
                    region_ind,
                    num_regions,
                    forced,
                    p_forced,
                    family_func,
                    prior_coef,
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
            num_vars <- sum(nonnull_col)
            xgamma <- cbind(x %*% region_z[, nonnull_col], forced)
        } else {
            xgamma <- forced
            nonnull_col <- rep(FALSE, num_regions)
            num_vars <- 0
        }
    } else if (hap && num_active > 1) {
        xgamma <- cbind(hap_bvs(x[, colnames(x)[z]], min_hap_freq = 0.02), forced)
    } else {
        xgamma <- cbind(x[, z, drop = FALSE], forced)
        num_vars <- num_active
    }

    # Fit glm
    fit <- suppressWarnings(glm_fit_custom(x = xgamma,
                                           y = y,
                                           nobs = n,
                                           nvars = intercept + p_forced + num_vars,
                                           weights = weights,
                                           mustart = mustart,
                                           m = m,
                                           offset = offset,
                                           family = family_func,
                                           control = control))

    fit$num_vars <- num_vars

    # compute coef
    if (num_vars > 0) {
        if (rare) {
            coef <- rep(0.0, num_regions)
            coef[nonnull_col] <- fit$coef[1:num_vars]
            fit$coef <- coef
        } else {
            fit$coef <- fit$coef[1:num_vars]
        }

        if (prior_coef[[1]] == "gprior") {
            fit$coef <- fit$coef * prior_coef[[2]] / (1 + prior_coef[[2]])
        }
    }

    # compute marginal log-likelihood
    if (prior_coef[[1]] == "gprior") {
        if (p_forced + num_vars > 0) {
            if (intercept){
                s_temp <- crossprod(y - mean(y)) - (prior_coef[[2]] / (1 + prior_coef[[2]])) * crossprod(y - mean(y), xgamma) %*% solve(crossprod(xgamma)) %*% crossprod(xgamma, y - mean(y))
            } else {
                s_temp <- crossprod(y) - (prior_coef[[2]] / (1 + prior_coef[[2]])) * crossprod(y, xgamma) %*% solve(crossprod(xgamma)) %*% crossprod(xgamma, y)
            }
            fit$marg_ll <- -(0.5 * num_vars) * log(prior_coef[[2]] + 1) - 0.5 * (2 * prior_coef[[3]] + n - 1) * log(2 * prior_coef[[4]] + s_temp)
        } else {
            if (intercept) {
                fit$marg_ll <- - 0.5 * (2 * prior_coef[[3]] + n - 1) * log(2 * prior_coef[[4]] + crossprod(y - mean(y)))
            } else {
                fit$marg_ll <- - 0.5 * (2 * prior_coef[[3]] + n - 1) * log(2 * prior_coef[[4]] + crossprod(y))
            }
        }
    } else {
        fit$marg_ll <- -(0.5 * fit$deviance + num_vars)
    }
    return(fit)
}





