
bvs_enumerate <- function(x,
                          y,
                          n,
                          p,
                          intercept,
                          family,
                          rare,
                          hap,
                          region_ind,
                          forced,
                          p_forced,
                          inform,
                          cov,
                          p_cov,
                          a1,
                          which_ind)
{

    # create all possible models
    num_models <- 2^p
    all_models <- t(expand.grid(lapply(vector("list", p), function(v) {c(FALSE, TRUE)})))

    # initialize glm parameters
    control <- do.call("glm.control", list())
    family_func <- do.call(family, list())
    weights <- rep(1.0, n)
    offset <- rep(0.0, n)
    if (family == "binomial") {
        family_func$mu.eta <- logit_mu_eta
        mustart <- (weights * y + 0.5) / (weights + 1)
        m <- weights * y
        if (any(abs(m - round(m)) > 0.001))
            warning("non-integer #successes in a binomial glm!")
    } else {
        family_func$mu.eta <- identity_mu_eta
        mustart <- y
    }

    # compute null deviance
    null_dev <- glm_fit_custom(x = forced,
                               y = y,
                               nobs = n,
                               nvars = intercept + p_forced,
                               family = family_func,
                               control = control,
                               weights = weights,
                               mustart = mustart,
                               m = m,
                               offset = offset)$dev
  
    # setup external data
    if (inform) {
        a0 <- qnorm((1 - 2^(-1 / nrow(cov))))
        mu <- a0 + as.matrix(cov) %*% as.matrix(a1)
        lprob_inc <- pnorm(0, mean = mu, lower.tail = FALSE, log.p = TRUE)
        lprob_ninc <- pnorm(0, mean = mu, lower.tail = TRUE, log.p = TRUE)
    }

    # initialize objects to hold results
    ll <- rep(NA, num_models)
    fitness <- rep(NA, num_models)
    logPrM <- rep(NA, num_models)
    active <- apply(all_models, 2, function(x) paste(which(x), collapse = "-"))
    active[1] <- "0"
    coef <- matrix(0, nrow = num_models, ncol = length(which_ind))
    alpha <- rep(a1, num_models)

    # fit all models
    for (i in seq.int(num_models)) {

        # get active variables
        z_current <- all_models[, i, drop = FALSE]
        num_active <- sum(z_current)

        # fit model
        fit_glm <- bvs_fit(z = z_current,
                           num_active = num_active,
                           y = y,
                           x = x,
                           n = n,
                           p = p,
                           intercept = intercept,
                           rare = rare,
                           hap = hap,
                           region_ind = region_ind,
                           forced = forced,
                           p_forced = p_forced,
                           family_func = family_func,
                           control = control,
                           weights = weights,
                           offset = offset,
                           mustart = mustart,
                           m = m)

        # get coef vector
        if (fit_glm$num_vars > 0) {
            if (rare) {
                coef[i, which_ind] <- fit_glm$coef
            } else {
                coef[i, z_current] <- fit_glm$coef
            }
        }

        # compute log-likelihood
        ll[i] <- 0.5 * fit_glm$deviance + fit_glm$num_vars

        # calculate prior on model
        # If inform ==TRUE use probit specification
        # and if mult.regions==TRUE use region level probit specification
        if (inform) {
            logPrM[i] <- sum(lprob_inc[z_current]) + sum(lprob_ninc[!z_current])
        } else {
            logPrM[i] <- BetaBinomial(p = p, pgamma = num_active)
        }

        # calculate fitness = ll - logPrM
        fitness[i] <- ll[i] - logPrM[i]
    }

    models_accepted <- list(model_id = 1:num_models,
                            active = active,
                            ll = ll,
                            coef = coef)

    results <- list(fitness = fitness,
                    logPrM = logPrM,
                    model_path = 1:num_models,
                    alpha = alpha,
                    models_accepted = models_accepted,
                    null_dev = null_dev)
}
