bvs_enumerate <- function(x,
                          y,
                          n,
                          p,
                          intercept,
                          family,
                          prior_model,
                          prior_coef,
                          rare,
                          hap,
                          region_ind,
                          forced,
                          p_forced,
                          maxk,
                          inform,
                          cov,
                          p_cov,
                          a1,
                          which_ind)
{

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

    # initialize model size beta-binomial parameters
    alpha_bb <- prior_model[[1]]
    beta_bb <- prior_model[[2]]

    # setup external data
    if (inform) {
        a0 <- qnorm((1 - 2^(-1 / nrow(cov))))
        mu <- a0 + as.matrix(cov) %*% as.matrix(a1)
        lprob_inc <- pnorm(0, mean = mu, lower.tail = FALSE, log.p = TRUE)
        lprob_ninc <- pnorm(0, mean = mu, lower.tail = TRUE, log.p = TRUE)
    }

    # initialize objects to hold results
    num_models <- sum(sapply(1:maxk, function(x) choose(p, x))) + 1
    if (num_models > 100000) {
        continue_analysis <- menu(c("Yes", "No"), title = paste0("Warning: Number of models (", num_models,") to compute is very large. Do you want to continue?"))
        if(continue_analysis == 2) {
            return(NULL)
        }
    }
    loglike <- rep(NA, num_models)
    logfitness <- rep(NA, num_models)
    logPrM <- rep(NA, num_models)
    active <- rep(NA, num_models)
    num_active <- rep(NA, num_models)
    coef <- matrix(0, nrow = num_models, ncol = length(which_ind))
    alpha <- rep(a1, num_models)

    # compute null model
    nullLogLike <- bvs_fit(z = NULL,
                           num_active = 0L,
                           y = y,
                           x = NULL,
                           n = n,
                           p = p,
                           intercept = intercept,
                           rare = rare,
                           hap = hap,
                           region_ind = region_ind,
                           forced = forced,
                           p_forced = p_forced,
                           family_func = family_func,
                           prior_coef = prior_coef,
                           control = control,
                           weights = weights,
                           offset = offset,
                           mustart = mustart,
                           m = m)$marg_ll


    active[1] <- "0"
    num_active[1] <- 0
    z_current <- rep(0, p)
    if (inform) {
        logPrM[1] <- sum(lprob_inc[z_current]) + sum(lprob_ninc[!z_current])
    } else {
        logPrM[1] <- computeModelPrior(0, p, alpha_bb, beta_bb)
    }
    loglike[1] <- nullLogLike
    logfitness[1] <- nullLogLike + logPrM[1]

    # loop through all models
    idx <- 2
    pb <- txtProgressBar(min = idx, max = num_models, style = 3)
    for (k in 1L:maxk) {
        # get indices for all unique models of size k
        models_size_k <- combn(1:p, k)

        # fit all models of size k
        for (i in 1L:NCOL(models_size_k)) {
            z_current <- 1:p %in% models_size_k[, i]
            active[idx] <- paste0(models_size_k[i], collapse = "-")
            num_active[idx] <- sum(z_current)

            # fit model
            fit_glm <- bvs_fit(z = z_current,
                               num_active = num_active[idx],
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
                               prior_coef = prior_coef,
                               control = control,
                               weights = weights,
                               offset = offset,
                               mustart = mustart,
                               m = m)

            # get coef vector
            if (fit_glm$num_vars > 0) {
                if (rare) {
                    coef[idx, which_ind] <- fit_glm$coef
                } else {
                    coef[idx, z_current] <- fit_glm$coef
                }
            }

            # get marginal log-likelihood
            loglike[idx] <- fit_glm$marg_ll

            # calculate prior on model
            # If inform ==TRUE use probit specification
            # and if mult.regions==TRUE use region level probit specification
            if (inform) {
                logPrM[idx] <- sum(lprob_inc[z_current]) + sum(lprob_ninc[!z_current])
            } else {
                logPrM[idx] <- computeModelPrior(num_active[idx], p, alpha_bb, beta_bb)
            }

            # calculate logfitness = loglike - logPrM
            logfitness[idx] <- loglike[idx] + logPrM[idx]

            # update progress
            idx <- idx + 1
            setTxtProgressBar(pb, idx)
        }

    }

    rownames(coef) <- 1:num_models
    models_accepted <- list(model_id = 1:num_models,
                            active = active,
                            num_active = num_active,
                            loglike = loglike,
                            coef = coef)

    results <- list(logfitness = logfitness,
                    logPrM = logPrM,
                    model_path = 1:num_models,
                    alpha = alpha,
                    models_accepted = models_accepted,
                    nullLogLike = nullLogLike)
}
