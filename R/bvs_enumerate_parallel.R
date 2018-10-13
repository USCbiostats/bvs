bvs_enumerate_parallel <- function(x,
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


    active <- "0"
    num_active <- 0
    z_current <- rep(0, p)
    if (inform) {
        logPrM <- sum(lprob_inc[z_current]) + sum(lprob_ninc[!z_current])
    } else {
        logPrM <- computeModelPrior(0, p, alpha_bb, beta_bb)
    }
    loglike <- nullLogLike
    logfitness <- nullLogLike + logPrM
    coef <- matrix(0, nrow = length(which_ind), ncol = 1)
    alpha <- rep(a1, num_models)

    # loop through all models
    if (num_models > 1) {
        for (k in 1L:maxk) {

            # get indices for all unique models of size k
            models_size_k <- combn(1:p, k)
            active <- c(active, apply(models_size_k, 2, function(x) paste0(x, collapse = "-")))
            num_active <- c(num_active, rep(k, NCOL(models_size_k)))

            # fit all models of size k
            resultk <- foreach(i = iter(models_size_k, by = "col"), .combine = 'cbind') %dopar% {

                z_current <- 1:p %in% i

                # fit model
                fit_glm <- bvs_fit(z = z_current,
                                   num_active = k,
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

                # get coef
                if (rare) {
                    coef_model <- rep(0.0, length(which_ind))
                    coef_model[which_ind] <- fit_glm$coef
                } else {
                    coef_model <- rep(0.0, p)
                    coef_model[z_current] <- fit_glm$coef
                }

                # get marginal log-likelihood
                loglike_model <- fit_glm$marg_ll

                # calculate prior on model
                # If inform ==TRUE use probit specification
                # and if mult.regions==TRUE use region level probit specification
                if (inform) {
                    logPrM_model <- sum(lprob_inc[z_current]) + sum(lprob_ninc[!z_current])
                } else {
                    logPrM_model <- computeModelPrior(k, p, alpha_bb, beta_bb)
                }

                # calculate logfitness = loglike - logPrM
                logfitness_model <- loglike_model + logPrM_model

                c(logfitness_model, loglike_model, logPrM_model, coef_model)
            }

            logfitness <- c(logfitness, unname(resultk[1, ]))
            loglike <- c(loglike, unname(resultk[2, ]))
            logPrM <- c(logPrM, unname(resultk[3, ]))
            coef <- cbind(coef, unname(resultk[-c(1:3), , drop = FALSE]))
        }
    }

    coef <- t(coef)
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
