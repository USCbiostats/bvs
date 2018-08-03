bvs_sample <- function(x,
                       y,
                       n,
                       p,
                       intercept,
                       family,
                       prior_coef,
                       rare,
                       hap,
                       region_ind,
                       num_regions,
                       forced,
                       p_forced,
                       inform,
                       cov,
                       p_cov,
                       which_ind,
                       iter)
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
            warning("non-integer # of successes in a binomial glm!")
    } else {
        family_func$mu.eta <- identity_mu_eta
        mustart <- y
    }

    # compute null marginal log like
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

    # initialize hash table to hold indices of all unique previous models
    models_fit <- new_table(1)
    iter <- iter + 1 # remove when done testing

    # initialize results objects
    logfitness <- rep(NA, iter)
    logPrM <- rep(NA, iter)
    model_path <- c(1, rep(NA, iter - 1))
    alpha <- matrix(0, nrow = iter, ncol = max(1, p_cov))
    active <- rep(NA, iter)
    num_active <- rep(NA, iter)
    coef <- matrix(0, nrow = iter, ncol = length(which_ind))
    loglike <- rep(NA, iter)

    # initialize vector that indicates currently active variables
    z_current <- rep(FALSE, p)
    z_current[sample(1:p, 5)] <- TRUE
    set_element_in_table(key = active[1] <- paste0(which(z_current), collapse = "-"),
                         value = 1,
                         table = models_fit)

    # fit initial model
    num_active[1] <- sum(z_current)

    fit_glm <- bvs_fit(z = z_current,
                       num_active = num_active[1],
                       y = y,
                       x = x,
                       n = n,
                       p = p,
                       intercept = intercept,
                       rare = rare,
                       hap = hap,
                       region_ind = region_ind,
                       num_regions = num_regions,
                       forced = forced,
                       p_forced = p_forced,
                       family_func = family_func,
                       prior_coef = prior_coef,
                       control = control,
                       weights = weights,
                       offset = offset,
                       mustart = mustart,
                       m = m)

    # get intial coef vector
    if (fit_glm$num_vars > 0) {
        if (rare) {
            coef[1, which_ind] <- fit_glm$coef
        } else {
            coef[1, z_current] <- fit_glm$coef
        }
    }

    # get marginal log likelihood
    loglike[1] <- fit_glm$marg_ll

    # initialize t_current
    lower <- rep(0, p)
    upper <- rep(Inf, p)
    lower[!z_current] <- -Inf
    upper[!z_current] <- 0
    t_current <- rtnorm(p, lower = lower, upper = upper)

    # calculate prior on model
    if (inform) {
        a0 <- qnorm(1 - 2^(-1 / p))
        v_hat <- solve(diag(1, p_cov) + crossprod(cov))
        a1_current <- rep(0, p_cov)
        mu <- a0 + cov %*% a1_current
        lprob_inc <- pnorm(0, mean = mu, lower.tail = FALSE, log.p = TRUE)
        lprob_ninc <- pnorm(0, mean = mu, lower.tail = TRUE, log.p = TRUE)
        logPrM[1] <- logPrM_current <- sum(lprob_inc[z_current]) + sum(lprob_ninc[!z_current])
    } else {
        logPrM[1] <- logPrM_current <- logBetaBinomial(p = p, pgamma = num_active[1])
    }

    # calculate log(fitness) = marginal_ll - logPrM
    logfitness[1] <- fitness_current <- loglike[1] + logPrM_current

    # loop through all trials, track progress with progress bar
    i <- 2
    idx <- 2
    pb <- txtProgressBar(min = i, max = iter, style = 3)

    while (i <= iter) {
        # sample new a1 if using prior or keep a1 = 0
        if (inform) {
            alpha_hat <- v_hat %*% crossprod(cov, t_current - a0)
            a1_current <- mvrnorm(1, mu = alpha_hat, Sigma = v_hat)

            # sample latent variable | current a1 and current model
            lower <- rep(0, p)
            upper <- rep(Inf, p)
            lower[!z_current] <- -Inf
            upper[!z_current] <- 0
            mu <- a0 + cov %*% a1_current
            t_current <- rtnorm(p, mean = mu, lower = lower, upper = upper)

            # calculate new prior
            lprob_inc <- pnorm(0, mean = mu, lower.tail = FALSE, log.p = TRUE)
            lprob_ninc <- pnorm(0, mean = mu, lower.tail = TRUE, log.p = TRUE)
            logPrM_current <- sum(lprob_inc[z_current]) + sum(lprob_ninc[!z_current])

            # update fitness | a1 for currently selected model
            fitness_current <- fitness_current - logPrM[i - 1] + logPrM_current
        }

        # sample model
        mutate_ind <- ceiling(runif(1) * p)
        testing <- runif(1) # remove when done testing
        z_current[mutate_ind] <- !z_current[mutate_ind]
        num_active_new <- sum(z_current)

        # check if model already fit
        if (num_active_new > 0) {
            active_new <- paste0(which(z_current), collapse = "-")
        } else {
            active_new <- "0"
        }
        model_ind <- get_element_from_table(active_new, models_fit)

        if (model_ind != 0) {
            model_propose <- model_ind
        } else {
            active[idx] <- active_new
            num_active[idx] <- num_active_new
            set_element_in_table(active_new, idx, models_fit)
            model_propose <- idx

            # fit model
            fit_glm <- bvs_fit(z = z_current,
                               num_active = num_active_new,
                               y = y,
                               x = x,
                               n = n,
                               p = p,
                               intercept = intercept,
                               rare = rare,
                               hap = hap,
                               region_ind = region_ind,
                               num_regions = num_regions,
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
            # get marginal log likelihood
            loglike[idx] <- fit_glm$marg_ll
            idx <- idx + 1
        }

        # calculate prior on model
        if (inform) {
            logPrM_new <- sum(lprob_inc[z_current]) + sum(lprob_ninc[!z_current])
        } else {
            logPrM_new <- logBetaBinomial(p = p, pgamma = num_active_new)
        }

        # calculate fitness = loglike + logPrM
        fitness_new <- loglike[model_propose] + logPrM_new

        #If fit ratio > runif(1), accept new model and z, else keep current model
        accept <- runif(1) <= min(1, exp(fitness_new - fitness_current))
        if (accept) {
            model_path[i] <- model_propose
            logfitness[i] <- fitness_current <- fitness_new
            logPrM[i] <- logPrM_current <- logPrM_new
        } else {
            model_path[i] <- model_path[i - 1]
            logfitness[i] <- fitness_current
            logPrM[i] <- logPrM_current
            z_current[mutate_ind] <- !z_current[mutate_ind]
        }

        if (inform) {
            alpha[i, ] <- a1_current
        }

        # Update progress
        i <- i + 1
        setTxtProgressBar(pb, i)
    }
    cat("\nDone\n")

    models_chosen <- unique(model_path)
    models_accepted <- list(model_id = models_chosen,
                            active = active[models_chosen],
                            num_active = num_active[models_chosen],
                            loglike = loglike[models_chosen],
                            coef = coef[models_chosen, , drop = FALSE])

    rownames(models_accepted$coef) <- models_chosen
    results <- list(logfitness = logfitness,
                    logPrM = logPrM,
                    model_path = model_path,
                    alpha = alpha,
                    models_accepted = models_accepted,
                    nullLogLike = nullLogLike)
}







