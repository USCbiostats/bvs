
bvs_sample <- function(x,
                       y,
                       n,
                       p,
                       family,
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
                       iter,
                       outfile,
                       status_file,
                       old_results)
{

    # initialize parameters
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
    if (inform) {
        a0 <- qnorm(1 - 2^(-1 / p))
    }

    # compute null deviance
    null_dev <- sum(family_func$dev.resids(y, mean(y), weights))

    # check if old results present
    if (length(old_results) > 0) {
        iter <- iter - ncol(old_results)
    } else {
        iter <- iter + 1
        z_current <- rep(FALSE, p)
        z_current[sample(1:p, 5)] <- TRUE
        if (inform) {
            a1_current <- rep(0, p_cov)
        }
    }

    # initialize results objects
    model_path <- c(1, rep(NA, iter - 1))
    coef <- matrix(0, nrow = iter, ncol = length(which_ind))
    ll <- rep(NA, iter)
    fitness <- rep(NA, iter)
    logPrM <- rep(NA, iter)
    alpha <- matrix(0, nrow = iter, ncol = max(1, p_cov))
    active <- rep(NA, iter)

    # initialize unordered map to hold all unique previous models
    models_fit <- new_table(1)
    set_element_in_table(active[1] <- paste0(which(z_current), collapse = "-"), 1, models_fit)

    # calculate current fitness
    lower <- rep(0, p)
    upper <- rep(Inf, p)
    lower[!z_current] <- -Inf
    upper[!z_current] <- 0
    t_current <- rtnorm(p, lower = lower, upper = upper)

    # fit initial model
    num_active <- sum(z_current)

    fit_glm <- bvs_fit(z = z_current,
                       num_active = num_active,
                       y = y,
                       x = x,
                       n = n,
                       p = p,
                       rare = rare,
                       hap = hap,
                       region_ind = region_ind,
                       num_regions = num_regions,
                       forced = forced,
                       p_forced = p_forced,
                       family_func = family_func,
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

    # compute log likelihood
    ll[1] <- 0.5 * fit_glm$deviance + fit_glm$num_vars

    # calculate prior on model
    if (inform) {
        v_hat <- solve(diag(1, p_cov) + crossprod(cov))
        mu <- a0 + cov %*% a1_current
        lprob_inc <- pnorm(0, mean = mu, lower.tail = FALSE, log.p = TRUE)
        lprob_ninc <- pnorm(0, mean = mu, lower.tail = TRUE, log.p = TRUE)
        logPrM[1] <- logPrM_current <- sum(lprob_inc[z_current]) + sum(lprob_ninc[!z_current])
    } else {
        logPrM[1] <- logPrM_current <- BetaBinomial(p = p, pgamma = num_active)
    }

    # calculate fitness = ll - logPrM
    fitness[1] <- fitness_current <- ll[1] - logPrM_current

    # loop through all trials
    pb <- txtProgressBar(min = 2, max = iter, style = 3)
    i <- 2
    idx <- 2

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

            # calculate new fitness | a1
            fitness_current <- fitness_current + logPrM[i - 1] - logPrM_current
        }

        # sample model
        mutate_ind <- ceiling(runif(1) * p)
        testing <- runif(1) # remove when done testing
        z_current[mutate_ind] <- !z_current[mutate_ind]
        num_active <- sum(z_current)

        # check if model already fit
        if (num_active > 0) {
            active_new <- paste0(which(z_current), collapse = "-")
        } else {
            active_new <- "0"
        }
        model_ind <- get_element_from_table(active_new, models_fit)

        if (model_ind != 0) {
            model_propose <- model_ind
        } else {
            active[idx] <- active_new
            set_element_in_table(active_new, idx, models_fit)
            model_propose <- idx

            # fit model
            fit_glm <- bvs_fit(z = z_current,
                               num_active = num_active,
                               y = y,
                               x = x,
                               n = n,
                               p = p,
                               rare = rare,
                               hap = hap,
                               region_ind = region_ind,
                               num_regions = num_regions,
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
                    coef[idx, which_ind] <- fit_glm$coef
                } else {
                    coef[idx, z_current] <- fit_glm$coef
                }
            }
            # compute log-likelihood
            ll[idx] <- 0.5 * fit_glm$deviance + fit_glm$num_vars
            idx <- idx + 1
        }

        # calculate prior on model
        if (inform) {
            logPrM_new <- sum(lprob_inc[z_current]) + sum(lprob_ninc[!z_current])
        } else {
            logPrM_new <- BetaBinomial(p = p, pgamma = num_active)
        }

        # calculate fitness = ll - logPrM
        fitness_new <- ll[model_propose] - logPrM_new

        #If fit ratio > runif(1), accept new model and z, else keep old results
        accept <- runif(1) <= min(1, 1 / exp(fitness_new - fitness_current))
        if (accept) {
            model_path[i] <- model_propose
            fitness[i] <- fitness_current <- fitness_new
            logPrM[i] <- logPrM_current <- logPrM_new
        } else {
            model_path[i] <- model_path[i - 1]
            fitness[i] <- fitness_current
            logPrM[i] <- logPrM_current
            z_current[mutate_ind] <- !z_current[mutate_ind]
        }

        if (inform) {
            alpha[i, ] <- a1_current
        }

        i <- i + 1
        # Update progress
        setTxtProgressBar(pb, i)
    }
    cat("\nDone\n")

    models_chosen <- unique(model_path)
    models_accepted <- list(model_id = models_chosen,
                           active = active[models_chosen],
                           ll = ll[models_chosen],
                           coef = coef[models_chosen, , drop = FALSE])

    results <- list(fitness = fitness,
                    logPrM = logPrM,
                    model_path = model_path,
                    alpha = alpha,
                    models_accepted = models_accepted,
                    null_dev = null_dev)
}







