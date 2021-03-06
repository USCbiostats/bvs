#' Calculates Posterior Summaries for BVS Methods
#'
#' @description This function calculates global, regional, and marginal Bayes Factors that give the strength
#'  of evidence of there being an association in the overall set of variants of interest,
#'  the individual genes of interest (if specified) and the individual variants of interest.
#'
#' @param object an object of class 'bvs'
#' @param burnin number of burn-in interations. Automatically set to zero if 'enumerate' method used.
#' @param prior_cov matrix of predictor-level prior covariates (required only if used in original model fit)
#' @param ... additional arguments as required by summary S3 object

#' @export
summary.bvs <- function(object, burnin = 1000, prior_cov = NULL, ...)
{

    # set burnin to zero if method is enumerate
    if (object$model_info$method == "enumerate" && burnin > 0) {
        burnin <- 0
    }

    # check whether prior data needed
    if (object$model_info$inform && is.null(prior_cov)) {
        stop("'bvs' object has inform = TRUE, but prior_cov = NULL. Please provide prior data matrix.")
    }

    # subset trials based on burnin
    if (burnin > 0) {
        models_sub <- which(object$models_accepted$model_id %in% object$model_path[-c(1:burnin)])
        a1 <- t(object$alpha[-c(1:burnin), , drop = FALSE])
    } else {
        models_sub <- 1:length(object$models_accepted$model_id)
        a1 <- t(object$alpha)
    }

    # get unique values for all accepted models / sort into order first accepted
    model_id <- object$models_accepted$model_id[models_sub]
    active <- object$models_accepted$active[models_sub]
    num_active <- object$models_accepted$num_active[models_sub]
    loglike <- object$models_accepted$loglike[models_sub]
    coef <- t(object$models_accepted$coef[models_sub, , drop = FALSE])
    logfitness <- object$logfitness[match(model_id, object$model_path)]
    logPrM <- object$logPrM[match(model_id, object$model_path)]

    # parameters
    p <- object$model_info$nvars
    p_cov <- object$model_info$nvars_prior
    alpha_bb <- object$model_info$prior_model$alpha
    beta_bb <- object$model_info$prior_model$beta
    num_coef <- nrow(coef)
    a0 <- qnorm(1 - 2^(-1 / p))

    # make sure null model is in the results
    null_ind <- which(active == "0")
    if (length(null_ind) == 0) {
        model_id <- c(0, model_id)
        active <- c("0", active)
        num_active <- c(0, num_active)
        loglike <- c(object$model_info$nullLogLike, loglike)
        coef <- cbind(0, coef)
        if (object$model_info$inform) {
            logPrM <- c(sum(pnorm(0, mean = rep(a0, p), lower.tail = TRUE, log.p = TRUE)), logPrM)
        } else {
            logPrM <- c(computeModelPrior(0, p, alpha_bb, beta_bb), logPrM)
        }
        loglike <- c(0.5 * object$model_info$null_dev, loglike)
        logfitness <- c(object$model_info$nullLogLike + logPrM[1], logfitness)
        null_ind <- 1
    }

    # Post expectation of alpha
    post_alpha <- rowMeans(a1)

    # If inform == TRUE average prior inclusion probability across multiple values of alpha!
    if (object$model_info$inform) {
        eta <- a0 + prior_cov %*% a1
        prior.marg <- apply(eta, 1, function(x) mean(pnorm(0, mean = x, lower.tail = FALSE)))
        lprob.inc = log(prior.marg)
        lprob.ninc = log(1 - prior.marg)
        null.prob.ninc = pnorm(0, mean = a0, lower.tail = TRUE, log.p = TRUE)
        null.prior.marg = 1 - exp(null.prob.ninc)

        # Global Priors
        prior.null = exp(sum(lprob.ninc))
        prior.alt = 1 - prior.null
    } else {
        null.prior.marg <- 1 / (p + 1)
        prior.null <- 0.5
        prior.alt <- 0.5
    }

    active_mat <- matrix(0, nrow = length(active), ncol = p)
    u_active <- lapply(active, function(x) as.integer(strsplit(x, "-")[[1]]))
    for (i in 1:nrow(active_mat)) {
        active_mat[i, u_active[[i]]] <- 1
    }

    if (object$model_info$inform) {
        new_lPrM <- apply(active_mat, 1, function(x) {sum(lprob.inc[x == 1]) + sum(lprob.ninc[x == 0])})
        logfitness <- loglike + new_lPrM
    } else {
        new_lPrM <- logPrM
    }

    # Calculate Posterior model probabilities
    PrMgivenD <- exp(logfitness - min(logfitness)) / sum(exp(logfitness - min(logfitness)))

    ##Global Posterior Prob & BF
    post.null <- PrMgivenD[null_ind]
    post.alt <- 1 - post.null
    global_bf <- post.alt / post.null

    # Marginal Inclusion probabilities
    post.marg <- drop(crossprod(active_mat, PrMgivenD))
    post.marg.n <- drop(crossprod(!active_mat, PrMgivenD))
    global_marg_bf <- exp(log(post.marg) - log(post.marg.n) - log(null.prior.marg) +  log(1 - null.prior.marg))

    # Gene Inclusion probabilities
    active_region <- NULL
    region_marg_bf <- NULL
    region.post.inc <- NULL
    if (object$model_info$num_regions > 1) {
        regions <- as.factor(object$model_info$regions)
        regions_mat <- model.matrix(~ regions - 1)
        active_region <- (active_mat %*% regions_mat) > 0
        active_region[] <- as.numeric(active_region)
        colnames(active_region) <- paste("region_", levels(regions), sep = "")
        probne0.r <- drop(PrMgivenD %*% active_region)
        probe0.r <- drop(PrMgivenD %*% !active_region)
        null.prior.marg.r <- (1 - null.prior.marg)^apply(regions_mat, 2, sum)
        region.post.inc <- probne0.r
        region_marg_bf <- exp(log(probne0.r) - log(probe0.r) - log(1 - null.prior.marg.r) + log(null.prior.marg.r))
    }

    # Posterior estimates of coefficients and sd
    post_coef <- drop(coef %*% PrMgivenD)
    names(post_coef) <- rownames(coef)

    # Add variable names
    names(global_marg_bf) <- object$model_info$varnames
    colnames(active_mat) <- object$model_info$varnames

    ##Save Results
    results <- list("global_bf" = global_bf,
                    "marg_bf" = global_marg_bf,
                    "post_coef" = post_coef,
                    region_level = list("region_bf" = region_marg_bf,
                                        "active_region" = active_region),
                    model_level = list("model_id" = model_id,
                                       "active" = active,
                                       "num_active" = num_active,
                                       "loglike" = loglike,
                                       "prior_prob" = exp(new_lPrM),
                                       "post_prob" = PrMgivenD,
                                       "coef" = t(coef),
                                       "active_mat" = active_mat),
                    "post_alpha" = post_alpha,
                    "model_info" = object$model_info)
    class(results) <- "summary.bvs"
    return(results)
}
