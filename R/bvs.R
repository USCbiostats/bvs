#' @useDynLib bvs, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Bayesian Variable Selection
#' @param y outcome variable
#' @param x predictor design matrix
#' @param forced (optional) \eqn{n x c} matrix of \emph{c} confounding variables that one wishes to adjust the
#' analysis for and that will be forced into every model.
#' @param family specifies error distribution for outcome variable, options include
#' \itemize{
#'    \item "gaussian"
#'    \item "binomial"
#' }
#' @param method specifies how to search the model space
#' \itemize{
#'    \item "enumerate": computes and summarizes all possible models in model space. Not recommended for problems where \eqn{p > 20}.
#'    \item "sample": performs basic basic Metropolis-Hastings algorithm to sample models from model space. For informative
#'    marginal inclusion probabilities, the algorithm also performs basic MCMC algorithm to sample the effects of predictor-level
#'    covariates (alpha).
#' }
#' @param rare if rare = TRUE, corresponds to the Bayesian Risk Index (BRI) algorithm of Quintana and Conti (2011) that constructs
#' a risk index based on the multiple rare variants within each model. The marginal likelihood of each model is then calculated
#' based on the corresponding risk index.
#' @param regions (optional) \eqn{p x 1} character or factor vector that identifies a user-defined region for each variant. If
#' rare = TRUE, then multiple region-specific risk indices are computed for each model.
#' @param prior_cov (optional) if method = "sample", a \eqn{p x q} matrix of \emph{q} predictor-level covariates that the user
#' wishes to incorporate into the estimation of the marginal inclusion probabilities using the iBMU algorithm.
#' @param a1 (optional) if method = "enumerate", a \eqn{q x 1} vector of specified effects of each predictor-level covariate.
#' @param hap (not yet implemented) if hap = TRUE, esimtate a set of haplotypes from the multiple variants within each moel and the marginal likelihood
#' of each model is calculated based on the set of haplotypes.
#' @param iter if method = "sample", the number of iterations to run the algorithm.
#' @param save_iter if method = "sample", the number of iterations between each checkpoint.  A checkpoint file is written
#' every save.iter iterations
#' @param outfile if method = "sample", character string giving the pathname of the checkpoint file to save the output of the algorithm.
#' @param status_file if method = "sample", character string giving the pathname of the file to write the status of the algorithm.
#' @param old_results if method = "sample", old output from sampleBVS that has been run for a subset of the total number of
#' iterations that the user wanted to run. if specified the sampling algorithm will start from the last sampled model in old.results.
#' To be used if sampleBVS has been interrupted for some reason.
#' @param control specifies 'bvs' control object.
#'
#' @import stats haplo.stats
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom MASS mvrnorm
#' @importFrom msm rtnorm

#' @export
bvs <- function(y,
                x,
                forced = NULL,
                family = c("gaussian", "binomial"),
                method = c("enumerate", "sample"),
                rare = FALSE,
                regions = NULL,
                prior_cov = NULL,
                a1 = 0,
                hap = FALSE,
                iter = 10000,
                save_iter = 0,
                outfile = NULL,
                status_file = NULL,
                old_results = NULL,
                control = list())
{

    # check hap
    if (hap) {
        stop("Error: haplotypes not yet implemented.")
    }

    # check family/method
    family <- match.arg(family)
    method <- match.arg(method)

    # check y
    y <- drop(y)
    if (family == "binomial") {
        if (is.factor(y))
            y <- y != levels(y)[1L]
        if (any(y < 0 | y > 1))
            stop("y values must be 0 <= y <= 1")
    }

    # check x and dim x/y
    n <- nrow(x)
    p <- ncol(x)

    y_len <- if (is.null(dim(y)))
                length(y)
            else
                dim(y)[1]

    if (n != y_len)
        stop(paste("Length of y (", y_len, ") not equal to number of rows of x (", n, ")", sep = ""))
    if (class(x) != "matrix")
        x <- as.matrix(x)
    if (!(typeof(x) %in% c("double", "numeric", "integer")))
        stop("x contains non-numeric values")

    # check regions / rare
    if (!is.null(regions)) {
        if (length(regions) != p)
            stop("Length of regions (", length(regions), ") not equal to number of variables (p)")
        num_regions <- length(unique(regions))
    } else {
        regions <- rep(1, p)
        num_regions <- 1
    }

    if (rare) {
        which_ind <- 1:num_regions
        if (num_regions > 1) {
            regions <- as.factor(regions)
            region_ind <- model.matrix( ~ regions - 1)
            colnames(region_ind) <- paste("region_", levels(regions), sep = "")
        } else {
            region_ind <- matrix(1, nrow = p)
            colnames(region_ind) <- "risk_index"
        }
    } else {
        which_ind <- 1:p
    }

    # check forced
    if (!is.null(forced)) {
        n_forced <- nrow(forced)
        p_forced <- ncol(forced)
        if (n_forced != y_len)
            stop(paste("Length of y (", y_len, ") not equal to number of rows of forced (", n_forced, ")", sep = ""))
        if (class(forced) != "matrix")
            forced <- as.matrix(forced)
        if (!(typeof(forced) %in% c("double", "numeric", "integer")))
            stop("forced contains non-numeric values")
    } else {
        p_forced <- 0
    }

    # check prior_cov
    if (!is.null(prior_cov)) {
        inform <- TRUE
        n_cov <- nrow(prior_cov)
        p_cov <- ncol(prior_cov)
        if (n_cov != p)
            stop(paste("Number of rows of prior_cov (", n_cov, ") not equal to number of columns of x (", p, ")", sep = ""))
        if (class(prior_cov) != "matrix")
            prior_cov <- as.matrix(prior_cov)
        if (!(typeof(prior_cov) %in% c("double", "numeric", "integer")))
            stop("prior_cov contains non-numeric values")
    } else {
        inform <- FALSE
        p_cov <- 0
    }

    # fit models by enumerating all combinations or M-H
    fit <- switch(method,
                  enumerate = bvs_enumerate(x, y, n, p, family, rare, hap, region_ind, forced,
                                            p_forced, inform, prior_cov, p_cov, a1, which_ind),

                  sample = bvs_sample(x, y, n, p, family, rare, hap, region_ind, num_regions, forced,
                                      p_forced, inform, prior_cov, p_cov, which_ind, iter, save_iter,
                                      status_file, old_results)
    )

    # add name attributes
    if (!is.null(colnames(x))) {
        varnames <- colnames(x)
    } else {
        varnames <- paste("coef", 1:ncol(x), sep = "")
    }

    if (rare) {
        colnames(fit$models_accepted$coef) <- colnames(region_ind)
    } else {
        colnames(fit$models_accepted$coef) <- varnames
    }

    if (inform) {
        if (!is.null(colnames(prior_cov))) {
            colnames(fit$alpha) <- colnames(prior_cov)
        } else {
            colnames(fit$alpha) <- paste("alpha", 1:ncol(prior_cov), sep = "")
        }
    }

    fit$model_info <- list(method = method,
                           nobs = n,
                           nvars = p,
                           varnames = varnames,
                           rare = rare,
                           regions = regions,
                           num_regions = num_regions,
                           inform = inform,
                           nvars_prior = p_cov,
                           hap = hap,
                           null_dev = fit$null_dev)
    fit$null_dev <- NULL
    fit$model.info
    class(fit) <- "bvs"
    return(fit)
}

