#' @useDynLib bvs, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Bayesian Variable Selection
#' @param y outcome variable
#' @param x predictor design matrix
#' @param forced (optional) \eqn{n x c} matrix of \emph{c} confounding variables that one wishes to adjust the
#' analysis for and that will be forced into every model.
#' @param intercept indicates whether models should include an intercept. Default = TRUE.
#' @param family specifies error distribution for outcome variable, options include
#' \itemize{
#'    \item "gaussian" (default)
#'    \item "binomial"
#' }
#' @param method specifies how to search the model space
#' \itemize{
#'    \item "sample" (default) : performs basic Metropolis-Hastings algorithm to sample models from model space. For informative
#'    marginal inclusion probabilities, the algorithm also performs basic MCMC algorithm to sample the effects of predictor-level
#'    covariates (alpha).
#'    \item "enumerate": computes and summarizes all possible models in model space. Not recommended for problems where \eqn{p > 20}.
#' }
#' @param prior_model specifies parameters for beta-binomial prior on model size. To specify, pass a list with the following elements
#' \itemize{
#'    \item \eqn{alpha} = numeric value for first shape parameter (default = 1)
#'    \item \eqn{beta} = numeric value for second shape parameter (default = p)
#' }
#' Example: \code{list(alpha = 1, beta = 2)}
#' @param prior_coef specifies prior for regression coefficients (only for use when family = "gaussian").
#' \itemize{
#'    \item "none" (default)
#'    \item "gprior": To specify the parameters for the g-prior, pass a list object with the following elements
#'    \itemize{
#'    \item "gprior"
#'    \item \eqn{g} = numeric value for sparsity parameter, recommended default = \eqn{max(n, p^2)}
#'    \item \eqn{alpha} = numeric value that specifies variance in g-prior, recommended default = 0.01
#'    \item \eqn{beta} = numeric value that specifies variance in g-prior, recommended default = 0.01
#'    }
#'    To use the default values, pass the string "default" for any parameter.
#'
#'    Example: \code{list("gprior", "default", 0.02, 0.02)}
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
#' @param control specifies 'bvs' control object.
#'
#' @import stats haplo.stats
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom MASS mvrnorm
#' @importFrom msm rtnorm
#' @importFrom TailRank dbb

#' @export
bvs <- function(y,
                x,
                forced = NULL,
                intercept = TRUE,
                family = c("gaussian", "binomial"),
                method = c("sample", "enumerate"),
                prior_model = list(alpha = 1, beta = p),
                prior_coef = list("none"),
                rare = FALSE,
                regions = NULL,
                prior_cov = NULL,
                a1 = 0,
                hap = FALSE,
                iter = 10000,
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
    if (length(dim(y)) > 1)
        stop("dim(y) is greater than 1, y should be a vector")
    if (family == "binomial") {
        if (is.factor(y))
            y <- y != levels(y)[1L]
        if (any(y < 0 | y > 1))
            stop("y values must be 0 <= y <= 1")
    } else {
        if(!(typeof(y) %in% c("double", "integer")))
            stop("y must be numeric")
    }

    # check x and dim x/y
    n <- nrow(x)
    p <- ncol(x)
    y_len <- length(y)

    if (n != y_len)
        stop(paste("Length of y (", y_len, ") not equal to number of rows of x (", n, ")", sep = ""))
    if (class(x) != "matrix")
        x <- as.matrix(x)
    if (!(typeof(x) %in% c("double", "numeric", "integer")))
        stop("x contains non-numeric values")

    # check prior(model) if using sample
    if (method == "sample") {
        if (prior_model[[1]] < 0.0) {
            stop("alpha must be > 0")
        }
        if (prior_model[[2]] < 0.0) {
            stop("beta must be > 0")
        }
    }

    # check prior(coef)
    if (prior_coef[[1]] == "gprior") {
        if (length(prior_coef) != 4) {
            stop(paste("number of arguments (", length(prior_coef), ") for gprior not correct"))
        }
        if (prior_coef[[2]] == "default") {
            prior_coef[[2]] <- max(n, p^2)
        } else if (prior_coef[[2]] < 0) {
            stop("g must be nonnegative")
        }
        if (prior_coef[[3]] == "default") {
            prior_coef[[3]] <- 0.01
        } else if (prior_coef[[3]] < 0) {
            stop("alpha must be nonnegative")
        }
        if (prior_coef[[4]] == "default") {
            prior_coef[[4]] <- 0.01
        } else if (prior_coef[[4]] < 0) {
            stop("beta must be nonnegative")
        }
    } else if (prior_coef[[1]] != "none") {
        stop("prior_coef argument not recognized")
    }

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

    # check forced + intercept
    if (!is.null(forced)) {
        n_forced <- nrow(forced)
        p_forced <- ncol(forced)
        if (n_forced != y_len)
            stop(paste("Length of y (", y_len, ") not equal to number of rows of forced (", n_forced, ")", sep = ""))
        if (class(forced) != "matrix")
            forced <- as.matrix(forced)

        if (!(typeof(forced) %in% c("double", "integer")))
            stop("forced contains non-numeric values")
        if (intercept) {
            forced <- cbind(1, forced)
        }
    } else {
        p_forced <- 0
        if (intercept) {
            forced <- matrix(1, nrow = n, ncol = 1)
        }
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
        if (!(typeof(prior_cov) %in% c("double", "integer")))
            stop("prior_cov contains non-numeric values")
    } else {
        inform <- FALSE
        p_cov <- 0
    }

    # fit models by enumerating all combinations or M-H
    fit <- switch(method,
                  enumerate = bvs_enumerate(x, y, n, p, intercept, family, prior_coef,
                                            rare, hap, region_ind, forced, p_forced,
                                            inform, prior_cov, p_cov, a1, which_ind),

                  sample = bvs_sample(x, y, n, p, intercept, family, prior_model, prior_coef, rare,
                                      hap, region_ind, num_regions, forced, p_forced,
                                      inform, prior_cov, p_cov, which_ind, iter)
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
                           prior_model = prior_model,
                           prior_coef = prior_coef,
                           nobs = n,
                           nvars = p,
                           intercept = intercept,
                           varnames = varnames,
                           rare = rare,
                           regions = regions,
                           num_regions = num_regions,
                           inform = inform,
                           nvars_prior = p_cov,
                           hap = hap,
                           nullLogLike = fit$nullLogLike)
    fit$nullLogLike <- NULL
    class(fit) <- "bvs"
    return(fit)
}

