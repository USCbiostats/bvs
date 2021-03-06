% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bvs.R
\name{bvs}
\alias{bvs}
\title{Bayesian Variable Selection}
\usage{
bvs(y, x, forced = NULL, intercept = TRUE, family = c("gaussian",
  "binomial"), method = c("sample", "enumerate"), prior_model = list(alpha =
  1, beta = p), prior_coef = list("none"), rare = FALSE, regions = NULL,
  prior_cov = NULL, a1 = 0, hap = FALSE, iter = 10000, maxk = 3,
  parallel = FALSE, control = list())
}
\arguments{
\item{y}{outcome variable}

\item{x}{predictor design matrix}

\item{forced}{(optional) \eqn{n x c} matrix of \emph{c} confounding variables that one wishes to adjust the
analysis for and that will be forced into every model.}

\item{intercept}{indicates whether models should include an intercept. Default = TRUE.}

\item{family}{specifies error distribution for outcome variable, options include
\itemize{
   \item "gaussian" (default)
   \item "binomial"
}}

\item{method}{specifies how to search the model space, options include
\itemize{
   \item "sample" (default) : performs basic Metropolis-Hastings algorithm to sample models from model space. For informative
   marginal inclusion probabilities, the algorithm also performs basic MCMC algorithm to sample the effects of predictor-level
   covariates (alpha).
   \item "enumerate": computes and summarizes all possible models in model space. Not recommended for problems where \eqn{p > 20}.
}}

\item{prior_model}{specifies parameters for beta-binomial prior on model size. To specify, pass a list with the following elements
\itemize{
   \item \eqn{alpha} = numeric value for first shape parameter (default = 1)
   \item \eqn{beta} = numeric value for second shape parameter (default = p)
}
Example: \code{list(alpha = 1, beta = 2)}}

\item{prior_coef}{specifies prior for regression coefficients (only for use when family = "gaussian"), options include
\itemize{
   \item "none" (default)
   \item "gprior": To specify the parameters for the g-prior, pass a list object with the following elements
   \itemize{
   \item "gprior"
   \item \eqn{g} = numeric value for sparsity parameter, recommended default = \eqn{max(n, p^2)}
   \item \eqn{alpha} = numeric value that specifies variance in g-prior, recommended default = 0.01
   \item \eqn{beta} = numeric value that specifies variance in g-prior, recommended default = 0.01
   }
   To use the default values, pass the string "default" for any parameter.

   Example: \code{list("gprior", "default", 0.02, 0.02)}
}}

\item{rare}{if rare = TRUE, corresponds to the Bayesian Risk Index (BRI) algorithm of Quintana and Conti (2011) that constructs
a risk index based on the multiple rare variants within each model. The marginal likelihood of each model is then calculated
based on the corresponding risk index.}

\item{regions}{(optional) \eqn{p x 1} character or factor vector that identifies a user-defined region for each variant. If
rare = TRUE, then multiple region-specific risk indices are computed for each model.}

\item{prior_cov}{(optional) if method = "sample", a \eqn{p x q} matrix of \emph{q} predictor-level covariates that the user
wishes to incorporate into the estimation of the marginal inclusion probabilities using the iBMU algorithm.}

\item{a1}{(optional) if method = "enumerate", a \eqn{q x 1} vector of specified effects of each predictor-level covariate.}

\item{hap}{(not yet implemented) if hap = TRUE, esimtate a set of haplotypes from the multiple variants within each moel and the marginal likelihood
of each model is calculated based on the set of haplotypes.}

\item{iter}{if method = "sample", the number of iterations to run the algorithm. Default = 1000.}

\item{maxk}{if method = "enumerate", the maximum model size (k) to consider when enumerating all possible models. Default = 3.}

\item{control}{specifies 'bvs' control object.}
}
\description{
Bayesian Variable Selection
}
