# computes log of Beta-Binomial Prior for model size
computeModelPrior <- function(pgamma, p, alpha, beta) {
    lbeta(pgamma + alpha, p - pgamma + beta) - lbeta(alpha, beta)
}
