# computes log of Beta-Binomial Prior
logBetaBinomial <- function(p, pgamma) {
    log(p) - lgamma(p + p + 1) + lgamma(1 + pgamma) + lgamma(p + p - pgamma)
}
