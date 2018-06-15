# Beta-Binomial Prior
BetaBinomial <- function(p, pgamma) {
    log(p) - lgamma(p + p + 1) + lgamma(1 + pgamma) + lgamma(p + p - pgamma)
}
