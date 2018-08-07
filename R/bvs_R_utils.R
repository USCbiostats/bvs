# computes log of Beta-Binomial Prior
logBetaBinomial <- function(pgamma, p, alpha, beta) {
    lbeta(pgamma + alpha, p - pgamma + beta) - lbeta(alpha, beta) + lchoose(p, pgamma)
}
