# computes log of Beta-Binomial Prior
logBetaBinomial <- function(p, pgamma) {
    log(p) - lgamma(p + p + 1) + lgamma(1 + pgamma) + lgamma(p + p - pgamma)
}

# creates function to compute prob Beta-Binomial Prior given defined parameters
create_model_prior <- function(p, alpha, beta) {
    if (alpha == "default") {
        a <- 1
    } else {
        a <- alpha
    }
    if (beta == "default") {
        b <- p
    } else {
        b <- beta
    }
    if (alpha == "default" && beta == "default") {
        logBetaBinomial <- function(pgamma) {
            log(p) - lgamma(p + p + 1) + lgamma(1 + pgamma) + lgamma(p + p - pgamma)
        }
    } else {
        logBetaBinomial <- function(pgamma) {
            lgamma(p + 1) - lgamma(pgamma + 1) - lgamma(p - pgamma + 1) + lgamma(pgamma + a) + lgamma(p - pgamma + b) - lgamma(p + a + b) + lgamma(a + b) - lgamma(a) - lgamma(b)
        }
    }
    return(logBetaBinomial)
}
