hap_bvs <- function(g, min_hap_freq = 0.02) {

    # convert G to allele matrix for haplo.stats
    g_alleles <- as.data.frame(codeG_AsAlleles(g))
    names(g_alleles) <- rep(names(g), 1, each=2)

    # estimate haplotypes
    hap_obj <- haplo.em(g_alleles, locus.label = names(g), miss.val = NA)

    ##Create Posterior Haplotype Matrix and Return Estimated Haplotypes
    x_post <- createPostHapMatrix(hap_obj, min_hap_freq = min_hap_freq)
    hap_index <- (1:length(hap_obj$hap.prob))[hap_obj$hap.prob >= min_hap_freq]
    haplotypes <- hap_obj$haplotype[hap_index,]
    x_haps <- as.data.frame(x_post[, 2:ncol(x_post)])
    names(x_haps) <- paste("Hap", hap_index[2:length(hap_index)], sep = "")
    return(x_haps)
}


codeG_AsAlleles <- function(x) {
    v <- {}
    for(m in 1:ncol(x)) {
        v <- cbind(v, codeAlleles(x[, m]))
    }
    v
}

codeAlleles <- function(x) {
    a1 <- ifelse(is.na(x), NA, ifelse(x <= 1, 0, 1))
    a2 <- ifelse(is.na(x), NA, ifelse(x >= 1, 1, 0))
    alleles <- cbind(a1, a2)
    alleles
}

createPostHapMatrix <- function(hap_obj, min_hap_freq = 0.025) {
    hap_prob <- hap_obj$hap.prob
    haplotypes <- hap_obj$haplotype

    hap1 <- hap_obj$hap1code
    hap2 <- hap_obj$hap2code
    indx <- hap_obj$indx.subj
    post <- hap_obj$post
    nreps <- as.vector(hap_obj$nreps)
    uhap  <- sort(unique(c(hap1,hap2)))

    which_haplo <- hap_obj$hap.prob >= min_hap_freq
    uhap <- uhap[which_haplo]
    x <- outer(hap1, uhap, "==") + outer(hap2, uhap, "==")

    n_subj <- length(unique(hap_obj$indx.subj))
    n_x <- ncol(x)
    x_post <- matrix(rep(NA, n_subj * n_x), ncol = n_x)

    for(j in 1:n_x){
        x_post[, j] <- tapply(x[, j] * post, indx, sum)
    }
    x_post
}
