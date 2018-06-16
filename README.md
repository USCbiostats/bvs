bvs: Bayesian Variable Selection
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
Introduction
============

The functions in this package focus on analyzing case-control association studies involving a group of genetic variants. In particular, we are intersted in modeling the outcome variable as a function of multivariate genetic profile using Bayesian model uncertainty and variable selection techniques. The package incorporates functions to analyze data sets involving common variants as well as externsions to model rare variants via the Bayesian Risk Index (Quintana and Conti, 2011), as well as haplotypes. Finally, allows the incorporation of biological information to inform the marginal inclusion probabilities via the iBMU (Quintana and Conti (submitted)).

Setup
=====

1.  For Windows users, install [RTools](https://cran.r-project.org/bin/windows/Rtools/) (not an R package)
2.  Install the R package [devtools](https://github.com/hadley/devtools)
3.  Install bvs package with the install\_github function
4.  Load the package

``` r
library(devtools)

# Master branch
install_github("USCbiostats/bvs")
```

Example: Rare Variant Analysis with the Bayesian Risk Index (BRI)
=================================================================

As a first example, below is code to run a rare variant analysis using the "sample" method on a binary outcome.

``` r
# load the rare data set, first column is outcome variable
load(RareData)

# run 2000 iterations of the MH algorithm (rare = TRUE will construct the BRI) 
rare_analysis <- bvs(y = RareData$case, 
                     x = as.matrix(RareData[, -1]),
                     family = "binomial",
                     method = "sample", 
                     rare = TRUE,
                     iter = 2000)

# summarize results of analysis
rare_results <- summary(rare_analysis)

# plot of top 25 models and top 10 variants
plot(rare_results, num_models = 25, num_snps = 10)
```
