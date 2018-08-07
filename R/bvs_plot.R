#' Plots for Top Variant and Region Inclusions

#' @description This function allows the user to create image plots of the top variants and top Regions
#' (any user specified set of variants such as pathways or genes) included in the top models.
#' Variants and Regions are ordered based on marginal BF and regional BF which are plotted on the right axis.
#' The width of the inclusion blocks are proportional to the posterior model probability that the variant or region is included in.
#'
#' @param x an object of class 'summary.bvs'
#' @param type specifies whether to plot the top variants ("s") or the top regions ("r")
#' @param num_models the number of top models to place on the x-axis
#' @param num_snps if type = "s", the number of the top variants to place on the y-axis
#' @param num_regions if type = "r", the number ofthe top regions to place on the y-axis
#' @param plot_coef only used for rare variant analysis when rare = TRUE and there are not
#' multiple regions. If plot_coef = TRUE, the log(OR) of the risk index for the top models is plotted on the x-axis
#' @param true_coef (optional) vector of the true odds ratios of each of the variants to plot on the
#' y-axis (i.e. if results are from a simulation)
#' @param regions (optional) string vector with the region name for each of the variants. By default, region
#' names are used from the 'summary.bvs' x. Using this argument will overwrite the names in the "summary.bvs" x.
#' @param prop_cases (optional) \eqn{p x 2} matrix giving the number of cases that have the
#' variant in column 1 and the number of controls with the variant in column 2.
#' If specified, these counts will be reported on the right axis under each variants marginal BF.
#' @param main optional string variable giving the title of the plot
#' @param ... additional arguments as required by plot S3 x

#' @importFrom graphics abline axis image par

#' @export
plot.summary.bvs <- function(x,
                             type = c("s", "r"),
                             num_models = 100,
                             num_snps = 20,
                             num_regions = 20,
                             plot_coef = FALSE,
                             true_coef = NULL,
                             regions = NULL,
                             prop_cases=NULL,
                             main = NULL, ...) {

    # check type
    type <- match.arg(type)
    if (type == "r" && is.null(x$region_level$active_region)) {
        stop("Error: type = 'r', but active_region = NULL in bvs.summary x")
    }

    active <- x$model_level$active_mat
    post_prob <- x$model_level$post_prob
    null_ind <- which(rowSums(active) == 0)
    null_post <- post_prob[null_ind]
    model_id <- x$model_level$model_id[-null_ind]
    active <- active[-null_ind, , drop = FALSE]
    post_prob <- post_prob[-null_ind]
    num_snps <- min(ncol(active), num_snps)
    num_models <- min(nrow(active) - 1, num_models)
    model_order <- order(post_prob, decreasing = TRUE)
    model_id <- model_id[model_order][1:num_models]
    post_prob <- post_prob[model_order]
    active <- active[model_order, , drop = FALSE]

    if (!is.null(x$region_level$active_region)) {
        active_region <- x$region_level$active_region[-null_ind, ][model_order, , drop = FALSE]
        regionnames <- colnames(active_region)
        if (!is.null(regions)) {
            warning("Note: Overwriting regions in bvs.summary x with values provided for regions = argument.")
        } else {
            regions <- x$model_info$regions
        }
    }

    if (plot_coef) {
        if (ncol(coef) != 1) {
            warning("Note: Coercing plot_coef to FALSE because ncol(coef) > 1.
                    To use plot_coef, check that your model has rare = TRUE and does not have multiple regions.")
            plot_coef <- FALSE
        } else {
            coef <- drop(x$model_level$coef[-null_ind, ][model_order, ])
        }
    }

    # create title, column / row labels, and matrix with color values
    if (type == "s") {
        nvar <- num_snps
        bf <- x$marg_bf
        order_top <- order(bf, decreasing = TRUE)[1:nvar]
        bf <- bf[order_top]
        rownms <- paste(colnames(active)[order_top], regions[order_top], sep = "\n")
        color_matrix <- active[1:num_models, order_top, drop = FALSE] + 2

        if (plot_coef) {
            prob_labels <- paste("Marg BF:", round(bf, 2), "\nTrue OR:", round(true_coef[order_top], 2), sep = "")
        } else if (!is.null(prop_cases)) {
            if (plot_coef) {
                warning("prop_cases cannot be used when plot_coef == TRUE, only plotting log(OR)")
            }
            prob_labels <- paste("Marg BF:", round(bf, 2), "\nCases: ", prop_cases[order_top, 1], " Controls: ", prop_cases[order_top, 2],sep = "")
        } else {
            prob_labels <- paste("Marg BF:", round(bf, 2))
        }

        if (is.null(main)) {
            main <- paste("SNP Inclusions of Top Models \nGlobal BF =", round(x$global_bf, 1))
        }
    } else {
        nvar <- min(length(regionnames), num_regions)
        bf_region <- x$region_level$region_bf
        order_top <- order(bf_region, decreasing = TRUE)[1:nvar]
        bf_region <- bf_region[order_top]
        rownms <- colnames(active_region)[order_top]
        active_region[active_region > 1] <- 1
        color_matrix <- active_region[1:num_models, order_top] + 2

        if (!is.null(true_coef)) {
            prob_labels <- paste("Region BF:", round(bf_region, 2), " \nTrue OR:", round(true_coef[order_top], 2), sep = "")
        } else {
            prob_labels <- paste("Region BF:", round(bf_region, 2))
        }

        if (is.null(main)) {
            main <- paste("Region Inclusions of Top Models \nGlobal BF =", round(x$global_bf, 1))
        }
    }

    # create plot
    keep.mar <- par(mar = c(5, 6, 4, 2) + 0.1)
    par(las = 1, mar = c(8, 12, 5, 12), ps = 10, font = 2)
    prob_axis <- post_prob[1:num_models] / sum(post_prob[1:num_models])
    prob_axis_cum <- cumsum(prob_axis)
    clr <- c("#FFFFFF", "#A020F0", "#0000CD")
    if (plot_coef) {
        xlab = "log(OR) (Models Ranked by Post. Prob)"
    } else {
        xlab = "Model ID (Models Ranked by Post. Prob)"
    }
    image(x = c(0, prob_axis_cum),
          y = 1:nvar,
          z = color_matrix,
          col = clr,
          xlab = xlab,
          ylab = "",
          xaxt = "n",
          yaxt ="n",
          xlim = c(0, 1),
          main = main)
    abline(v = prob_axis_cum[prob_axis > 0.01], col = "white")
    xat <- (prob_axis_cum + c(0, prob_axis_cum[-num_models])) / 2
    if (plot_coef) {
        beta.labels <- round(coef, 1)
        if (num_models > 5) {
            beta.labels[6:num_models] <- NA
        }
        axis(1, at = xat,labels = beta.labels)
    } else {
        axis(1, at = xat, labels = model_id)
    }
    axis(2, at = 1:nvar, labels = rownms)
    axis(4, at = 1:nvar, labels = prob_labels)
    par(mar = keep.mar)
}
