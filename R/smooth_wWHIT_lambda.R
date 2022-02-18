# ' @param d integer, the order of difference of Whittaker smoother.

# Functions for working with the V-curve
#' lambda_vcurve
#'
#' @inheritParams check_input
#' @param lg_lambdas numeric vector, log10(lambda) candidates. The optimal `lambda`
#' will be optimized from `lg_lambda`.
#' @param plot logical. If `TRUE`, V-curve will be plotted.
#' @param plot logical. If `TRUE`, the optimized `lambda` will be printed on the 
#' console.
#' 
#' @param ... ignored.
#' 
#' @keywords internal
#' @export
lambda_vcurve <- function(y, w,
    # d = 2,
    lg_lambdas = seq(0.1, 5, 0.1),
    plot = FALSE, 
    verbose = FALSE, ...)
{
    d = 2
    fits = pens = NULL
    for (lla in lg_lambdas) {
        # param$lambda <- 10^lla
        # z    <- do.call(smooth_wWHIT, param)$zs %>% last()
        z    = whit2(y, 10 ^ lla, w)
        fit  = log(sum(w * (y - z) ^ 2))
        pen  = log(sum(diff(z, diff = d) ^2))
        fits = c(fits, fit)
        pens = c(pens, pen)
    }

    # Construct V-curve
    dfits   = diff(fits)
    dpens   = diff(pens)
    llastep = lg_lambdas[2] - lg_lambdas[1]
    v       = sqrt(dfits ^ 2 + dpens ^ 2) / (log(10) * llastep)

    nla     = length(lg_lambdas)
    lamids  = (lg_lambdas[-1] + lg_lambdas[-nla]) / 2
    k       = which.min(v)
    opt_lambda  = 10 ^ lamids[k]

    if (plot) {
        par(mfrow = c(2, 1), mar = c(2.5, 2.5, 1, 0.2),
            mgp = c(1.3, 0.6, 0), oma = c(0, 0, 0.5, 0))

        ylim = c(0, max(v))
        plot(lamids, v, type = 'l', col = 'blue', ylim = ylim,
           xlab = 'log10(lambda)')
        points(lamids, v, pch = 16, cex = 0.5, col = 'blue' )
        abline(h = 0, lty = 2, col = 'gray')
        abline(v = lamids[k], lty = 2, col = 'gray', lwd = 2)
        title(sprintf("v-curve, lambda = %5.2f", opt_lambda))
        grid()
        # plot_input(INPUT, wmin = 0.2)
        plot(y, type = "l")
        z <- whit2(y, opt_lambda, w)
        lines(z, col = "red", lwd = 0.5)
        # colors <- c("blue", "red")
        # lines(vc$fit$t, last(vc$fit), col = "blue", lwd = 1.2)
        # lines(vc$fit$t, vc$fit$ziter2, col = "red" , lwd = 1.2)
    }
    if (verbosee) {
        fprintf("The optimized `lambda` by V-curve theory is %.2f", opt_lambda)
    }
    list(lambda = opt_lambda, vcurve = data.table(lg_lambda = lamids, v = v))
    # opt_lambda
}

#' @rdname lambda_vcurve
#' @export
lambda_vcurve_jl <- function(y, w, lg_lambda_min = 0.1, lg_lambda_max = 3) {
    JuliaCall::julia_call("phenofit.lambda_vcurve", y, w,
        lg_lambda_min = lg_lambda_min,
        lg_lambda_max = lg_lambda_max)
}

#' @rdname lambda_vcurve
#' @export
lambda_cv_jl <- function(y, w, d = 2, lg_lambda_min = 0.1, lg_lambda_max = 3) {
    JuliaCall::julia_call("phenofit.lambda_cv", y, w,
        lg_lambda_min = lg_lambda_min,
        lg_lambda_max = lg_lambda_max)
}

#' V-curve theory to optimize Whittaker parameter `lambda`.
#'
#' @description
#' V-curve is used to optimize Whittaker parameter `lambda`.
#' This function is not for users!!!
#'
#' Update 20180605 add weights updating to whittaker lambda selecting
#'
#' @inheritParams season
#' @inheritParams smooth_wWHIT
#' @inheritParams lambda_vcurve
#' @param ... ignored.
#' 
#' @keywords internal
#' @seealso lambda_vcurve
#' @examples
#' data("CA_NS6"); d = CA_NS6
#' INPUT <- check_input(d$t, d$y, d$w, nptperyear = 23,
#'     maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
#'
#' r <- v_curve(INPUT, lg_lambdas = seq(0, 3, 0.1), plot = TRUE)
#' @export
v_curve = function(
    INPUT,
    lg_lambdas = seq(0, 5, by = 0.005),
    # d = 2,
    # wFUN = wTSM, iters=2,
    plot = FALSE,
    ...)
{
    cal_coef <- function(y){
        y <- y[!is.na(y)]
        list(mean = mean(y),
            sd = sd(y),
            kurtosis = kurtosis(y, type = 2),
            skewness = skewness(y, type = 2))
    }

    # Compute the V-cure
    y <- INPUT$y
    w <- INPUT$w
    # w <- w*0 + 1
    nptperyear <- INPUT$nptperyear
    if (length(unique(y)) == 0) return(NULL)

    # , wFUN = wFUN, iters=iters
    param <- c(INPUT, nptperyear = nptperyear,
        second = FALSE, lambda=NA)

    # param$lambdas <- lambda
    # fit <- do.call(smooth_wWHIT, param)
    # d_sm <- fit %$% c(ws, zs) %>% as.data.table() %>% cbind(t = INPUT$t, .)
    l_opt = lambda_vcurve(y, w, 
        # d = d,
        lg_lambdas = lg_lambdas, plot = plot) # optimal lambda

    z <- whit2(y, l_opt$lambda, w)
    d_sm <- data.table(t = INPUT$t, z) # smoothed
    # result of v_curve
    vc <- list(lambda = l_opt$lambda,
               vmin = l_opt$vcurve$v %>% min(),
        fit = d_sm, optim = l_opt$vcurve)

    vc$coef_all <- cal_coef(INPUT$y) # %>% as.list()
    vc
}

## All year togather
# Call:
# lm_fun(formula = lambda ~ mean + sd + cv + skewness, data = d,
#     na.action = na.exclude)

# Residuals:
#     Min      1Q  Median      3Q     Max
# -2.4662 -0.4267  0.1394  0.4144  2.5824

# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
# (Intercept)  0.831120   0.021897  37.956  < 2e-16 ***
# mean         1.599970   0.089914  17.794  < 2e-16 ***
# sd          -4.094027   0.168844 -24.247  < 2e-16 ***
# cv          -0.035160   0.008459  -4.157 3.25e-05 ***
# skewness    -0.063533   0.007966  -7.976 1.62e-15 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.5851 on 15135 degrees of freedom
# Multiple R-squared:  0.1572,    Adjusted R-squared:  0.1569
# F-statistic: 705.5 on 4 and 15135 DF,  p-value: < 2.2e-16

#          term    estimate   std.error  statistic       p.value
# 1 (Intercept)  0.83112003 0.021897138  37.955647 3.306152e-301
# 2        mean  1.59997023 0.089914112  17.794428  4.039120e-70
# 3          sd -4.09402663 0.168843860 -24.247412 1.876008e-127
# 4          cv -0.03516045 0.008458787  -4.156677  3.246863e-05
# 5    skewness -0.06353256 0.007965580  -7.975886  1.620571e-15

## Three year fitting result
# Call:
# lm_fun(formula = lambda ~ (mean + sd + cv + skewness + kurtosis),
#     data = d, na.action = na.exclude)

# Residuals:
#     Min      1Q  Median      3Q     Max
# -2.6601 -0.4564  0.0490  0.4418  2.7551

# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
# (Intercept)  0.817783   0.010569  77.379  < 2e-16 ***
# mean         1.803588   0.042016  42.926  < 2e-16 ***
# sd          -4.263469   0.081937 -52.033  < 2e-16 ***
# cv          -0.038240   0.004041  -9.462  < 2e-16 ***
# skewness    -0.066914   0.003762 -17.785  < 2e-16 ***
# kurtosis     0.011289   0.001506   7.496 6.62e-14 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.6639 on 90639 degrees of freedom
# Multiple R-squared:  0.1541,    Adjusted R-squared:  0.1541
# F-statistic:  3303 on 5 and 90639 DF,  p-value: < 2.2e-16

#          term    estimate   std.error  statistic      p.value
# 1 (Intercept)  0.81778299 0.010568590  77.378628 0.000000e+00
# 2        mean  1.80358830 0.042016401  42.925816 0.000000e+00
# 3          sd -4.26346883 0.081937136 -52.033413 0.000000e+00
# 4          cv -0.03823967 0.004041253  -9.462331 3.080333e-21
# 5    skewness -0.06691403 0.003762368 -17.785083 1.216689e-70
# 6    kurtosis  0.01128888 0.001505916   7.496350 6.621350e-14
