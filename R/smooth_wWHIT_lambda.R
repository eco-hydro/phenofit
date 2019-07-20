# Functions for working with the V-curve

#' @keywords internal
#' @importFrom spam diag.spam
v_point = function(y, w = 0 * y + 1, lambda = 100, d = 2) {
    # Compute the value of the normalized V-curve for one value of lambda
    # Prepare for smoothing
    n = length(y)
    E = diag.spam(n)
    D = diff(E, diff = d)
    P = t(D) %*% D

    # Smooth for  log-lambdas to the left and to the right
    z     = whit2(y, lambda, w)
    pz    = P %*% z #D' * D * z
    zgrad = lambda * log(10) * whit2(- pz/ w, lambda, w) #whit2(- pz * lambda, lambda, 1)
    # zgrad1 = whit2(-lambda * pz, lambda, w)

    fit   = sum(w * (y - z) ^ 2)
    dlfit = 2 * sum(-zgrad * w * (y - z)) / fit
    pen   = sum(z * pz)
    dlpen = 2 * sum(pz * zgrad) / pen

    # Take distance
    v = sqrt(dlfit ^ 2 + dlpen ^ 2)
    return(v)
}

# sometimes not converge
# v_opt = function(y, w = 0 * y + 1, d = 2, lambdas = c(0, 4), tol = 0.01) {
#     # Locate the optimal value of log10(lambda) with optimizer
#     # Specify bounds of search range for log10(lambda) in paramter 'lambdas'

#     v_fun = function(lla, y, w, d) v_point(y, w, 10 ^ lla, d)
#     op = optimize(v_fun, lambdas, y, w, d, tol = tol)
#     return(op$minimum)
# }

#' v_curve
#'
#' V-curve is used to optimize Whittaker parameter lambda.
#' Update 20180605 add weights updating to whittaker lambda selecting
#'
#' @inheritParams season
#' @inheritParams wWHIT
#' @param lg_lambdas `lg` lambda vectors of Whittaker parameter.
#' @param d Difference order.
#' @param IsPlot Boolean. Whether to plot figure?
#' 
#' @keywords internal
#' @examples
#' library(phenofit)
#' data("MOD13A1")
#' 
#' dt <- tidy_MOD13.gee(MOD13A1$dt)
#' st <- MOD13A1$st
#' 
#' sitename <- dt$site[1]
#' d     <- dt[site == sitename, ] # get the first site data
#' sp    <- st[site == sitename, ] # station point
#' # global parameter
#' IsPlot = TRUE
#' nptperyear = 23
#' 
#' dnew     <- add_HeadTail(d, nptperyear = nptperyear) # add one year in head and tail
#' INPUT    <- check_input(dnew$t, dnew$y, dnew$w, nptperyear, 
#'                         maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
#' # INPUT$y0 <- dnew$y   # raw time-series, for visualization
#' 
#' lg_lambdas <- seq(0, 3, 0.1)
#' r <- v_curve(INPUT, lg_lambdas, d = 2, IsPlot = TRUE)
#' @export
v_curve = function(INPUT, lg_lambdas, d = 2, IsPlot = FALSE,
    wFUN = wTSM, iters=2)
{
    # Compute the V-cure
    y <- INPUT$y
    w <- INPUT$w
    # w <- w*0 + 1
    nptperyear <- INPUT$nptperyear

    if (length(unique(y)) == 0) return(NULL)

    param <- c(INPUT, nptperyear = nptperyear, wFUN = wFUN, iters=iters,
        second = FALSE, lambda=NA)

    fits = pens = NULL
    for (lla in lg_lambdas) {
        # param$lambda <- 10^lla
        # z    <- do.call(wWHIT, param)$zs %>% dplyr::last()

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
    lambda  = 10 ^ lamids[k]

    # param$lambdas <- lambda
    # fit <- do.call(wWHIT, param)
    # d_sm <- fit %$% c(ws, zs) %>% as.data.table() %>% cbind(t = INPUT$t, .)

    z    <- whit2(y, lambda, w)
    d_sm <- data.table(t = INPUT$t, z)

    # result of v_curve
    vc <- list(lambda = lambda, vmin = v[k], 
        fit = d_sm, optim = data.table(lg_lambda = lamids, v = v)) 
    
    cal_COEF <- function(y){
        y <- y[!is.na(y)]
        list(mean = mean(y),
            sd = sd(y),
            kurtosis = kurtosis(y, type = 2),
            skewness = skewness(y, type = 2))
    }

    vc$coef_all <- cal_COEF(INPUT$y)  #%>% as.list()

    if (IsPlot) {
        par(mfrow = c(2, 1), mar = c(2.5, 2.5, 1, 0.2),
            mgp = c(1.3, 0.6, 0), oma = c(0, 0, 0.5, 0))

        ylim = c(0, max(v))
        plot(lamids, v, type = 'l', col = 'blue', ylim = ylim,
           xlab = 'log10(lambda)')
        points(lamids, v, pch = 16, cex = 0.5, col = 'blue' )
        abline(h = 0, lty = 2, col = 'gray')
        abline(v = lamids[k], lty = 2, col = 'gray', lwd = 2)
        title(sprintf("v-curve, lambda = %5.2f", lambda))
        grid()

        plot_input(INPUT, wmin = 0.2)
        colors <- c("blue", "red")

        lines(vc$fit$t, last(vc$fit), col = "blue", lwd = 1.2)
        # lines(vc$fit$t, vc$fit$ziter2, col = "red" , lwd = 1.2)
    }
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
