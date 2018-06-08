# Functions for working with the V-curve

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
v_opt = function(y, w = 0 * y + 1, d = 2, llas = c(0, 4), tol = 0.01) {
    # Locate the optimal value of log10(lambda) with optimizer
    # Specify bounds of search range for log10(lambda) in paramter 'llas'
    
    v_fun = function(lla, y, w, d) v_point(y, w, 10 ^ lla, d)
    op = optimize(v_fun, llas, y, w, d, tol = tol)
    return(op$minimum)
}

#' update 20180605 add weights updating to whittaker lambda selecting
#' @export
v_curve = function(INPUT, nptperyear, llas,  d = 2, show = F, 
    wFUN = wTSM, iters=2) {
    # Compute the V-cure
    y <- INPUT$y
    w <- INPUT$w
    
    param <- c(INPUT, nptperyear = nptperyear, wFUN = wFUN, iters=iters, 
        second = FALSE, lambdas=NA)
    # whitsmw2(y, w, ylu, nptperyear, wFUN = wTSM, iters=1, lambdas=1000,

    fits = pens = NULL
    for (lla in llas) {
        # param$lambdas <- 10^lla
        # z    <- do.call(whitsmw2, param) %>% dplyr::last()
        z    = whit2(y, 10 ^ lla, w)
        fit  = log(sum(w * (y - z) ^ 2))
        pen  = log(sum(diff(z, diff = d) ^2))
        fits = c(fits, fit)
        pens = c(pens, pen)
    }

    # Construct V-curve
    dfits   = diff(fits)
    dpens   = diff(pens)
    llastep = llas[2] - llas[1]
    v       = sqrt(dfits ^ 2 + dpens ^ 2) / (log(10) * llastep)

    nla     = length(llas)
    lamids  = (llas[-1] + llas[-nla]) / 2
    k       = which.min(v)
    lambda  = 10 ^ lamids[k]
    
    # param$lambdas <- lambda
    # z    <- do.call(whitsmw2, param) %>% dplyr::last()
    z       = whit2(y, lambda, w)

    if (show) {
        ylim = c(0, max(v))
        plot(lamids, v, type = 'l', col = 'blue', ylim = ylim,
           xlab = 'log10(lambda)')
        points(lamids, v, pch = 16, cex = 0.5, col = 'blue' )
        abline(h = 0, lty = 2, col = 'gray')
        abline(v = lamids[k], lty = 2, col = 'gray', lwd = 2)
        title('V-curve')
        grid()
    }
    return(list(z = z, llas = lamids, lambda = lambda, v = v, vmin = v[k]))
}

#' Initial lambda value of whittaker taker
#' @export
init_lambda  <- function(y){
    y        <- y[!is.na(y)] #rm NA values
    mean     <- mean(y)
    sd       <- sd(y)
    cv       <- sd/mean
    skewness <- skewness(y, type = 2)
    kurtosis <- kurtosis(y, type = 2)
    # lambda was transformed by log10
    # lambda   <- 0.555484 + 1.671514*mean - 3.434064*sd - 0.052609*skewness + 0.009057*kurtosis
    # lambda   <- 0.555465 + 1.501239*mean - 3.204295*sd - 0.031902*skewness # Just three year result
    
    lambda <- 0.831120 + 1.599970*mean - 4.094027*sd - 0.035160*cv - 0.063533*skewness # all year 
    # lambda <- 0.817783 + 1.803588*mean - 4.263469*sd - 0.038240*cv - 0.066914*skewness - 0.011289*kurtosis  #3y
    return(10^lambda)
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

#' kurtosis
#' 
#' Inherit from package `e1071`
#' 
#' @export
kurtosis <- function (x, na.rm = FALSE, type = 3) {
    if (any(ina <- is.na(x))) {
        if (na.rm) {
            x <- x[!ina]            
        } else {
            return(NA)   
        }
    }
    if (!(type %in% (1:3))) stop("Invalid 'type' argument.")
    n <- length(x)
    x <- x - mean(x)
    r <- n * sum(x^4)/(sum(x^2)^2)

    y <- if (type == 1) {
        r - 3
    } else if (type == 2) {
        if (n < 4) stop("Need at least 4 complete observations.")
        ((n + 1) * (r - 3) + 6) * (n - 1)/((n - 2) * (n - 3))
    } else{
        r * (1 - 1/n)^2 - 3
    }
    y
}

#' @export
skewness <- function (x, na.rm = FALSE, type = 3) {
    if (any(ina <- is.na(x))) {
        if (na.rm) {
            x <- x[!ina]            
        } else {
            return(NA)   
        }
    }
    if (!(type %in% (1:3))) stop("Invalid 'type' argument.")
    n <- length(x)
    x <- x - mean(x)
    y <- sqrt(n) * sum(x^3)/(sum(x^2)^(3/2))
    
    if (type == 2) {
        if (n < 3) stop("Need at least 3 complete observations.")
        y <- y * sqrt(n * (n - 1))/(n - 2)
    } else if (type == 3){
        y <- y * ((1 - 1/n))^(3/2)
    }
    y
}