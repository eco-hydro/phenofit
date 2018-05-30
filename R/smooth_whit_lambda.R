# Functions for working with the V-curve

# require(spam)
# require(ptw)

whit <- function (y, lambda, w = rep(1, ny), iters = 2) {
    ny <- length(y)
    if (all(is.na(y))){
        return(list(z = y, w = w))
    }

    for (i in 1:iters){
        z  <- whit2(y, lambda, w)
        re <- z - y
        if (i < iters){
            w <- wBisquare(re, w) #update weights
        }
    }
    return(list(z = z, w = w))
}

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

v_curve = function(y, w = 0 * y + 1, llas,  d = 2, show = F) {
  # Compute the V-cure
  fits = pens = NULL
  for (lla in llas) {
    z = whit2(y, 10 ^ lla, w)
  	fit = log(sum(w * (y - z) ^ 2))
	  pen = log(sum(diff(z, diff = d) ^2))
	  fits = c(fits, fit)
	  pens = c(pens, pen)
  }

  # Construct V-curve
  dfits   = diff(fits)
  dpens   = diff(pens)
  llastep = llas[2] - llas[1]
  v      = sqrt(dfits ^ 2 + dpens ^ 2) / (log(10) * llastep)
  nla    = length(llas)
  lamids = (llas[-1] + llas[-nla]) / 2
  k      = which.min(v)
  lambda = 10 ^ lamids[k]
  z      = whit2(y, lambda, w)

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
    kurtosis <- kurtosis(y, type = 2)
    skewness <- skewness(y, type = 2)
    # lambda   <- 0.555484 + 1.671514*mean - 3.434064*sd - 0.052609*skewness + 0.009057*kurtosis
    lambda   <- 0.555465 + 1.501239*mean - 3.204295*sd - 0.031902*skewness # Just three year result
    return(lambda)
}

## All year togather
# Call:
# lm(formula = lambda ~ kurtosis + mean + sd + skewness, data = stat)

# Residuals:
#      Min       1Q   Median       3Q      Max 
# -2.33795 -0.33385  0.04024  0.37228  2.08361 

# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.555484   0.016632  33.399  < 2e-16 ***
# kurtosis     0.009057   0.002746   3.299 0.000973 ***
# mean         1.671514   0.038983  42.878  < 2e-16 ***
# sd          -3.434064   0.118632 -28.947  < 2e-16 ***
# skewness    -0.052609   0.008332  -6.314 2.79e-10 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.5506 on 15074 degrees of freedom
# Multiple R-squared:  0.2118,  Adjusted R-squared:  0.2116 
# F-statistic:  1013 on 4 and 15074 DF,  p-value: < 2.2e-16

## Three year fitting result
# Call:
# lm(formula = lambda ~ mean + sd + skewness, data = stat, na.action = na.exclude)

# Residuals:
#      Min       1Q   Median       3Q      Max 
# -2.77002 -0.40751  0.01012  0.42447  2.02898 

# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.555465   0.007000  79.357   <2e-16 ***
# mean         1.501239   0.017132  87.628   <2e-16 ***
# sd          -3.204295   0.044921 -71.332   <2e-16 ***
# skewness    -0.031902   0.003406  -9.367   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.6188 on 90265 degrees of freedom
#   (4 observations deleted due to missingness)
# Multiple R-squared:  0.1441,  Adjusted R-squared:  0.144 
# F-statistic:  5064 on 3 and 90265 DF,  p-value: < 2.2e-16

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