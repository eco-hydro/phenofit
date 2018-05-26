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

bisquare <- function(re, w, ...){
    re_abs <- abs(re)
    sc <- 6 * median(re_abs, na.rm = T)

    w_new <- (1 - (re/sc)^2)^2

    if (!missing(w)){
        I <- which(re > 0)
        if (length(I) > 0){
            w_new[I] <- w[I]#keep the original weights of positive bias points
        }
    }
    w_new[re >= sc] <- 0
    # constrain growing VI: enlarge the positive bias values, reduce the weights
    #   of negative bias values, as TIMESAT
    # diff = 2 * (yfit(i) - y(i)) / yfitstd;
    # w <- wfact * wfit(i) * exp( - ydiff ^ 2);
    return(w_new)
}
