#' Weighted HANTS SMOOTH
#'
#' Weighted HANTS smoother
#' 
#' @template INPUT
#' @param nf number of frequencies to be considered above the zero frequency
#' @param ylu [low, high] of time-series y (curve fitting values are constrained
#' in the range of \code{ylu}.
#' @param nptperyear number of images per year
#' @param periodlen length of the base period, measured in virtual samples
#'           (days, dekads, months, etc.). nptperyear in timesat.
#' @template WEIGHT
#' @param ... Additional parameters are passed to \code{wFUN}.
#' 
#' @author
#' Wout Verhoef, NLR, Remote Sensing Dept. June 1998
#' Mohammad Abouali (2011), Converted to MATLAB
#' Dongdong Kong (2018), introduced to R and modified into weighted model.
#' @return
#' \itemize{
#'    \item \code{ws} weights of every iteration
#'    \item \code{zs} curve fittings of every iteration
#' }
#' @export
# Modified:
#   Apply suppression of high amplitudes for near-singular case by
#   adding a number delta to the diagonal elements of matrix A,
#   except element (1,1), because the average should not be affected
#
#   Output of reconstructed time series in array yr
#
#   Change call and input arguments to accommodate a base period length (nptperyear)
#   All frequencies from 1 (base period) until nf are included
# RESULT:
# amp returned array of amplitudes, first element is the average of
#         the curve
# phi returned array of phases, first element is zero
# yr array holding reconstructed time series
wHANTS <- function(y, t, w, nf = 3, ylu, periodlen = 365, nptperyear,
                   wFUN = wTSM, iters = 2, wmin = 0.1, ...){
    if (is.Date(t)) t <- as.numeric(t - t[1])

    n    <- length(y)
    ncol <- min(2*nf+1, n)

    mat <- array(0, dim = c(n, ncol) )
    amp <- numeric(nf+1)
    phi <- numeric(nf+1)
    yz  <- numeric(n)

    ang <- 2*pi*(0:(periodlen-1))/periodlen
    cs  <- cos(ang)
    sn  <- sin(ang)

    mat[, 1] =1.0
    #  f*2*pi*[0:nptperyear-1]/nptperyear mod replace it
    for (i in 1:nf){
        I = 1 + (i*(t - 1)) %% periodlen
        mat[, 2*i  ] = cs[I]
        mat[, 2*i+1] = sn[I]
    }

    fits <- list()
    ws   <- list()
    for (i in 1:iters){
        ws[[i]] <- w
        za = t(mat) %*% (w*y)

        A = t(mat*w) %*% mat # mat' * diag(w) * mat
        # % A = A + diag(ones(nr,1))*delta
        # % A(1,1) = A(1,1) - delta
        b <- solve(A, za) # coefficients
        z <- mat %*% b; z <- z[, 1]
        # w = wFUN(y, yr, w, 0.5, i, nptperyear) #%wfact = 0.5
        wnew <- wFUN(y, z, w, i, nptperyear, ...)
        # \code{check_fit} has constrained ylu
        # print(unique(wnew - w))
        # w <- wnew
        # w[z < ylu[1] | z > ylu[2] ] <- wmin

        z <- check_fit(z, ylu) # very necessary
        fits[[i]] <- z
    }

    amp[1]   = b[1]
    phi[1]   = 0.0
    i   <- seq(2, ncol, 2)
    ifr <- (i+2)/2
    ra  <- b[i]
    rb  <- b[i+1]
    amp[ifr] = sqrt(ra*ra+rb*rb)

    dg    <- 180.0/pi
    phase <- atan2(rb, ra)*dg
    phase[phase<0] = phase[phase<0] + 360
    phi[ifr] = phase

    fits %<>% set_names(paste0('ziter', 1:iters))
    ws   %<>% set_names(paste0('witer', 1:iters))
    list(ws = ws, zs = fits)
    # list(fit = fits, amp = amp, phi = phi) #quickly return
}
