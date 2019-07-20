#' init_param
#'
#' Initialize parameters of double logistic function
#'
#' @inheritParams FitDL.Zhang
#'
#' @keywords internal
#' @examples
#' library(phenofit)
#' # simulate vegetation time-series
#' fFUN = doubleLog.Beck
#' par  = c(
#'     mn  = 0.1,
#'     mx  = 0.7,
#'     sos = 50,
#'     rsp = 0.1,
#'     eos = 250,
#'     rau = 0.1)
#' t    <- seq(1, 365, 8)
#' tout <- seq(1, 365, 1)
#' y <- fFUN(par, t)
#'
#' l_param <- init_param(y, t)
#' @export
init_param <- function(y, t, w){
    if (any(is.na(y)))
        stop("NA in the time series are not allowed: fill them with e.g. na.approx()")
    if (missing(w)) w <- rep(1, length(y))

    # t      <- t - t[1] # (20190103) seq_along(y)  
    w_min  <- 0.5 # weights greater than w_min are treated as good values.
    # fixed 2018-07-25, If have no enough good points, then set w_min=0
    if (sum(w >= w_min)/length(y) < .4) w_min <- 0

    mx     <- max(y[w >= w_min], na.rm = TRUE)
    mn     <- min(y[w >= w_min], na.rm = TRUE)
    avg    <- mean(y, na.rm = TRUE)

    doy.mx <- t[which.max(y)]
    # fixed 06 March, 2018; Avoid doy.mx not in the range of doy
    # doy    <- quantile(t, c(0.25, 0.75), na.rm = TRUE),
    doy  <- c((doy.mx + first(t))/2, (last(t) + doy.mx) /2)
    t1   <- (doy.mx - doy[1])/3 # adjust for doy[1]
    t2   <- (doy[2] - doy.mx)/3 # adjust for doy[2]

    # if (doy[1] >= doy.mx) doy[1] <- (doy.mx - first(t))/2 + first(t)
    # if (doy[2] <= doy.mx) doy[2] <- (last(t) - doy.mx) /2 + doy.mx
    ampl   <- mx - mn
    deltaY <- ampl*0.1
    half   <- (max(t) - min(t))/2
    deltaT <- half/4
    
    k      <- 4/half*2.67 #approximate value
    # k limits: about 0.004 - 0.2
    # kmin <- 4 / (half * 5), half = 200, k = 0.004
    # kmax <- 4 / (half / 5), half = 100, k = 0.2

    # parameters limit
    lims = list(
        t0  = c(doy.mx - deltaT, doy.mx + deltaT),
        mn  = c(mn - deltaY    , mn + deltaY),
        mx  = c(mx - deltaY*2  , mx + deltaY*2),
        r   = c(k/1.2, k*5),
        sos = c(min(t)         , doy.mx + deltaT),
        eos = c(doy.mx - deltaT, max(t))
    )
    res <- listk( mx, mn, ampl, doy, doy.mx,
        deltaT, deltaY, half, t1, t2,
        k, lims)
    res
}
