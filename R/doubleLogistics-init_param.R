#' init_param
#' 
#' Initialize parameters of double logistic function
#' 
#' @inheritParams check_input
#' @param e The object returned by [init_param()]
#' @param type integer, `1` or `-1`
#' -  `1`: trough-to-trough curve fitting
#' - `-1`: peak-to-peak curve fitting
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
init_param <- function(y, t, w, type = 1L){
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

    doy.mx = if (type == 1L) t[which.max(y)] else t[which.min(y)]
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
    # constrain parameter in a reasonable range
    tmax = max(t)
    lims = list(
        t0  = c(doy.mx - deltaT, doy.mx + deltaT),
        mn  = c(mn - deltaY    , mn + deltaY),
        mx  = c(mx - deltaY*2  , mx + deltaY*2),
        r   = c(k/1.2, k*5),
        # r   = sort(c(k/1.2, k*5) * type),                                     # >= v0.3.5
        sos = c(min(t)         , pmin(doy.mx + deltaT), tmax),
        eos = c(doy.mx - deltaT, tmax)
    )
    # plot(t, y, main = "init_param")
    if (type == -1) {
        # swap value mx and mn, this is suit for AG and Zhang
        tmp = mx
        mx  = mn
        mn  = tmp

        tmp     = lims$mx
        lims$mx = lims$mn
        lims$mn = tmp
    }
    res <- listk( mx, mn, ampl, doy, doy.mx,
        deltaT, deltaY, half, t1, t2,
        k, lims)
    res
}

#' @rdname init_param
#' @export
init_Zhang <- function(e, type = 1L, ...) {
    # if (missing(w)) w <- rep(1, length(y))
    # e <- init_param(y, t, w, type = type)
    prior <- with(e, rbind(
        c(doy.mx, mn, mx, doy[1], k, doy[2], k),
        c(doy.mx, mn, mx, doy[1] + t1, k * 2, doy[2] - t1, k * 2),
        c(doy.mx, mn, mx, doy[1] - t1, k, doy[2] + t2, k)
    ))
    param_lims <- e$lims[c("t0", "mn", "mx", "sos", "r", "eos", "r")]
    lower <- sapply(param_lims, `[`, 1)
    upper <- sapply(param_lims, `[`, 2)
    # lower[["r"]] %>% multiply_by(1/3)
    # upper[["r"]] %>% multiply_by(3)
    listk(prior, lower, upper)
}

#' @rdname init_param
#' @export
init_AG <- function(e, type = 1L, ...) {
    # if (missing(w)) w <- rep(1, length(y))
    # e <- init_param(y, t, w, type = type)
    prior <- with(e, rbind(
        c(doy.mx, mn, mx, 1 / half, 2, 1 / half, 2),
        # c(doy.mx, mn, mx, 0.2*half, 1  , 0.2*half, 1),
        # c(doy.mx, mn, mx, 0.5*half, 1.5, 0.5*half, 1.5),
        c(doy.mx, mn, mx, 1 / (0.8 * half), 3, 1 / (0.8 * half), 3)
    ))
    # c("t0", "mn", "mx", "rsp", "a3", "rau", "a5")
    # referenced by TIMESAT
    lower <- with(e$lims, c(t0[1], mn[1], mx[1], 1 / (1.4 * e$half), 2, 1 / (1.4 * e$half), 2))
    upper <- with(e$lims, c(t0[2], mn[2], mx[2], 1 / (0.1 * e$half), 6, 1 / (0.1 * e$half), 6))
    listk(prior, lower, upper)
}

#' @rdname init_param
#' @export
init_AG2 <- function(e, type = 1L, ...) {
    # if (missing(w)) w <- rep(1, length(y))
    # e <- init_param(y, t, w, type = type)
    prior <- with(e, rbind(
        c(doy.mx, mn, mn, mx, 1 / half, 2, 1 / half, 2),
        # c(doy.mx, mn, mx, 0.2*half, 1  , 0.2*half, 1),
        # c(doy.mx, mn, mx, 0.5*half, 1.5, 0.5*half, 1.5),
        c(doy.mx, mn, mn, mx, 1 / (0.8 * half), 3, 1 / (0.8 * half), 3)
    ))
    # referenced by TIMESAT
    lower <- with(e$lims, c(t0[1], mn[1], mn[1], mx[1], 1 / (1.4 * e$half), 2, 1 / (1.4 * e$half), 2))
    upper <- with(e$lims, c(t0[2], mn[2], mn[2], mx[2], 1 / (0.1 * e$half), 6, 1 / (0.1 * e$half), 6))
    listk(prior, lower, upper)
}

#' @rdname init_param
#' @export
init_Beck <- function(e, type = 1L, ...) {
    # if (missing(w)) w <- rep(1, length(y))
    # e <- init_param(y, t, w, type = type)
    prior <- with(e, rbind(
        c(mn, mx, doy[1], k, doy[2], k),
        c(mn, mx, doy[1] + t1, k * 2, doy[2] - t2, k * 2)
    ))
    param_lims <- e$lims[c("mn", "mx", "sos", "r", "eos", "r")]
    lower <- sapply(param_lims, `[`, 1)
    upper <- sapply(param_lims, `[`, 2)
    listk(prior, lower, upper)
}

#' @rdname init_param
#' @export
init_Elmore <- function(e, type = 1L, ...) {
    # if (missing(w)) w <- rep(1, length(y))
    # e <- init_param(y, t, w, type = type)
    # doy_q  <- quantile(t, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
    # TODO: remove bad one
    # generally m7 < 0.001
    prior <- with(e, rbind(
        c(mn, mx - mn, doy[1] + t1, k * 2.5, doy[2] - t2, k * 2.5, 0.002),
        # c(mn, mx - mn, doy[1]   , k*1.25 , doy[2]   , k*1.25 , 0.002),
        # c(mn, mx - mn, doy[1]   , k*1    , doy[2]   , k      , 0.001),
        c(mn, mx - mn, doy[1] - t1, k * 0.25, doy[2] + t2, k * 0.25, 0.001)
    ))
    # c("mn", "mx", "sos", "rsp", "eos", "rau", "m7")
    # TODO: m7 might need to be adjust
    param_lims <- e$lims[c("mn", "mx", "sos", "r", "eos", "r")]
    lower <- c(sapply(param_lims, `[`, 1), 0)
    upper <- c(sapply(param_lims, `[`, 2), 1)
    listk(prior, lower, upper)
}

#' @rdname init_param
#' @export
init_Gu <- function(e, type = 1L, ...) {
    # if (missing(w)) w <- rep(1, length(y))
    # e <- init_param(y, t, w, type = type)
    a <- e$ampl * type
    b1 <- 0.1
    b2 <- 0.1
    c1 <- 1
    c2 <- 1
    prior <- with(e, rbind(
        # c(mn, a, a, doy[1]-t1, k/2 , doy[2]+t2, k/2, 1  , 1),
        c(mn, a, a, doy[1], k, doy[2], k, 2, 2),
        # c(mn, a, a, doy[1]   , k*2 , doy[2]   , k*2, 3  , 3),
        c(mn, a, a, doy[1] + t1, k * 2, doy[2] - t2, k * 2, 0.5, 0.5),
        c(mn, a, a, doy[1] + t1, k * 3, doy[2] - t2, k * 3, 5, 5)
    ))

    # y0 + (a1/(1 + exp(-(t - t1)/b1))^c1) - (a2/(1 + exp(-(t - t2)/b2))^c2)
    # c('y0', 'a1', 'a2', 'sos', 'rsp', 'eos', 'rau', 'c1', 'c2')
    if (type == 1L) {
        param_lims <- e$lims[c("mn", "mx", "mx", "sos", "r", "eos", "r")]
    } else if (type == -1L) {
        param_lims <- e$lims[c("mn", "mn", "mn", "sos", "r", "eos", "r")]
        # mn and mx had been wrapped
        A_lims = rev(-param_lims[["mn"]])
        param_lims[[2]] = A_lims # a1
        param_lims[[3]] = A_lims # a2
    }
    
    lower <- c(sapply(param_lims, `[`, 1), 0, 0)
    upper <- c(sapply(param_lims, `[`, 2), Inf, Inf)
    listk(prior, lower, upper)
}

#' @rdname init_param
#' @export
init_Klos <- function(e, type = 1L, ...) {
    # if (missing(w)) w <- rep(1, length(y))
    # e <- init_param(y, t, w, type = type)
    a1 <- 0
    a2 <- 0    # ok
    b1 <- e$mn # ok
    b2 <- 0    # ok
    c <- 0.2 * e$mx # ok
    ## very slightly smoothed spline to get reliable maximum
    # tmp <- smooth.spline(y, df = 0.5 * length(y))#, find error: 20161104, fix tomorrow

    prior <- with(e, {
        B1 <- 4 / (doy.mx - doy[1])
        B2 <- 3.2 / (doy[2] - doy.mx)
        m1 <- doy[1] + 0.5 * (doy.mx - doy[1])
        m2 <- doy.mx + 0.5 * (doy[2] - doy.mx)
        m1.bis <- doy[1]
        m2.bis <- doy[2]
        q1 <- 0.5 # ok
        q2 <- 0.5 # ok
        v1 <- 2 # ok
        v2 <- 2 # ok

        rbind(
            c(a1, a2, b1, b2, c, B1, B2, m1, m2, q1, q2, v1, v2),
            c(a1, a2, b1, 0.01, 0, B1, B2, m1, m2.bis, q1, 1, v1, 4),
            c(a1, a2, b1, b2, c, B1, B2, m1.bis, m2, q1, q2, v1, v2),
            c(a1, a2, b1, b2, c, B1, B2, m1, m2.bis, q1, q2, v1, v2),
            c(a1, a2, b1, b2, c, B1, B2, m1.bis, m2, q1, q2, v1, v2)
        )
    })
    listk(prior, lower = -Inf, upper = Inf)
}
