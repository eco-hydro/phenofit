# .normalize and .backnormalize function referenced by phenopix package
.normalize     <- function(x, sf) (x-sf[1])/(sf[2]-sf[1])
.backnormalize <- function(x, sf) (x+sf[1]/(sf[2]-sf[1]))*(sf[2]-sf[1])

#' check_input
#'
#' If just replace NA values with a const missing value, it is inappropriate for
#' middle growing season points. If interpolating all values by na.approx, it is
#' unsuitable for large number continous missing segments, e.g. in the start
#' or end of growing season. So a maxgap parameter with a default value of
#' `nptperyear/4` will be a good choice.
#' @param t A numeric vector
#' @param y A numeric vector
#' @param w A numeric vector
#' @param missval Double, which is used to replace NA values in y. If missing,
#' the default vlaue is `ylu[1] - diff(ylu)/10`.
#' @param maxgap Integer, nptperyear/4 will be a suitable value. If continuous
#' missing value numbers less than maxgap, then interpolate those NA values by
#' zoo::na.approx; If false, then replace those NA values with a const value,
#' `ylu[1]`
#'
#' @export
check_input <- function(t, y, w, missval, maxgap = 10, trim = TRUE, alpha = 0.01){
    # alpha <- 0.02
    ylu   <- quantile(y[w == 1], c(alpha/2, 1 - alpha), na.rm = T) #only consider good value
    # ylu    <- range(y, na.rm = T)
    ylu[1] <- pmax(ylu[1], 0)
    # NAN values check
    if (missing(missval)){
        missval <- ylu[1] #- diff(ylu)/10
        # missval <- -1
    }
    if (missing(w)) w <- rep(1, length(y))

    if (trim){
        I_trim    <- y < ylu[1] | y > ylu[2]
        w[I_trim] <- 0
    }
    w[is.na(y)] <- 0
    y <- na.approx(y, maxgap = maxgap, na.rm = FALSE)#just replaced with missval performance bad
    # if still have na values, whether should use mapgap to constrain interpolating range?
    y[is.na(y)] <- missval
    list(t = t, y = y, w = w, ylu = ylu)#quickly return
}

#' check_fit
#' @param ylu limits of y value, [min, max]
#' @export
check_fit <- function(yfit, ylu){
    I_max <- yfit > ylu[2]
    I_min <- yfit < ylu[1]
    yfit[I_max] <- ylu[2]
    yfit[I_min] <- ylu[1]
    return(yfit)
}

#' values out of ylu, set to be na and interpolate it.
check_fit2 <- function(y, ylu){
    I <- which(y < ylu[1] | y > ylu[2])
    if (length(I) > 0){
        n    <-length(y)
        y[I] <- NA
        y <- na.approx(y, na.rm = F)
        # if still have na values in y
        I_nona <- which(!is.na(y)) # not NA id
        if (length(I_nona) != n){
            # na values must are in tail or head now
            iBegin <- first(I_nona)
            iEnd   <- last(I_nona)
            if (iBegin > 2) y[1:iBegin] <- y[iBegin]
            if (iEnd < n)   y[iEnd:n]   <- y[iEnd]
        }
    }
    return(y)
}

# https://rpubs.com/aaronsc32/bisection-method-r
bisection <- function(f, a, b, tolx = 1e-7, toly = 1e-3, maxit = 1000, trace = T,
    ln = TRUE) {
    # mid  <- (a+b)/2
    mid  <- ifelse(ln, exp( (log(a) + log(b))/2), (a + b)/2)
    fmid <- f(mid)
    y0   <- f(a)
    y1   <- f(b)

    i <- 1
    if(sign(y0) == sign(y1)){
        stop("Starting vaules are not unsuitable!")
    } else {
        while(abs(f(mid)) > toly){
            if (i > maxit) stop('Max iteration reached!')

            if(sign(y1) == sign(ymid)){b <- mid}else{a <- mid}
            mid  <- ifelse(ln, exp( (log(a) + log(b))/2), (a + b)/2)
            ymid <- f(mid)
            y0   <- f(a)
            y1   <- f(b)
            #Print the value obtained in each iteration next line
            if (trace) cat(sprintf('i = %3d: fmid = %.4fn\n', i, fmid))
            i<-i+1
        }
    }
    return(mid)
}

# ' Weigthed Whittaker Smoother
# '
# ' With the help of sparseMatrix in Matrix package
# ' @references
# ' Eilers, P.H., 2003. A perfect smoother. Analytical chemistry, 75(14), pp.3631-3636.
# ' @import Matrix
# ' @export
# ' whittaker smoother


# 
#' using cubic spline function to avoid the difficult in setting parameter
#' lambda in smooth.spline
#'
#' cubic spline is inappropriate for daily inputs. Its smooth is not enough.
#' On the contrary, smooth.spline with a low freedom can smooth well.
#'
#' @export
splinefit <- function(x, t = index(x), tout = t, plot = FALSE, df.factor = 0.06, ...){
    # xpred.out <- spline(t, x, xout = tout)$y %>% zoo(., tout)
    n <- length(x)
    # if n < 40, means x was satellite VI
    # if n > 40, means daily data
    df.factor <- ifelse (n <= 46, 1/3, df.factor)
    freedom   <- pmax(df.factor * n, 15)
    fit       <- smooth.spline(t, x, df = freedom)
    xpred.out <- predict(fit, tout)$y %>% zoo(., tout)
    structure(list(data = list(x = x, t = t),
        pred = xpred.out, par = NULL, fun = NULL), class = "phenofit")
}

#' @export
FitDL.Zhang <- function(x, t = index(x), tout = t, optimFUN = I_optimx,
                        method = 'nlm', ...){
    # print("debug")
    e <- FIT_check(x, t)
    deltaY <- ampl/3
    deltaT <- (max(t) - min(t))/4
    k      <- deltaY/deltaT
    k_max  <- 2

    FUN    <- doubleLog.zhang
    prior  <- rbind(
        c(doy.mx, mn, mx, doy[1], k*2, doy[2], k*2),
        c(doy.mx + deltaT/2, mn, mx, doy[1], k*3, doy[2], k*3),
        c(doy.mx, mn, mx, doy[1], k*4, doy[2], k*4))

    lower  <- c(doy.mx - deltaT, mn - deltaY, mx - deltaY, min(t), 0 , doy.mx-deltaT, 0)
    upper  <- c(doy.mx + deltaT, mn + deltaY, mx + deltaY, doy.mx+deltaT, k_max, max(t),   k_max)

    optim_pheno(prior, FUN, x, t, tout, optimFUN, method, lower = lower, upper = upper, ...)#quick return
}

#' @export
FitAG <- function(x, t = index(x), tout = t, FUN, optimFUN = I_optimx,
    method = 'nlminb', ...){
    e <- FIT_check(x, t)
    # print(ls.str(envir = e))
    deltaY     <- ampl/3
    deltaT     <- (max(t) - min(t))/3
    half       <- (max(t) - min(t))/2

    FUN <- doubleAG
    prior <- rbind(
        c(doy.mx, mn, mx, half, 2, half, 2),
        # c(doy.mx, mn, mx, 0.2*half, 1  , 0.2*half, 1),
        # c(doy.mx, mn, mx, 0.5*half, 1.5, 0.5*half, 1.5),
        c(doy.mx, mn, mx, 0.8*half, 3  , 0.8*half, 3))
    # referenced by TIMESAT
    lower      <- c(doy.mx - deltaT, mn - deltaY, mx - deltaY, 0.1*half, 2, 0.1*half, 2)
    upper      <- c(doy.mx + deltaT, mn + deltaY, mx + deltaY, 1.4*half, 8, 1.4*half, 8)

    optim_pheno(prior, FUN, x, t, tout, optimFUN, method, lower = lower, upper = upper, ...)#quick return
}

#'
#' Fitting double logistics, asymmetric gaussian functions
#'
#' @param x input vegetation index time-series.
#' @param t the corresponding doy(day of year) of x.
#' @param tout the output curve fitting time-series time steps.
#' @param optimFUN optimization function to solve curve fitting functions'
#' parameters. It's should be `optimx_fun`, or `optim_p`.
#' @param method method passed to `optimx` or `optim` function.
#' @param ... other paraters passed to optimFUN, such as weights.
#'
#' @return list(pred, par, fun)
#'
#' @examples
#' FitDL.Beck  (x, t, tout, optimFUN = optim_p, pfun = p_nlminb)
#' FitDL.Elmore(x, t, tout, optimFUN = optim_p, pfun = p_nlminb)
#' FitDL.Gu    (x, t, tout, optimFUN = optim_p, pfun = p_nlminb)
#' FitDL.Klos  (x, t, tout, optimFUN = optim_p, pfun = p_optim, method = 'BFGS')
#'
#' FitDL.Zhang (x, t, tout, optimFUN = optim_p, pfun = p_nlm)
#' FitAG (x, t, tout, optimFUN = optim_p, pfun = p_nlminb)
#' @export
FitDL.Beck <- function(x, t = index(x), tout = t, optimFUN = I_optimx,
    method = 'nlminb', ...) {
    if (any(is.na(x)))
        stop("NA in the time series are not allowed: fill them with e.g. na.approx()")

    n    <- length(x)
    avg  <- mean(x, na.rm = TRUE)
    mx   <- max(x, na.rm = TRUE)
    mn   <- min(x, na.rm = TRUE)
    ampl <- mx - mn
    # doy  <- quantile(t, c(0.25, 0.75), na.rm = TRUE)
    doy.mx <- t[which.max(x)]
    doy <- c((doy.mx - first(t))/2 + first(t),
             (last(t) - doy.mx) /2 + doy.mx)

    FUN   <- doubleLog.beck
    deltaY <- ampl/3
    deltaT <- (max(t) - min(t))/3

    prior <- rbind(
        c(mn, mx, doy[1], 0.5, doy[2], 0.5),
        c(mn, mx, doy[1]+deltaT, 0.8, doy[2]-deltaT, 0.8))

    # mn + (mx - mn)*(1/(1 + exp(-rsp*(t - sos))) + 1/(1 + exp(rau*(t - eos))))
    # attr(doubleLog.beck, 'par') <- c("mn", "mx", "sos", "rsp", "eos", "rau")
    # attr(doubleLog.beck, 'formula') <- expression(mn + (mx - mn)*(1/(1 + exp(-rsp*(t - sos))) + 1/(1 + exp(rau*(t - eos)))))

    lower <- c(mn - deltaY, mx - deltaY, min(t), 0, doy.mx - deltaT, 0)
    upper <- c(mn + deltaY, mx + deltaY, doy.mx + deltaT, 1, max(t), 1)
    optim_pheno(prior, FUN, x, t, tout, optimFUN, method, lower = lower, upper = upper, ...)#return
}

#' @export
FitDL.Elmore <- function(x, t = index(x), tout = t, optimFUN = I_optimx,
    method = 'nlminb', ...) {
    if (any(is.na(x)))
        stop("NA in the time series are not allowed: fill them with e.g. na.approx()")
    n    <- length(x)
    avg  <- mean(x, na.rm = TRUE)
    mx   <- max(x, na.rm = TRUE)
    mn   <- min(x, na.rm = TRUE)
    ampl <- mx - mn
    doy  <- quantile(t, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
    doy.mx <- t[which.max(x)]
    # doy <- c((doy.mx - first(t))/2 + first(t),
    #          (last(t) - doy.mx) /2 + doy.mx)

    half   <- (max(t) - min(t))/2

    FUN   <- doubleLog.elmore
    prior <- rbind(
        c(mn, mx - mn, doy[2], half*0.1, doy[4], half*0.1, 0.002),
        c(mn, mx - mn, doy[2], half*0.2, doy[5], half*0.2, 0.002),
        c(mn, mx - mn, doy[1], half*0.5, doy[4], half*0.5, 0.05),
        c(mn, mx - mn, doy[1], half*0.8, doy[5], half*0.8, 0.1))
    # xpred <- m1 + (m2 - m7*t)*((1/(1 + exp((m3l - t)/m4l))) - (1/(1 + exp((m5l - t)/m6l))))
    deltaY <- ampl/3
    deltaT <- (max(t) - min(t))/3
    lower  <- c(mn - deltaY, mx - deltaY, min(t), 0  , doy.mx-deltaT,   0, 0)
    upper  <- c(mn + deltaY, mx + deltaY, doy.mx+deltaT, 200, max(t), 200, Inf)

    optim_pheno(prior, FUN, x, t, tout, optimFUN, method, lower = lower, upper = upper, ...)#return
}

#  @reference https://github.com/kongdd/phenopix/blob/master/R/FitDoubleLogGu.R
#' @export
FitDL.Gu <- function(x, t = index(x), tout = t, optimFUN = I_optimx, method = "nlminb", ...) {
    if (any(is.na(x)))
        stop("NA in the time series are not allowed: fill them with e.g. na.approx()")

    n    <- length(x)
    avg  <- mean(x, na.rm = TRUE)
    mx   <- max(x, na.rm = TRUE)
    mn   <- min(x, na.rm = TRUE)
    ampl <- mx - mn
    # doy  <- quantile(t, c(0.25, 0.75), na.rm = TRUE)
    doy.mx <- t[which.max(x)]
    doy <- c((doy.mx - first(t))/2 + first(t),
             (last(t) - doy.mx) /2 + doy.mx)

    deltaY <- ampl/3
    deltaT <- (max(t) - min(t))/4

    a1 <- ampl
    a2 <- ampl
    t1 <- doy[1] + 0.5 * (doy.mx - doy[1])
    t2 <- doy.mx + 0.5 * (doy[2] - doy.mx)
    b1 <- 0.1
    b2 <- 0.1
    c1 <- 1
    c2 <- 1

    FUN <- doubleLog.gu
    prior <- rbind(
        c(mn, a1, a2, t1, 0.05, t2, 0.05, 1, 1),
        c(mn, a1, a2, t1,  0.1, t2 , 0.1, 2, 2),
        c(mn, a1, a2, t1,  0.2, t2 , 0.2, 3, 3),
        c(mn, a1, a2, doy[1], 0.5, t2    , 0.5, 0.5, 0.5),
        c(mn, a1, a2, t1    , 0.8, doy[2], 0.8,   5, 5))
    # y0 + (a1/(1 + exp(-(t - t1)/b1))^c1) - (a2/(1 + exp(-(t - t2)/b2))^c2)
    lower <- c(mn - deltaY, mx - deltaY, mx - deltaY, min(t), 0, doy.mx - deltaT, 0, 0, 0)
    upper <- c(mn + deltaY, mx + deltaY, mx + deltaY, doy.mx + deltaT, 1, max(t), 1, Inf, Inf)
    optim_pheno(prior, FUN, x, t, tout, optimFUN, method, lower = lower, upper = upper, ...)#return
}

#' @export
FitDL.Klos <- function(x, t = index(x), tout = t, optimFUN = I_optimx,
    method = 'BFGS', ...) {
    if (any(is.na(x)))
        stop("NA in the time series are not allowed: fill them with e.g. na.approx()")

    n    <- length(x)
    avg  <- mean(x, na.rm = TRUE)
    mx   <- max(x, na.rm = TRUE)
    mn   <- min(x, na.rm = TRUE)
    ampl <- mx - mn
    # doy  <- quantile(t, c(0.25, 0.75), na.rm = TRUE)
    doy.mx <- t[which.max(x)]
    doy <- c((doy.mx - first(t))/2 + first(t),
             (last(t) - doy.mx) /2 + doy.mx)

    a1 <- 0
    a2 <- 0  #ok
    b1 <- mn #ok
    b2 <- 0  #ok
    c  <- 0.2 * max(x)  # ok
    ## very slightly smoothed spline to get reliable maximum
    # tmp <- smooth.spline(x, df = 0.5 * length(x))#, find error: 20161104, fix tomorrow

    B1 <- 4/(doy.mx - doy[1])
    B2 <- 3.2/(doy[2] - doy.mx)
    m1 <- doy[1] + 0.5 * (doy.mx - doy[1])
    m2 <- doy.mx + 0.5 * (doy[2] - doy.mx)
    m1.bis <- doy[1]
    m2.bis <- doy[2]
    q1 <- 0.5  # ok
    q2 <- 0.5  # ok
    v1 <- 2    # ok
    v2 <- 2    # ok

    FUN <- doubleLog.klos
    prior <- rbind(
        c(a1, a2, b1, b2, c, B1, B2, m1, m2, q1, q2, v1, v2),
        c(a1, a2, b1, 0.01, 0, B1, B2, m1, m2.bis, q1, 1, v1, 4),
        c(a1, a2, b1, b2, c, B1, B2, m1.bis, m2, q1, q2, v1, v2),
        c(a1, a2, b1, b2, c, B1, B2, m1, m2.bis, q1, q2, v1, v2),
        c(a1, a2, b1, b2, c, B1, B2, m1.bis, m2, q1, q2, v1, v2))
    optim_pheno(prior, FUN, x, t, tout, optimFUN, method, ...)#quickly return
}

#' @export
FIT_check <- function(x, t){
    if (any(is.na(x)))
        stop("NA in the time series are not allowed: fill them with e.g. na.approx()")
    mx   = max(x, na.rm = TRUE)
    mn   = min(x, na.rm = TRUE)

    doy.mx <- t[which.max(x)]
    # fixed 06 March, 2018; Avoid doy.mx not in the range of doy
    # doy    <- quantile(t, c(0.25, 0.75), na.rm = TRUE),
    doy <- c((doy.mx - first(t))/2 + first(t),
             (last(t) - doy.mx) /2 + doy.mx)
    # if (doy[1] >= doy.mx) doy[1] <- (doy.mx - first(t))/2 + first(t)
    # if (doy[2] <= doy.mx) doy[2] <- (last(t) - doy.mx) /2 + doy.mx
    
    list2env(list(
        n    = length(x),
        avg  = mean(x, na.rm = TRUE),
        mx   = mx,
        mn   = mn,
        ampl = mx - mn,
        doy  = doy,
        doy.mx = doy.mx
    ), envir = parent.frame())
}

