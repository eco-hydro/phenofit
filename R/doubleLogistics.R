#' Double logistics functions
#'
#' Define double logistics, piecewise logistics and many other functions to
#' curve fit VI time-series
#' \itemize{
#'    \item \code{Logistic} The traditional simplest logistic function. It can
#'      be only used in half growing season, i.e. vegetation green-up or senescence
#'      period.
#'    \item \code{doubleLog.Zhang} Piecewise logistics, Zhang Xiaoyang, RME, 2003.
#'    \item \code{doubleAG} Asymmetric Gaussian.
#'    \item \code{doubleLog.Beck} Beck logistics.
#'    \item \code{doubleLog.Gu} Gu logistics.
#'    \item \code{doubleLog.Elmore} Elmore logistics.
#'    \item \code{doubleLog.Klos} Klos logistics.
#' }
#'
#' All of those function have \code{par} and \code{formula} attributes for the
#' convenience for analytical D1 and D2
#' 
#' @param par A vector of parameters
#' @param t A \code{Date} or numeric vector
#' @references
#' Peter M. Atkinson, et al., 2012, RSE, 123:400-417
#'
#' @rdname logistics
#' @export
Logistic <- function(par, t){
    mn   <- par[1]
    mx   <- par[2]
    sos  <- par[3]
    rsp  <- par[4]
    pred <- (mx - mn)/(1 + exp(-rsp*(t - sos))) + mn
    # pred <- c/(1 + exp(a + b * t)) + d
    return(pred)
}
attr(Logistic, 'name')    <- 'Logistic'
attr(Logistic, 'par')     <- c("mn", "mx", "sos", "rsp")
attr(Logistic, 'formula') <- expression((mx - mn)/(1 + exp(-rsp*(t - sos))) + mn)

#' @rdname logistics
#' @export
doubleLog.Zhang <- function(par, t){
    t0  <- par[1]
    mn  <- par[2]
    mx  <- par[3]
    sos <- par[4]
    rsp <- par[5]
    eos <- par[6]
    rau <- par[7]

    if (t0 - sos <= 1 || t0 - eos >= -1) return(rep(NA, length(t)))
    # In order to make sure the shape of S curve, should be satisfy:
    # t0 < eos, t0 > sos

    # xpred1 <- (mx - mn)/(1 + exp(-rsp*(t - sos))) + mn
    # xpred2 <- (mx - mn)/(1 + exp(rau*(t - eos))) + mn
    # pred  <- xpred1*(t <= t0) + xpred2*(t > t0)

    # above expressions cost 1.5 times as the below
    pred <- (mx - mn)/c(1 + exp(-rsp*(t[t <= t0] - sos)),
                         1 + exp( rau*(t[t >  t0] - eos))) + mn
    return(pred)
}
attr(doubleLog.Zhang, 'name')    <- 'doubleLog.Zhang'
attr(doubleLog.Zhang, 'par')     <- c("t0", "mn", "mx", "sos", "rsp", "eos", "rau")
# piecewise function
attr(doubleLog.Zhang, 'formula') <- expression( (mx - mn)/(1 + exp(-rsp*(t - sos))),
                                                (mx - mn)/(1 + exp( rau*(t - eos))) )

#' @rdname logistics
#' @export
doubleAG <- function(par, t){
    t0  <- par[1]
    mn  <- par[2]
    mx  <- par[3]
    rsp <- par[4]
    a3  <- par[5]
    rau <- par[6]
    a5  <- par[7]

    pred <- mn + (mx - mn)*exp(- c( ((t0 - t[t <= t0])*rsp) ^a3,
                    ((t[t >  t0] - t0)*rau) ^a5) )
    return(pred)
    # pred <- (mx - mn)/c(1 + exp(-rsp*(t[t <= t0] - sos)),
    #     1 + exp(-rau*(t[t >  t0] - eos))) + mn
}
# a3, a5 should be greater than 1
attr(doubleAG, 'name')    <- 'doubleAG'
attr(doubleAG, 'par')     <- c("t0", "mn", "mx", "rsp", "a3", "rau", "a5")
attr(doubleAG, 'formula') <- expression( mn + (mx - mn)*exp(- ((t0 - t)*rsp) ^a3 ),
                                         mn + (mx - mn)*exp(- ((t - t0)*rau) ^a5 ))

#' @rdname logistics
#' @export
doubleLog.Beck <- function(par, t) {
    mn  <- par[1]
    mx  <- par[2]
    sos <- par[3]
    rsp <- par[4]
    eos <- par[5]
    rau <- par[6]
    # if (sos >= eos) return(rep(9999, length(t)))
    try(if (eos < sos) return(rep(99, length(t))), silent = T)

    pred <- mn + (mx - mn)*(1/(1 + exp(-rsp*(t - sos))) + 1/(1 + exp(rau*(t - eos))) - 1)
    return(pred)
}
attr(doubleLog.Beck, 'name') <- 'doubleLog.Beck'
attr(doubleLog.Beck, 'par') <- c("mn", "mx", "sos", "rsp", "eos", "rau")
attr(doubleLog.Beck, 'formula') <- expression(mn + (mx - mn)*(1/(1 + exp(-rsp*(t - sos))) + 1/(1 + exp(rau*(t - eos)))) - 1)

#' @rdname logistics
#' @export
doubleLog.Elmore <- function(par, t) {
    mn  <- par[1]
    mx  <- par[2]
    # m3  <- par[3]
    # m4  <- par[4]
    # m5  <- par[5]
    # m6  <- par[6]
    # m3l <- m3/m4
    # m4l <- 1/m4
    # m5l <- m5/m6
    # m6l <- 1/m6
    sos  <- par[3] # SOS
    rsp  <- par[4] # 1/rsp
    eos  <- par[5] # EOS
    rau  <- par[6] # 1/rau
    m7   <- par[7]
    pred <- mn + (mx - m7*t)*( 1/(1 + exp(-rsp*(t-sos))) - 1/(1 + exp(-rau*(t-eos))) )
    return(pred)
}

attr(doubleLog.Elmore, 'name')    <- 'doubleLog.Elmore'
attr(doubleLog.Elmore, 'par')     <- c("mn", "mx", "sos", "rsp", "eos", "rau", "m7")
attr(doubleLog.Elmore, 'formula') <- expression( mn + (mx - m7*t)*( 1/(1 + exp(-rsp*(t-sos))) - 1/(1 + exp(-rau*(t-eos))) ) )

#' @rdname logistics
#' @export
doubleLog.Gu <- function(par, t) {
    y0  <- par[1]
    a1  <- par[2]
    a2  <- par[3]
    sos <- par[4]
    rsp <- par[5]
    eos <- par[6]
    rau <- par[7]
    # t1  <- par[4]
    # t2  <- par[5]
    # b1  <- par[6]
    # b2  <- par[7]
    c1  <- par[8]
    c2  <- par[9]
    # pred <- y0 + (a1/(1 + exp(-(t - t1)/b1))^c1) - (a2/(1 + exp(-(t - t2)/b2))^c2)
    pred <- y0 + (a1/(1 + exp(-rsp*(t - sos)))^c1) - (a2/(1 + exp(-rau*(t - eos)))^c2)
    return(pred)
}
attr(doubleLog.Gu, 'name')    <- 'doubleLog.Gu'
attr(doubleLog.Gu, 'par')     <- c('y0', 'a1', 'a2', 'sos', 'rsp', 'eos', 'rau', 'c1', 'c2')
attr(doubleLog.Gu, 'formula') <- expression(y0 + (a1/(1 + exp(-rsp*(t - sos)))^c1) - (a2/(1 + exp(-rau*(t - eos)))^c2))

#' @rdname logistics
#' @export
doubleLog.Klos <- function(par, t) {
    a1 <- par[1]
    a2 <- par[2]
    b1 <- par[3]
    b2 <- par[4]
    c  <- par[5]
    B1 <- par[6]
    B2 <- par[7]
    m1 <- par[8]
    m2 <- par[9]
    q1 <- par[10]
    q2 <- par[11]
    v1 <- par[12]
    v2 <- par[13]
    pred <- (a1*t + b1) + (a2*t^2 + b2*t + c) * (1/(1 + q1 * exp(-B1 * (t - m1)))^v1
        - 1/(1 + q2 * exp(-B2 * (t - m2)))^v2)
    return(pred)
}
attr(doubleLog.Klos, 'name')    <- 'doubleLog.Klos'
attr(doubleLog.Klos, 'par')     <- c('a1', 'a2', 'b1', 'b2', 'c', 'B1', 'B2',
    'm1', 'm2', 'q1', 'q2', 'v1', 'v2')
attr(doubleLog.Klos, 'formula') <- expression((a1*t + b1) + (a2*t^2 + b2*t + c) * (1/(1 + q1 * exp(-B1 * (t - m1)))^v1
        - 1/(1 + q2 * exp(-B2 * (t - m2)))^v2))

.qr.solve <- function(a, b, tol = 1e-07, LAPACK = TRUE) {
    if (!is.qr(a)) a <- qr(a, tol = tol, LAPACK = LAPACK)
    nc <- ncol(a$qr)
    nr <- nrow(a$qr)
    if (a$rank != min(nc, nr)) stop("singular matrix 'a' in solve")
    if (missing(b)) {
        if (nc != nr) stop("only square matrices can be inverted")
        b <- diag(1, nc)
    }
    res <- qr.coef(a, b)
    res[is.na(res)] <- 0
    res
}
# vc <- .qr.solve(opt$hessian)
# npar <- nrow(vc)
# s2 <- opt.df$cost[best]^2 / (n - npar)
# std.errors <- sqrt(diag(vc) * s2)     # standard errors
# return: stdError=std.error

#' Common goal function of those curve fitting methods
#'
#' @inheritParams optim_pheno
#' @inheritParams Logistic
#'
#' @param fun A curve fitting function, see \code{\link{Logistic}}.
#' @param ... other parameters will be ignored.
#'
#' @return RMSE root mean square error of curve fitting values.
#' @export
f_goal <- function(
    par, y, t,
    fun = c(doubleLog.Elmore,
           doubleLog.Beck,
           doubleLog.Klos,
           doubleLog.Gu,
           doubleLog.Zhang,
           doubleAG),
    w, ylu, ...) {
    # FUN <- match.fun(fun)
    if (!all(is.finite(par))) return(9999)

    pred <- fun(par, t = t)
    # If have no finite values, return 9999
    if (!all(is.finite(pred))) return(9999) # for Klos fitting
    if (missing(w)) w <- rep(1, length(y))

    if (!missing(ylu)){
        # points out of ylu should be punished!
        w[pred < ylu[1] | pred > ylu[2]] <- 0
        # pred   <- check_fit(pred, ylu)
    }
    SSE  <- sum((y - pred)^2 * w)
    RMSE <- sqrt(SSE/length(y))
    NSE  <- sum((y - pred)^2 * w)/sum((y - mean(pred))^2)

    # 1. better handle low and high values simulation
    # xpred_2 <- sqrt(xpred_2)
    # x_2     <- sqrt(x_2)
    # xpred_2 <- log(xpred_2+1)
    # x_2     <- log(x_2+1)
    # xpred_2 <- 1/pred          # inverse NSE
    # x_2     <- 1/y

    # xpred_2 <- pred - mean(y)
    # x_2     <- y - mean(y)
    # NSE2 <- sum((x_2 - xpred_2)^2 * w)/sum((x_2 - mean(x_2))^2) #NSE

    # const <- ylu[2]
    # xpred_2 <- pred - ylu[1]; xpred_2[xpred_2 < 0] <- const
    # x_2     <- y     - ylu[1]; x_2[x_2 < 0] <- const
    return(RMSE)
}

# attach gradient and hessian analytical function to curve fitting functions
.dls <- lapply(
    c("doubleLog.Beck", "doubleLog.Elmore", "doubleLog.Gu",
      "doubleLog.Klos", "doubleLog.Zhang", "doubleAG"),
    function (FUN){
        # FUN <- deparse(substitute(fun))
        fun <- get(FUN)
        attr(fun, 'gradient') <- gradf_t(fun) # gradient
        attr(fun, 'hessian')  <- hessf_t(fun) # hessian
        # print(environment(fun))
        # print(fun)
        # print(FUN)
        assign(FUN, fun, envir = environment(fun)) #environment("namespace:phenofit"))#
        # fun
    })