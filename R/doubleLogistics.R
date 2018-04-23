#'
#' Define double logistics, piecewise logistics and many other functions used to
#' curve fit VI time-series
#'
#' @usage
#' doubleLog.beck(par, t)
#' doubleLog.gu(par, t)
#' doubleLog.elmore(par, t)
#' doubleLog.klos(par, t)
#' doubleLog.zhang(par, t) #piecewise
#' doubleAG(par, t)        #piecewise
#'
#' All of those function have `par` and `formula` attributes for the convenience for
#' analytical D1 and D2
#' @export
doubleLog.beck <- function(par, t) {
    mn  <- par[1]
    mx  <- par[2]
    sos <- par[3]
    rsp <- par[4]
    eos <- par[5]
    rau <- par[6]
    # if (sos >= eos) return(rep(9999, length(t)))
    # try(if (eos < sos) return(rep(99, length(t))), silent = T)
    xpred <- mn + (mx - mn)*(1/(1 + exp(-rsp*(t - sos))) + 1/(1 + exp(rau*(t - eos))) - 1)
    return(xpred)
}
attr(doubleLog.beck, 'name') <- 'doubleLog.beck'
attr(doubleLog.beck, 'par') <- c("mn", "mx", "sos", "rsp", "eos", "rau")
attr(doubleLog.beck, 'formula') <- expression(mn + (mx - mn)*(1/(1 + exp(-rsp*(t - sos))) + 1/(1 + exp(rau*(t - eos)))) - 1)

#' @export
doubleLog.gu <- function(par, t) {
    y0  <- par[1]
    a1  <- par[2]
    a2  <- par[3]
    sos  <- par[4]
    eos  <- par[5]
    rsp  <- par[6]
    rau  <- par[7]
    # t1  <- par[4]
    # t2  <- par[5]
    # b1  <- par[6]
    # b2  <- par[7]
    c1  <- par[8]
    c2  <- par[9]
    # xpred <- y0 + (a1/(1 + exp(-(t - t1)/b1))^c1) - (a2/(1 + exp(-(t - t2)/b2))^c2)
    xpred <- y0 + (a1/(1 + exp(-rsp*(t - sos)))^c1) - (a2/(1 + exp(-rau*(t - eos)))^c2)
    return(xpred)
}
attr(doubleLog.gu, 'name')    <- 'doubleLog.gu'
attr(doubleLog.gu, 'par')     <- c('y0', 'a1', 'a2', 'sos', 'eos', 'rsp', 'rau', 'c1', 'c2')
attr(doubleLog.gu, 'formula') <- expression(y0 + (a1/(1 + exp(-rsp*(t - sos)))^c1) - (a2/(1 + exp(-rau*(t - eos)))^c2))

#' @export
doubleLog.elmore <- function(par, t) {
    m1  <- par[1]
    m2  <- par[2]
    # m3  <- par[3]
    # m4  <- par[4]
    # m5  <- par[5]
    # m6  <- par[6]
    # m3l <- m3/m4
    # m4l <- 1/m4
    # m5l <- m5/m6
    # m6l <- 1/m6
    m3l  <- par[3]
    m4l  <- par[4]
    m5l  <- par[5]
    m6l  <- par[6]
    m7   <- par[7]
    xpred <- m1 + (m2 - m7*t)*((1/(1 + exp((m3l - t)/m4l))) - (1/(1 + exp((m5l - t)/m6l))))
    return(xpred)
}

attr(doubleLog.elmore, 'name')    <- 'doubleLog.elmore'
attr(doubleLog.elmore, 'par')     <- c("m1", "m2", "m3l", "m4l", "m5l", "m6l", "m7")
attr(doubleLog.elmore, 'formula') <- expression(m1 + (m2 - m7*t)*((1/(1 + exp((m3l - t)/m4l))) - (1/(1 + exp((m5l - t)/m6l)))))

#' @export
doubleLog.klos <- function(par, t) {
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
    xpred <- (a1*t + b1) + (a2*t^2 + b2*t + c) * (1/(1 + q1 * exp(-B1 * (t - m1)))^v1
        - 1/(1 + q2 * exp(-B2 * (t - m2)))^v2)
    return(xpred)
}
attr(doubleLog.klos, 'name')    <- 'doubleLog.klos'
attr(doubleLog.klos, 'par') <- c('a1', 'a2', 'b1', 'b2', 'c', 'B1', 'B2',
    'm1', 'm2', 'q1', 'q2', 'v1', 'v2')
attr(doubleLog.klos, 'formula') <- expression((a1*t + b1) + (a2*t^2 + b2*t + c) * (1/(1 + q1 * exp(-B1 * (t - m1)))^v1
        - 1/(1 + q2 * exp(-B2 * (t - m2)))^v2))

#' only fit part of growing season NDVI, before or after peak NDVI
#' @export
Log.zhang <- function(par, t){
    mn <- par[1]
    mx <- par[2]
    sos <- par[3]
    rsp <- par[4]
    xpred <- (mx - mn)/(1 + exp(-rsp*(t - sos))) + mn
    # xpred <- c/(1 + exp(a + b * t)) + d
    return(xpred)
}
attr(Log.zhang, 'name')    <- 'Log.zhang'
attr(Log.zhang, 'par')     <- c("mn", "mx", "sos", "rsp")
attr(Log.zhang, 'formula') <- expression((mx - mn)/(1 + exp(-rsp*(t - sos))) + mn)

#'
#' simplest piecewise logistics function
#'
#' fit whole growing season NDVI. First introduced to phenology by
#'  Zhang Xiaoyang, RME, 2003
#'
#' @export
doubleLog.zhang <- function(par, t){
    t0 <- par[1]
    mn <- par[2]
    mx <- par[3]
    sos <- par[4]
    rsp <- par[5]
    eos <- par[6]
    rau <- par[7]

    if (t0 < sos | t0 > eos) return(rep(NA, length(t)))
    # In order to make sure the shape of S curve, should be satisfy:
    # t0 < eos, t0 > sos
    
    # xpred1 <- (mx - mn)/(1 + exp(-rsp*(t - sos))) + mn
    # xpred2 <- (mx - mn)/(1 + exp(rau*(t - eos))) + mn
    # xpred  <- xpred1*(t <= t0) + xpred2*(t > t0)

    # above expressions cost 1.5 times as the below
    xpred <- (mx - mn)/c(1 + exp(-rsp*(t[t <= t0] - sos)),
                         1 + exp(-rau*(t[t >  t0] - eos))) + mn
    return(xpred)
}
attr(doubleLog.zhang, 'name')    <- 'doubleLog.zhang'
attr(doubleLog.zhang, 'par')     <- c("t0", "mn", "mx", "sos", "rsp", "eos", "rau")
# piecewise function
attr(doubleLog.zhang, 'formula') <- expression( (mx - mn)/(1 + exp(-rsp*(t - sos))),
                                                (mx - mn)/(1 + exp(-rau*(t - eos))) )
# global local function
#' Asymmetric Gaussian
#'
#' @references
#' Peter M. Atkinson, et al., 2012, RSE, 123:400-417
#' @export
doubleAG <- function(par, t){
    t0 <- par[1]
    mn <- par[2]
    mx <- par[3]
    a2 <- par[4]
    a3 <- par[5]
    a4 <- par[6]
    a5 <- par[7]

    xpred <- mn + (mx - mn)*exp(- c( ((t0 - t[t <= t0])/a2) ^a3,
                    ((t[t >  t0] - t0)/a4) ^a5) )
    return(xpred)
    # xpred <- (mx - mn)/c(1 + exp(-rsp*(t[t <= t0] - sos)),
    #     1 + exp(-rau*(t[t >  t0] - eos))) + mn
}
# a3, a5 should be greater than 1
attr(doubleAG, 'name')    <- 'doubleAG'
attr(doubleAG, 'par') <- c("t0", "mn", "mx", "a2", "a3", "a4", "a5")
# piecewise function
attr(doubleAG, 'formula') <- expression( mn + (mx - mn)*exp(- ((t0 - t)/a2) ^a3 ),
                                         mn + (mx - mn)*exp(- ((t - t0)/a4) ^a5 ))
# doubleAG_grad <- function(par, t){
#     t0 <- par[1]
#     mn <- par[2]
#     mx <- par[3]
#     a2 <- par[4]
#     a3 <- par[5]
#     a4 <- par[6]
#     a5 <- par[7]

#     .expr1 <- mx - mn
#     .expr3 <- (t0 - t[t < t0])/a2
#     .expr6 <- exp(-.expr3^a3)
#     # .value <- mn + .expr1 * .expr6
#     grad_a <- .expr1 * (.expr6 * (.expr3^(a3 - 1) * (a3 * (1/a2))))

#     .expr3 <- (t[t > t0] - t0)/a4
#     .expr6 <- exp(-.expr3^a4)
#     # .value <- mn + .expr1 * .expr6
#     grad_b <- -(.expr1 * (.expr6 * (.expr3^(a4 - 1) * (a4 * (1/a4)))))
#     c(grad_a, grad_b)
# }

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

#' .error
#' @export
.error <- function(
    par, x, t,
    fun = c(doubleLog.elmore,
           doubleLog.beck,
           doubleLog.klos,
           doubleLog.gu,
           doubleLog.zhang,
           doubleAG),
    w, ylu, ...) {
    # FUN <- match.fun(fun)
    if (!all(is.finite(par))) return(9999)

    xpred <- fun(par, t = t)
    # If have no finite values, return 9999
    if (!all(is.finite(xpred))) return(9999) # for Klos fitting
    if (missing(w)) w <- rep(1, length(x))

    if (!missing(ylu)){
        # points out of ylu should be punished!
        w[xpred < ylu[1] | xpred > ylu[2]] <- 0
        # xpred   <- check_fit(xpred, ylu)
    }
    SSE  <- sum((x - xpred)^2 * w)
    RMSE <- sqrt(SSE/length(x))
    NSE  <- sum((x - xpred)^2 * w)/sum((x - mean(xpred))^2) 

    # 1. better handle low and high values simulation
    # xpred_2 <- sqrt(xpred_2)
    # x_2     <- sqrt(x_2)
    # xpred_2 <- log(xpred_2+1)
    # x_2     <- log(x_2+1)
    # xpred_2 <- 1/xpred          # inverse NSE
    # x_2     <- 1/x
    
    # xpred_2 <- xpred - mean(x)
    # x_2     <- x - mean(x)
    # NSE2 <- sum((x_2 - xpred_2)^2 * w)/sum((x_2 - mean(x_2))^2) #NSE
    
    # const <- ylu[2]
    # xpred_2 <- xpred - ylu[1]; xpred_2[xpred_2 < 0] <- const
    # x_2     <- x     - ylu[1]; x_2[x_2 < 0] <- const    
    return(RMSE)
}

# attach gradient and hessian analytical function to curve fitting functions
.dls <- lapply(c("doubleLog.beck", "doubleLog.elmore", "doubleLog.gu",
                 "doubleLog.klos", "doubleLog.zhang", "doubleAG"),
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
