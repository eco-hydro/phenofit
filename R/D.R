# ' @rdname derivative
# ' @export
hess.fFIT <- function(fit, tout = NULL){
    FUN <- get(fit$fun, mode = 'function')
    grad(function(t) grad(FUN, t, par= fit$par), tout)
}


# ' Get gradient and hessian by \code{grad} in \code{numDeriv} package.
# '
# ' @param fit A curve fitting object returned by \code{curvefit}.
# ' @param tout A vector of time steps at which the function can be predicted.
# '
# ' @examples
# ' FUN <- doubleLog.Beck
# '
# ' @rdname derivative
# ' @export
grad.fFIT <- function(fit, tout = NULL){
    FUN <- get(fit$fun, mode = 'function')
    grad(FUN, tout, par = fit$par)

}

#' @title D
#' @name D
#'
#' @description Get derivative of `phenofit` object.
#' `D1` first order derivative, `D2` second order derivative, n
#' `curvature` curvature.
#'
#' @details If `fit$fun` has no gradient function or `smoothed.spline = TRUE`,
#' time-series smoothed by spline first, and get derivatives at last.
#' If `fit$fun` exists and `analytical = TRUE`, `smoothed.spline`
#' will be ignored.
#'
#' @param fit A curve fitting object returned by `curvefit`, with the object of:
#' - `par`: parameters of curve fitting function
#' - `fun`: curve fitting function name, e.g., "doubleLog_AG"
#' - `zs`: predicted values, vector or data.frame
#'
#' @param analytical If true, `numDeriv` package `grad` and `hess`
#' will be used; if false, `D1` and `D2` will be used.
#' @param smoothed.spline Whether apply `smooth.spline` first?
#' @param ... Other parameters will be ignored.
#'
#' @keywords internal
#' @return
#' \itemize{
#' \item der1 First order derivative
#' \item der2 Second order derivative
#' \item k    Curvature
#' }
#'
#' @example R/examples/ex-D1.R
#' @rdname D
NULL

#' @rdname D
#' @export
D1 <- function(fit, t = NULL, analytical = FALSE, smoothed.spline = FALSE, ...) UseMethod('D1', fit)

#' @rdname D
#' @export
D2 <- function(fit, t = NULL, analytical = FALSE, smoothed.spline = FALSE, ...) UseMethod("D2", fit)

#' @keywords internal
#' @rdname D
#' @export
D1.fFIT <- function(fit, t = NULL, analytical = FALSE, smoothed.spline = FALSE, ...){
    if (is.null(t)) t = fit$tout
    pred <- last2(fit$zs)
    # t    <- fit$tout
    par  <- fit$par

    FUN <- get(fit$fun, mode = 'function')
    D1  <- attr(FUN, 'gradient')# first order derivative, D1 was 6 times faster
                                # than grad, and 20 times faster then diff
    # the derivate of curve fitting time-series

    if (analytical && !is.null(D1)) smoothed.spline <- FALSE

    if (is.null(D1) || smoothed.spline) {
        # 1. Numerical solution
        spline.eq <- smooth.spline(pred, df = length(pred)*0.1)
        der1      <- predict(spline.eq, d = 1)$y
        # der1 <- diff(pred)/diff(t)
    } else if (analytical){
        # real analytical
        der1 <- D1(par, t)[, 1] # the default option
    } else {
        # numerical approximation by package `numDeriv`
        der1 <- grad.fFIT(fit, t)
    }

    der1[is.infinite(der1)] <- NA
    #rule = 2, means y range out of xlim also could get a approximate value
    if (any(is.na(der1))) der1 %<>% na.approx(rule = 2)
    return(der1)
}

#' @keywords internal
#' @rdname D
#' @export
D2.fFIT <- function(fit, t = NULL, analytical = FALSE, smoothed.spline = FALSE, ...){
    if (is.null(t)) t = fit$tout
    pred <- last2(fit$zs)
    # t    <- fit$tout
    par  <- fit$par

    FUN  <- get(fit$fun, mode = 'function')
    D1   <- attr(FUN, 'gradient') # first order derivative, D1 was 6 times faster
                                  # than grad, and 20 times faster then diff
    D2   <- attr(FUN, 'hessian')  # second order derivative

    if (is.null(D1) || smoothed.spline) {
        # 1. Numerical solution
        spline.eq <- smooth.spline(pred, df = length(pred))
        der1      <- predict(spline.eq, d = 1)$y
        der2      <- predict(spline.eq, d = 2)$y
    } else if (analytical){
        # real analytical
        der1 <- D1(par, t)[, 1]
        der2 <- D2(par, t)[, 1, 1]
    } else {
        # numerical approximation
        der1 <- grad.fFIT(fit, t)
        der2 <- hess.fFIT(fit, t)
    }

    ## in case for NA values
    der1[is.infinite(der1)] <- NA
    der2[is.infinite(der2)] <- NA
    if (any(is.na(der1))) der1 %<>% na.approx(rule = 2)
    if (any(is.na(der2))) der2 %<>% na.approx(rule = 2)

    return(list(der1 = der1, der2 = der2))
}

#' @rdname D
#' @export
curvature <- function(fit, t = NULL, analytical = FALSE, smoothed.spline = FALSE, ...) UseMethod('curvature', fit)

#' @keywords internal
#' @rdname D
#' @export
curvature.fFIT <- function(fit, t = NULL, analytical = FALSE, smoothed.spline = FALSE, ...){
    derivs <- D2.fFIT(fit, t, analytical, smoothed.spline)
    k      <- derivs$der2 / (1 + derivs$der1 ^ 2) ^ (3 / 2)
    c(derivs, list(k = k))
}

#' @importFrom zoo na.spline
rm_spike <- function(y, times = 3, halfwin = 1, maxgap = 4) {
    # 强化除钉值模块, 20191127
    std <- sd(y, na.rm = TRUE)
    # ymov <- cbind(y[c(1, 1:(n - 2), n-1)], y[c(2, 3:n, n)]) %>% rowMeans(na.rm = TRUE)
    # # ymov2 <- movmean(y, 1)
    # halfwin <- ceiling(nptperyear/36) # about 10-days
    ymov2 <- movmean(y, halfwin = halfwin)
    # which(abs(y - ymean) > std) & w <= w_critical
    #  | abs(y - ymov2) > 2*std
    I_spike <- which(abs(y - ymov2) > times * std) # 95.44% interval, `(1- 2*pnorm(-2))*100`
    # print(I_spike)
    if (length(I_spike) > 0) {
        y[I_spike] <- NA # missval
        y = na.spline(y, maxgap = maxgap, na.rm = FALSE)
    }
    y
}
