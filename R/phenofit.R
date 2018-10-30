#' statistic
#'
#' statistics to evaluate curve fitting performance
#' @param fit The phenofit object
#' @export
statistic.phenofit <- function(fit){
    t     <- fit$data$t
    tout  <- fit$tout
    ti    <- intersect(t, tout)
    I_org <- match(ti, t)
    I_sim <- match(ti, tout)

    Y_obs  <- fit$data$y[I_org]
    Y_sim  <- last(fit$fits)[I_sim]

    I <- which(!(is.na(Y_obs) | is.na(Y_sim)))
    # n_obs <- length(Y_obs)

    Y_sim <- Y_sim[I]
    Y_obs <- Y_obs[I]

    R      <- NA
    pvalue <- NA
    n      <- length(I)

    if (n < 2) return(c(RMSE = NA, NSE = NA, R = R, pvalue = pvalue, n = n))    

    tryCatch({
        cor.obj <- cor.test(Y_obs, Y_sim, use = "complete.obs")
        R       <- cor.obj$estimate[[1]]
        pvalue  <- cor.obj$p.value
    }, error = function(e){
        message(sprintf('[statistic] %s', e$message))
    })

    RMSE <- sqrt(sum((Y_obs - Y_sim)^2)/n)
    NSE  <- 1 - sum((Y_sim - Y_obs)^2) / sum((Y_obs - mean(Y_obs))^2)
    return(c(RMSE = RMSE, NSE = NSE, R = R, pvalue = pvalue, n = n))
}

#' Gradient (grad) and Hessian (hess) based on \code{numDeriv} package
#' 
#' @param fit Curve fitting result returned by \code{curvefit}.
#' @param tout A vector of time steps at which the function can be predicted
#' 
#' @rdname derivative
#' @export
grad.phenofit <- function(fit, tout){
    FUN <- get(fit$fun, mode = 'function')
    grad(FUN, tout, par= fit$par)
}

#' @rdname derivative
#' @export
hess.phenofit <- function(fit, tout){
    FUN <- get(fit$fun, mode = 'function')
    grad(function(t) grad(FUN, t, par= fit$par), tout)
}

#' @export
D1 <- function(fit, ...) UseMethod('D1', fit)

#' @export
D2 <- function(fit, ...) UseMethod('D2', fit)

#' @export
curvature <- function(fit, ...) UseMethod('curvature', fit)

#' @export
print.phenofit <- function(x, ...){
    FUN <- get(x$fun, mode = 'function')
    cat(sprintf(" formula:\t%s\n", attr(FUN, "formula") ))
    cat("pars:\n")
    print(x$par)
}

#' @rdname derivative
#' @param numDeriv If true, \code{numDeriv} package \code{grad} and \code{hess} 
#' will be used; if false, \code{D1} and \code{D2} will be used.
#' @param smspline If smspline = TRUE or attr(fit$fun, 'gradient') == null,
#'                 it will use smooth.spline to get der1 and der2.
#' @param ... Other parameters are ignored.
#' 
#' @export
D1.phenofit <- function(fit, numDeriv = FALSE, smspline = FALSE, ...){
    pred <- last(fit$fits)                 #zoo obj
    t    <- fit$tout
    par  <- fit$par

    FUN <- get(fit$fun, mode = 'function')
    D1  <- attr(FUN, 'gradient')# first order derivative, D1 was 6 times faster
                                     # than grad, and 20 times faster then diff
    # the derivate of curve fitting time-series
    if (numDeriv){
        der1 <- grad.phenofit(fit, t)
    }else{
        if (is.null(D1) || smspline) {
            # numerical solution
            spline.eq <- smooth.spline(pred, df = length(pred))
            der1      <- predict(spline.eq, d = 1)$y
            # der1 <- diff(pred)/diff(t)
        } else {
            der1 <- D1(par, t)[, 1]
        }
    }
    der1[is.infinite(der1)] <- NA
    #rule = 2, means y range out of xlim also could get a approximate value
    if (any(is.na(der1))) der1 %<>% na.approx(rule = 2)
    return(der1)
}

#' @rdname derivative
#' @export
D2.phenofit <- function(fit, numDeriv = FALSE, smspline = FALSE, ...){
    pred <- last(fit$fits)
    t    <- fit$tout
    par  <- fit$par

    FUN  <- get(fit$fun, mode = 'function')
    D1   <- attr(FUN, 'gradient') # first order derivative, D1 was 6 times faster
                                  # than grad, and 20 times faster then diff
    D2   <- attr(FUN, 'hessian')  # second order derivative

    if (numDeriv){
        der1 <- grad.phenofit(fit, t)
        der2 <- hess.phenofit(fit, t)
    }else{
        if (is.null(D2)  || smspline) {
            spline.eq <- smooth.spline(pred, df = length(pred))
            der1 <- predict(spline.eq, d = 1)$y
            der2 <- predict(spline.eq, d = 2)$y
        } else{
            der1 <- D1(par, t)[, 1]
            der2 <- D2(par, t)[, 1, 1]
        }
    }
    ## in case for NA values
    der1[is.infinite(der1)] <- NA
    der2[is.infinite(der2)] <- NA
    if (any(is.na(der1))) der1 %<>% na.approx(rule = 2)
    if (any(is.na(der2))) der2 %<>% na.approx(rule = 2)

    return(list(der1 = der1, der2 = der2))
}

#' @rdname derivative
#' @export
curvature.phenofit <- function(fit, numDeriv = FALSE, smspline = FALSE, ...){
    derivs <- D2.phenofit(fit, numDeriv = numDeriv, smspline = smspline)
    k      <- derivs$der2 / (1 + derivs$der1 ^ 2) ^ (3 / 2)
    return(list(k = k, der1 = derivs$der1, der2 = derivs$der2))
}

#' plot function for phenofit class
#'
#' plot curve fitting VI, gradient (first order difference D1), hessian (D2),
#' curvature (k) and the change rate of curvature(der.k)
#'
#' @inheritParams grad.phenofit
#' @export
plot.phenofit <- function(x, ...){
    name <- deparse(substitute(x))

    pred <- last(x$fits)                 #zoo obj
    t    <- x$tout

    pop    <- t[which.max(pred)]
    derivs <- curvature(x)

    # e <- environment(D1)
    # der1_diff <- c(NA, diff(pred))
    # microbenchmark::microbenchmark(
    #     gradient <- numDeriv::grad(x$fun,  t, par = x$par),
    #     der1     <- D2(x),
    #     diff1 <- c(NA, NA, diff(diff(values)))
    # )
    # summary(der1 - gradient)
    # summary(c(NA, diff1) - gradient)
    # diff <- der1 - der1_diff
    # der <- data.frame(der1_diff, der1 = der1, diff = diff)

    multi_p <- max(par()$mfrow) == 1
    op <- ifelse(multi_p, par(mfrow = c(2, 4), oma = c(0, 0, 1, 0)), par())
    plot(x$data$t, x$data$y,
        type= "b", pch = 20, cex = 1.3, col = "grey",
        main = "curve fitting VI", xlab = "Index", ylab = "VI")
    lines(t, pred); grid()
    abline(v = pop, col ="green")

    maxd_der1 <- t[which.max(derivs$der1)]
    mind_der1 <- t[which.min(derivs$der1)]

    plot(t, derivs$der1, main = "D1"); grid()
    abline(v = maxd_der1, col ="blue")
    abline(v = mind_der1, col ="red")
    abline(v = pop, col ="green")

    plot(t, derivs$der2, main = "D2"); grid()
    plot(t, derivs$k, main = "k")    ; grid()
    abline(v = maxd_der1, col ="blue")
    abline(v = mind_der1, col ="red")
    abline(v = pop, col ="green")
    abline(v = pop + 20, col ="green", lty = 2)
    # plot(diff(der1_diff), main = "diff2")

    # k <- derivs$k
    PhenoTrs(x, IsPlot = TRUE, trs = 0.2)
    PhenoTrs(x, IsPlot = TRUE, trs = 0.5)
    metrics <- PhenoGu(x, IsPlot = TRUE)
    metrics <- PhenoKl(x, IsPlot = TRUE)

    # show figure title
    op <- par(new = TRUE, mfrow = c(1, 1), oma = rep(0, 4), mar = rep(0, 4))
    plot(0, axes = F, type = "n", xaxt = "n", yaxt = "n") #
    text(1, 1, name, font = 2, cex = 1.3)
    par(op)
}
