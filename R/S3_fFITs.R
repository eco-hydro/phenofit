#' @name fFIT
#' @title S3 class of fine curve fitting object.
#'
#' @description
#' `fFIT` is returned by [optim_pheno()].
#'
#' @format
#' * `tout`: Corresponding doy of prediction
#' * `zs`: curve fitting values of every iteration
#' * `ws`: weight of every iteration
#' * `par`: Optimized parameter of fine curve fitting method
#' * `fun`: The name of fine curve fitting function.
#'
#' @keywords internal
NULL

#' @export
print.fFIT <- function(x, ...) {
    FUN <- get(x$fun, mode = "function")
    cat(sprintf(" formula:\t%s\n", attr(FUN, "formula")))
    cat("pars:\n")
    print(x$par)
}

#' @name fFITs
#' @title S3 class of multiple fine curve fittings object.
#'
#' @examples
#' library(phenofit)
#' # simulate vegetation time-series
#' fFUN = doubleLog.Beck
#' par  = c( mn  = 0.1 , mx  = 0.7 , sos = 50 , rsp = 0.1 , eos = 250, rau = 0.1)
#' 
#' t    <- seq(1, 365, 8)
#' tout <- seq(1, 365, 1)
#' y <- fFUN(par, t)
#'
#' methods <- c("AG", "Beck", "Elmore", "Gu", "Zhang") # "Klos" too slow
#' fit <- curvefit(y, t, tout, methods)
#'
#' # plot
#' plot(fit)
#' @keywords internal
NULL

#' plot function for phenofit class
#'
#' plot curve fitting VI, gradient (first order difference D1), hessian (D2),
#' curvature (k) and the change rate of curvature(der.k)
#'
#' @inheritParams get_pheno
#' @param x Fine curve fitting object [fFITs()] returned by [curvefit()].
#' @param show.pheno whether to plot phenological metrics.
#' @param ... other parameters to [curvature()].
#' 
#' @seealso [curvature()]
#' @rdname fFITs
#' @keywords internal
#' @export
plot.fFITs <- function(x, method, show.pheno = TRUE, ...){
    if (missing(method)) {
        methods <- names(x$model)
    } else {
        methods <- method
    }
    nmeth <- length(methods)

    t <- x$tout
    # Need to cal derivatives first
    for (i in seq_along(methods)){
        method <- methods[i]
        fFIT <- x$model[[method]]
        pred   <- last2(fFIT$zs)

        pos    <- t[which.max(pred)] # if type = 1, should be reversed
        derivs <- curvature(fFIT, t, ...)

        # plot for every method
        multi_p <- max(par()$mfrow) == 1
        
        mfrow = if (show.pheno) c(2, 4) else c(2, 2)
        op = par(mfrow = mfrow, oma = c(0, 1, 2, 0), mar = c(3, 3, 1, 1), mgp = c(2, 0.6, 0))
        # op <- ifelse(multi_p, , par())

        plot(x$data$t, x$data$y,
            type= "b", pch = 20, cex = 1.3, col = "grey",
            ylab = "Fine fitting result")
        lines(t, pred); grid()
        abline(v = pos, col ="green")

        max_der1 <- t[which.max(derivs$der1)]
        min_der1 <- t[which.min(derivs$der1)]

        plot(t, derivs$der1, type = "l", ylab = "D1"); grid()
        abline(v = max_der1, col ="blue")
        abline(v = min_der1, col ="red")
        abline(v = pos, col ="green")
        abline(h = 0, col = "grey", lty = 2)

        plot(t, derivs$der2, type = "l", ylab = "D2"); grid()
        abline(h = 0, col = "grey", lty = 2)
        
        plot(t, derivs$k, type = "l", ylab = "K'"); grid()
        abline(v = max_der1, col ="blue")
        abline(v = min_der1, col ="red")
        abline(v = pos, col ="green")
        abline(v = pos + 20, col ="green", lty = 2)
        abline(h = 0, col = "grey", lty = 2)
        # plot(diff(der1_diff), main = "diff2")

        if (show.pheno) {
            # k <- derivs$k
            PhenoTrs(fFIT, IsPlot = TRUE, trs = 0.2)
            PhenoTrs(fFIT, IsPlot = TRUE, trs = 0.5)
            metrics <- PhenoGu(fFIT, IsPlot = TRUE)
            metrics <- PhenoKl(fFIT, IsPlot = TRUE)
        }
        par(op)
        # show figure title
        op = par(new = TRUE, mfrow = c(1, 1), oma = rep(0, 4), mar = rep(0, 4))
        plot(0, axes = F, type = "n", xaxt = "n", yaxt = "n") #
        text(1, 1, method, font = 2, cex = 1.3)
        par(op)
    }
}

#' @inheritParams get_GOF
#' @rdname fFITs
#' @keywords internal
#' @export
plot.fFIT <- function(x, data, ...) {
    t <- x$tout
    pred <- last2(x$zs)

    plot(data$t, data$y,
        type = "b", pch = 20, cex = 1.3, col = "grey",
        ylab = "Fine fitting result"
    )
    lines(t, pred)
    grid()

    legend('topright', c('Original', 'Fine fitting'), 
        col = c("grey", "black"),
        lty = c(0, 1), pch = c(16, NA))
}

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
