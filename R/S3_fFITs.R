#' @name fFITs
#' @title S3 class of multiple fine curve fittings object.
#'
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
#' methods <- c("AG", "Beck", "Elmore", "Gu", "Zhang") # "Klos" too slow
#' fFITs <- curvefit(y, t, tout, methods)
#'
#' # plot
#' plot(fFITs)
NULL


#' plot function for phenofit class
#'
#' plot curve fitting VI, gradient (first order difference D1), hessian (D2),
#' curvature (k) and the change rate of curvature(der.k)
#'
#' @inheritParams PhenoExtract
#' @param x Fine curve fitting object \code{\link{fFITs}} returned by
#' \code{\link{curvefit}}.
#' @param ... ignored.
#'
#' @rdname fFITs
#' @export
plot.fFITs <- function(x, method, ...){

    if (missing(method)) {
        methods <- names(x$fFIT)
    } else {
        methods <- method
    }

    nmeth   <- length(methods)

    t    <- x$tout
    # Need to cal derivatives first
    for (i in seq_along(methods)){
        method <- methods[i]
        fFIT   <- x$fFIT[[i]]

        pred   <- last(fFIT$zs)

        pop    <- t[which.max(pred)]
        derivs <- curvature(fFIT)

        # plot for every method
        multi_p <- max(par()$mfrow) == 1
        op <- ifelse(multi_p, par(mfrow = c(2, 4), oma = c(0, 0, 1, 0)), par())

        # browser() # need to fix data error
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
        PhenoTrs(fFIT, IsPlot = TRUE, trs = 0.2)
        PhenoTrs(fFIT, IsPlot = TRUE, trs = 0.5)
        metrics <- PhenoGu(fFIT, IsPlot = TRUE)
        metrics <- PhenoKl(fFIT, IsPlot = TRUE)

        # show figure title
        op <- par(new = TRUE, mfrow = c(1, 1), oma = rep(0, 4), mar = rep(0, 4))
        plot(0, axes = F, type = "n", xaxt = "n", yaxt = "n") #
        text(1, 1, method, font = 2, cex = 1.3)
        par(op)
    }

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
}


#' GOF_fFITs
#'
#' Good-of-fitting (GOF) to fine curve fitting results.
#' 
#' @param fFITs \code{fFITs} object returned by \code{\link{curvefit}}.
#'
#' @return
#' \describe{
#'   \item{meth}{The name of fine curve fitting method}
#'   \item{RMSE}{Root Mean Square Error}
#'   \item{NSE}{Nash-Sutcliffe model efficiency coefficient}
#'   \item{R}{Pearson-Correlation}
#'   \item{pvalue}{pvalue of \code{R}}
#'   \item{n}{The number of observations}
#' }
#'
#' @references
#' [1]. https://en.wikipedia.org/wiki/Nash-Sutcliffe_model_efficiency_coefficient \cr
#' [2]. https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
#' 
#' @seealso \code{\link{curvefit}}
#' 
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
#' methods <- c("AG", "Beck", "Elmore", "Gu", "Zhang") # "Klos" too slow
#' fFITs <- curvefit(y, t, tout, methods)
#'
#' GOF_fFITs(fFITs)
#' @export
GOF_fFITs <- function(fFITs){
    methods <- names(fFITs$fFIT)
    nmeth   <- length(methods)

    t     <- fFITs$data$t
    tout  <- fFITs$tout
    ti    <- intersect(t, tout)
    I_org <- match(ti, t)
    I_sim <- match(ti, tout)

    Y_obs  <- fFITs$data$y[I_org]
    # Y_sim  <- map(x$fFIT, ~last(.$fits))

    # The following script assume that tout in every method is equal length.
    # calculate statistic for every meth
    info <- ldply(fFITs$fFIT, function(fFIT){
        Y_sim <- last(fFIT$zs)[I_sim]
        I <- which(!(is.na(Y_obs) | is.na(Y_sim)))

        # n_obs <- length(Y_obs)
        Y_sim2 <- Y_sim[I]
        Y_obs2 <- Y_obs[I]

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

        c(RMSE = RMSE, NSE = NSE, R = R, pvalue = pvalue, n = n)
    }, .id = "meth")
    info
}

# GOFs_fFITs <- function(fit){
# }
