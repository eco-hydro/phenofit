#' get_GOF
#'
#' Goodness-of-fitting (GOF) of fine curve fitting results.
#' 
#' @param fit Object returned by `curvefits`.
#' @param fFITs `fFITs` object returned by [curvefit()].
#'
#' @return
#' \describe{
#'   \item{meth}{The name of fine curve fitting method}
#'   \item{RMSE}{Root Mean Square Error}
#'   \item{NSE}{Nash-Sutcliffe model efficiency coefficient}
#'   \item{R}{Pearson-Correlation}
#'   \item{pvalue}{pvalue of `R`}
#'   \item{n}{The number of observations}
#' }
#'
#' @references
#' 1. https://en.wikipedia.org/wiki/Nash-Sutcliffe_model_efficiency_coefficient \cr
#' 2. https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
#' 
#' @seealso [curvefit()]
#' 
#' @example inst/examples/ex-get_fitting_param_GOF.R
#' @export
get_GOF <- function(fit){
    ldply(fit, get_GOF.fFITs, .id = "flag") %>% data.table()
}

#' @rdname get_GOF
#' @export
get_GOF.fFITs <- function(fFITs){
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
            cor.obj <- cor.test(Y_obs2, Y_sim2, use = "complete.obs")
            R       <- cor.obj$estimate[[1]]
            pvalue  <- cor.obj$p.value
        }, error = function(e){
            message(sprintf('[statistic] %s', e$message))
        })

        RMSE <- sqrt(sum((Y_obs2 - Y_sim2)^2)/n)
        NSE  <- 1 - sum((Y_obs2 - Y_sim2)^2) / sum((Y_obs2 - mean(Y_obs2))^2)

        c(RMSE = RMSE, NSE = NSE, R = R, pvalue = pvalue, n = n)
    }, .id = "meth")
    info
}
