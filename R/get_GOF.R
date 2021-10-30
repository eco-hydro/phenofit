#' get_GOF
#'
#' Goodness-of-fitting (GOF) of fine curve fitting results.
#'
#' @inheritParams get_param
#' 
#' @return
#' - `meth`: The name of fine curve fitting method
#' - `RMSE`: Root Mean Square Error
#' - `NSE` : Nash-Sutcliffe model efficiency coefficient
#' - `R`   : Pearson-Correlation
#' - `R2`  : determined coefficient
#' - `pvalue`: pvalue of `R`
#' - `n`   : The number of observations
#'
#' @references
#' 1. https://en.wikipedia.org/wiki/Nash-Sutcliffe_model_efficiency_coefficient \cr
#' 2. https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
#'
#' @seealso [curvefit()]
#'
#' @example inst/examples/ex-get_fitting_param_GOF.R
#' @export
#'
#' @example inst/examples/ex-get_fitting_param_GOF.R
#' @export
get_GOF <- function(x) UseMethod("get_GOF", x)

#' @rdname get_GOF
#' @export
get_GOF.list <- function(x){
    map_df(x, get_GOF.fFITs, .id = "flag") %>% data.table()
}

#' @rdname get_GOF
#' @export
get_GOF.fFITs <- function(x){
    models = x$model
    methods <- names(models)
    nmeth   <- length(methods)

    t     <- x$data$t
    tout  <- x$tout
    ti    <- intersect(t, tout)
    I_org <- match(ti, t)
    I_sim <- match(ti, tout)

    Y_obs  <- x$data$y[I_org]
    # Y_sim  <- map(x$fFIT, ~last(.$fits))

    # The following script assume that tout in every method is equal length.
    # calculate statistic for every meth
    info <- map_df(models, function(fFIT){
        Y_sim <- last(fFIT$zs)[I_sim]
        I <- which(!(is.na(Y_obs) | is.na(Y_sim)))

        Y_sim2 <- Y_sim[I]
        Y_obs2 <- Y_obs[I]

        indexes = c("R2", "NSE", "R", "RMSE", "pvalue", "n_sim")
        GOF(Y_obs2, Y_sim2)[indexes] %>% as.list() %>% as.data.table()
    }, .id = "meth")
    info
}
