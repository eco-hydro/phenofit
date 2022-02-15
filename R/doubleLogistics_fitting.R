# #' using cubic spline function to avoid the difficult in setting parameter
# #' lambda in smooth.spline
# #'
# #' cubic spline is inappropriate for daily inputs. Its smooth is not enough.
# #' On the contrary, smooth.spline with a low freedom can smooth well.
# #'
# #' @export
# splinefit <- function(y, t = index(y), tout = t, plot = FALSE, df.factor = 0.06, ...){
#     # xpred.out <- spline(t, y, xout = tout)$y %>% zoo(., tout)
#     n <- length(y)
#     # if n < 40, means y was satellite VI
#     # if n > 40, means daily data
#     df.factor <- ifelse (n <= 46, 1/3, df.factor)
#     freedom   <- pmax(df.factor * n, 15)
#     fit       <- smooth.spline(t, y, df = freedom)
#     xpred.out <- predict(fit, tout)$y %>% zoo(., tout)
#     structure(list(data = list(y = y, t = t),
#         pred = xpred.out, par = NULL, fun = NULL), class = "phenofit")
# }

#' @name FitDL
#' @title Fine fitting
#' 
#' @description Fine curve fitting function is used to fit vegetation 
#' time-series in every growing season.
#'
#' @param y input vegetation index time-series.
#' @param t the corresponding doy(day of year) of y.
#' @param tout the time of output curve fitting time-series.
#' @param method method passed to `optimx` or `optim` function.
#' @param w weights
#' @param ... other paraters passed to [optim_pheno()].
#' @inheritParams init_param
#' 
#' @return
#' - `tout`: The time of output curve fitting time-series.
#' - `zs`  : Smoothed vegetation time-series of every iteration.
#' - `ws`  : Weights of every iteration.
#' - `par` : Final optimized parameter of fine fitting.
#' - `fun` : The name of fine fitting.
#' 
#' @references
#' 1. Beck, P.S.A., Atzberger, C., Hogda, K.A., Johansen, B., Skidmore, A.K.,
#'      2006. Improved monitoring of vegetation dynamics at very high latitudes:
#'      A new method using MODIS NDVI. Remote Sens. Environ.
#'      https://doi.org/10.1016/j.rse.2005.10.021.
#' 2. Elmore, A.J., Guinn, S.M., Minsley, B.J., Richardson, A.D., 2012.
#'      Landscape controls on the timing of spring, autumn, and growing season
#'      length in mid-Atlantic forests. Glob. Chang. Biol. 18, 656-674.
#'      https://doi.org/10.1111/j.1365-2486.2011.02521.x. \cr
#' 
#' 3. Gu, L., Post, W.M., Baldocchi, D.D., Black, TRUE.A., Suyker, A.E., Verma,
#'      S.B., Vesala, TRUE., Wofsy, S.C., 2009. Characterizing the Seasonal Dynamics
#'      of Plant Community Photosynthesis Across a Range of Vegetation Types,
#'      in: Noormets, A. (Ed.), Phenology of Ecosystem Processes: Applications
#'      in Global Change Research. Springer New York, New York, NY, pp. 35-58.
#'      https://doi.org/10.1007/978-1-4419-0026-5_2. \cr
#'
#' 4. https://github.com/cran/phenopix/blob/master/R/FitDoubleLogGu.R
#' @example R/examples/ex-FitDL.R
NULL

#' @rdname FitDL
#' @export
FitDL.Zhang <- function(y, t = index(y), tout = t, 
    method = 'nlm', w, type = 1L, ...)
{
    if (missing(w)) w <- rep(1, length(y))
    e <- init_param(y, t, w, type = type)
    
    sFUN    <- "doubleLog.Zhang"
    p = init_Zhang(e, type = type)
    optim_pheno(p$prior, sFUN, y, t, tout, method, w, lower = p$lower, upper = p$upper, ...)
}

#' @rdname FitDL
#' @export
FitDL.AG <- function(y, t = index(y), tout = t, 
    method = 'nlminb', w, type = 1L, ...)
{
    if (missing(w)) w <- rep(1, length(y))
    e <- init_param(y, t, w, type = type)
    # print(ls.str(envir = e))
    sFUN <- "doubleLog.AG"
    p = init_AG(e, type = type)
    optim_pheno(p$prior, sFUN, y, t, tout, method, w, lower = p$lower, upper = p$upper, ...)
}

# background value 非对称高斯分布
#' @rdname FitDL
#' @export
FitDL.AG2 <- function(y, t = index(y), tout = t, 
    method = 'nlminb', w, type = 1L, ...)
{
    if (missing(w)) w <- rep(1, length(y))
    e <- init_param(y, t, w, type = type)
    # print(ls.str(envir = e))
    sFUN <- "doubleLog.AG2"
    p = init_AG2(e, type = type)
    optim_pheno(p$prior, sFUN, y, t, tout, method, w, lower = p$lower, upper = p$upper, ...)
}

#' @rdname FitDL
#' @export
FitDL.Beck <- function(y, t = index(y), tout = t, 
    method = 'nlminb', w, type = 1L, ...)
{
    if (missing(w)) w <- rep(1, length(y))
    e <- init_param(y, t, w, type = type)

    sFUN   <- "doubleLog.Beck"
    p = init_Beck(e, type = type)
    optim_pheno(p$prior, sFUN, y, t, tout, method, w, lower = p$lower, upper = p$upper, ...)
}

#' @rdname FitDL
#' @export
FitDL.Elmore <- function(y, t = index(y), tout = t, 
    method = 'nlminb', w, type = 1L, ...) 
{
    e <- init_param(y, t, w, type = type)

    # doy_q  <- quantile(t, c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
    sFUN   <- "doubleLog.Elmore"
    
    # TODO: remove bad one 
    p = init_Elmore(e, type = type)
    optim_pheno(p$prior, sFUN, y, t, tout, method, w, lower = p$lower, upper = p$upper, ...)
}

# c(mn, mx - mn, doy[2], half*0.1, doy[4], half*0.1, 0.002),
# c(mn, mx - mn, doy[2], half*0.2, doy[5], half*0.2, 0.002),
# c(mn, mx - mn, doy[1], half*0.5, doy[4], half*0.5, 0.05),
# c(mn, mx - mn, doy[1], half*0.8, doy[5], half*0.8, 0.1))

#' @rdname FitDL
#' @export
FitDL.Gu <- function(y, t = index(y), tout = t, 
    method = "nlminb", w, type = 1L, ...)
{
    if (missing(w)) w <- rep(1, length(y))
    e <- init_param(y, t, w, type = type)

    sFUN  <- "doubleLog.Gu"
    p = init_Gu(e, type = type)
    optim_pheno(p$prior, sFUN, y, t, tout, method, w, lower = p$lower, upper = p$upper, ...)
}

#' @rdname FitDL
#' @export
FitDL.Klos <- function(y, t = index(y), tout = t, 
    method = 'BFGS', w, type = 1L, ...)
{
    if (missing(w)) w <- rep(1, length(y))
    e <- init_param(y, t, w, type = type)

    sFUN <- "doubleLog.Klos"
    p = init_Klos(e, type = type)
    optim_pheno(p$prior, sFUN, y, t, tout, method, w, lower = p$lower, upper = p$upper, ...)
}
