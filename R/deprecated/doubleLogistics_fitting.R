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
#' @param w weights
#' 
#' @param tout the time of output curve fitting time-series.
#' @param method method passed to `optimx` or `optim` function.
#' @param ... other paraters passed to [optim_pheno()].
#'

NULL
# c(mn, mx - mn, doy[2], half*0.1, doy[4], half*0.1, 0.002),
# c(mn, mx - mn, doy[2], half*0.2, doy[5], half*0.2, 0.002),
# c(mn, mx - mn, doy[1], half*0.5, doy[4], half*0.5, 0.05),
# c(mn, mx - mn, doy[1], half*0.8, doy[5], half*0.8, 0.1))
