#' Weight updating functions
#'
#' \itemize{
#'    \item \code{wSELF} weigth are not changed and return the original.
#'    \item \code{wTSM} weight updating method in TIMESAT.
#'    \item \code{wBisquare} Bisquare weight update method. wBisquare has been
#' modified to emphasis on upper envelope.
#'    \item \code{wChen} Chen et al., (2004) weight updating method.
#'    \item \code{wBeck} Beck et al., (2006) weigth updating method. wBeck need
#' sos and eos input. The function parameter is different from others. It is
#' still not finished.
#' }
#'
#' @inheritParams wHANTS
#' @param yfit Numeric vector curve fitting values.
#' @param iter iteration of curve fitting.
#' @param wfact weight adaptation factor (0-1), equal to the reciprocal of
#' 'Adaptation strength' in TIMESAT.
#' @param ... other parameters are ignored.
#'
#' @return wnew Numeric Vector, adjusted weights.
#'
#' @references
#' [1]. Per J\"onsson, P., Eklundh, L., 2004. TIMESAT - A program for analyzing
#'     time-series of satellite sensor data. Comput. Geosci. 30, 833-845.
#'     https://doi.org/10.1016/j.cageo.2004.05.006. \cr
#' [2]. https://au.mathworks.com/help/curvefit/smoothing-data.html#bq_6ys3-3 \cr
#' [3]. Garcia, D., 2010. Robust smoothing of gridded data in one and higher
#' dimensions with missing values. Computational statistics & data analysis,
#' 54(4), pp.1167-1178. \cr
#' [4]. Chen, J., J\"onsson, P., Tamura, M., Gu, Z., Matsushita, B., Eklundh, L.,
#'      2004. A simple method for reconstructing a high-quality NDVI time-series
#'      data set based on the Savitzky-Golay filter. Remote Sens. Environ. 91,
#'      332-344. https://doi.org/10.1016/j.rse.2004.03.014. \cr
#' [5]. Beck, P.S.A., Atzberger, C., Hogda, K.A., Johansen, B., Skidmore, A.K.,
#'      2006. Improved monitoring of vegetation dynamics at very high latitudes:
#'      A new method using MODIS NDVI. Remote Sens. Environ.
#'      https://doi.org/10.1016/j.rse.2005.10.021 \cr
#' [6]. https://github.com/kongdd/phenopix/blob/master/R/FitDoubleLogBeck.R
#'
#' @rdname wFUN
#' @export
wSELF <- function(y, yfit, w, ...){w}

#' @author
#' wTSM is implemented by Per J\"onsson, Malm\"o University, Sweden
#' \email{per.jonsson@ts.mah.se} and Lars Eklundh, Lund University, Sweden
#' \email{lars.eklundh@nateko.lu.se}. And Translated into Rcpp by Dongdong Kong,
#' 01 May 2018.
#'
#' @rdname wFUN
#' @export
wTSM <- function(y, yfit, w, iter = 2, nptperyear, wfact = 0.5, ...){
    wTSM_cpp(y, yfit, w, iter, nptperyear, wfact)
}

#' @rdname wFUN
#' @export
wBisquare <- function(y, yfit, w, ..., wmin = 0.2){
    
    if (missing(w)) w  <- rep(1, length(y))
    wnew <- w

    # Update to avoid decrease the weights ungrowing season points too much
    # This idea is also occur in wTSM and wBeck;
    # 2018-07-25
    ylu <- range(yfit, na.rm = T)
    A = diff(ylu)
    re     <- yfit - y
    re_abs <- abs(re)

    sc     <- 6 * median(re_abs, na.rm = T)

    # only decrease the weights of growing season `& yfit > 1/3*A`.
    I_pos        <- which(re > 0 & re < sc & (yfit > 0.3*A + ylu[1]) )
    wnew[I_pos]  <- (1 - (re_abs[I_pos]/sc)^2)^2 * w[I_pos]
    # have a problem, in this way, original weights will be ignored.

    # I_zero       <- which(re >= sc) # update 20180723
    I_zero       <- which(abs(re) >= sc) # update 20180924, positive bias outlier
    
    wnew[I_zero] <- wmin
    wnew[wnew < wmin] <- wmin

    # constrain growing VI: enlarge the positive bias values, reduce the weights
    #   of negative bias values, as TIMESAT
    # diff = 2 * (yfit(i) - y(i)) / yfitstd;
    # w <- wfact * wfit(i) * exp( - ydiff ^ 2);
    return(wnew)
}

#' @rdname wFUN
#' @export
wChen <- function(y, yfit, w, ..., wmin = 0.2){
    if (missing(w)) w  <- rep(1, length(y))
    wnew <- w

    re     <- yfit - y
    re_abs <- abs(re)

    d_max  <- max(re_abs, na.rm = T) #6 * median(re, na.rm = T)

    I_pos  <- re > 0
    wnew[ I_pos]      <- (1 - re_abs[I_pos] / d_max) * w[I_pos]
    wnew[wnew < wmin] <- wmin
    return(wnew)
}

# ' #@export
#' @rdname wFUN
wBeck <- function(y, yfit, w, ...){
    # get optimized parameters
    t   <- seq_along(y)
    sos <- opt$par[3]
    eos <- opt$par[5]

    m  <- lm(c(0, 100) ~ c(sos, eos))
    tr <- coef(m)[2] * t + coef(m)[1]
    tr[tr < 0] <- 0
    tr[tr > 100] <- 100

    # estimate weights
    res <- xpred - x
    weights <- 1/((tr * res + 1)^2)
    weights[res > 0 & res <= 0.01] <- 1
    weights[res < 0] <- 4
}
