#' wSELF
#' Weigths are unchanged and return the original.
wSELF <- function(y, z, w, ...){w}

#' wTSM
#'
#' Weight updating method in TIMESAT
#'
#' @author Per J\"onsson, Malm\"o University, Sweden \email{per.jonsson@ts.mah.se} \cr
#'     Lars Eklundh, Lund University, Sweden \email{lars.eklundh@nateko.lu.se} \cr
#'     Translate by Dongdong Kong, 01 May 2018.
#'
#' @references
#' [1]. Per J\"onsson, P., Eklundh, L., 2004. TIMESAT - A program for analyzing
#'     time-series of satellite sensor data. Comput. Geosci. 30, 833-845.
#'     https://doi.org/10.1016/j.cageo.2004.05.006
wTSM <- function(y, yfit, w, iters = 2, nptperyear, wfact = 0.5, ...){
    wTSM_cpp(y, yfit, w, iters, nptperyear, wfact)
}

#' Bisquare weight update method
#' 
#' Bisquare has been modified to emphasis on upper envelope.
#'
#' @references
#' [1]. https://au.mathworks.com/help/curvefit/smoothing-data.html#bq_6ys3-3 \cr
#' [2]. Garcia, D., 2010. Robust smoothing of gridded data in one and higher
#' dimensions with missing values. Computational statistics & data analysis,
#' 54(4), pp.1167-1178.
#' @export
wBisquare <- function(y, yfit, w, ...){
    if (missing(w)) w  <- rep(1, length(y))
    wnew <- w

    re     <- yfit - y
    re_abs <- abs(re)

    sc <- 6 * median(re_abs, na.rm = T)

    I_pos        <- which(re > 0 & re < sc)
    wnew[I_pos]  <- (1 - (re_abs[I_pos]/sc)^2)^2 * w[I_pos]

    I_zero       <- which(re_abs >= sc)
    wnew[I_zero] <- 0
    # constrain growing VI: enlarge the positive bias values, reduce the weights
    #   of negative bias values, as TIMESAT
    # diff = 2 * (yfit(i) - y(i)) / yfitstd;
    # w <- wfact * wfit(i) * exp( - ydiff ^ 2);
    return(wnew)
}

#' 
#' Chen et al., (2004) weights updating method
#' 
#' @references
#' [1]. Chen, J., Jönsson, P., Tamura, M., Gu, Z., Matsushita, B., Eklundh, L., 
#'      2004. A simple method for reconstructing a high-quality NDVI time-series 
#'      data set based on the Savitzky-Golay filter. Remote Sens. Environ. 91, 
#'      332–344. https://doi.org/10.1016/j.rse.2004.03.014
wChen <- function(y, yfit, w, ...){
    if (missing(w)) w  <- rep(1, length(y))
    wnew <- w
    
    re     <- yfit - y
    re_abs <- abs(re)

    d_max  <- max(re_abs, na.rm = T) #6 * median(re, na.rm = T)

    I_pos  <- re > 0
    wnew[ I_pos]  <- (1 - re_abs / d_max) * w[I_pos]
    
    return(wnew)
}

#' 
#' Beck et al., (2006) weigths updating method
#' 
#' @references
#' [1]. Beck, P.S.A., Atzberger, C., Høgda, K.A., Johansen, B., Skidmore, A.K., 
#'      2006. Improved monitoring of vegetation dynamics at very high latitudes: 
#'      A new method using MODIS NDVI. Remote Sens. Environ. 
#'      https://doi.org/10.1016/j.rse.2005.10.021 \cr
#' [2]. https://github.com/kongdd/phenopix/blob/master/R/FitDoubleLogBeck.R
wBECK <- function(y, yfit, w, ...){
    # get optimized parameters
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
