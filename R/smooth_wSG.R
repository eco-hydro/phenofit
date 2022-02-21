sgolayS <- function(frame, d){
    outer((-(frame-1)/2):((frame-1)/2), 0:d, "^")
}

#' Weighted Savitzky-Golay
#'
#' @inheritParams smooth_wHANTS
#' @param frame Savitzky-Golay windows size
#' @param d polynomial of degree. When d = 1, it becomes moving average.
#' @param ylu (optional) `[low, high]` value of time-series y (curve fitting values
#' are constrained in the range of ylu.
#' @inherit smooth_wHANTS return
#' 
#' @references
#' 1. Chen, J., J\"onsson, P., Tamura, M., Gu, Z., Matsushita, B., Eklundh, L.,
#'      2004. A simple method for reconstructing a high-quality NDVI time-series
#'      data set based on the Savitzky-Golay filter. Remote Sens. Environ. 91,
#'      332-344. https://doi.org/10.1016/j.rse.2004.03.014. \cr
#' 2. https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
#'
#' @examples
#' library(phenofit)
#' data("MOD13A1")
#' dt <- tidy_MOD13(MOD13A1$dt)
#' d <- dt[site == "AT-Neu", ]
#'
#' l <- check_input(d$t, d$y, d$w, nptperyear=23)
#' r_wSG <- smooth_wSG(l$y, l$w, l$ylu, nptperyear = 23, iters = 2)
#' @export
smooth_wSG <- function(y, w, nptperyear, ylu, wFUN = wTSM, iters = 2,
                   frame = floor(nptperyear/5)*2 + 1, d=2, ...){
    if (all(is.na(y)) || sum(y, na.rm = TRUE) == 0) {
        return(list(zs = list(ziter1 = y * NA), ws = list(witer1 = y*NA)))
    }
    if (missing(w)) w <- rep(1, length(y))
    if (missing(ylu)) {
        alpha = 0.01
        ylu = quantile(y, c(alpha / 2, 1 - alpha), na.rm = TRUE)
    }
    w[is.na(y)] = 0
    w[y < ylu[1] | y > ylu[2]] = 0

    halfwin <- floor((frame - 1) / 2)
    yiter <- y
    fits  <- list()
    ws    <- list()

    for (i in 1:iters){
        ws[[i]] <- w
        z <- rcpp_wSG(yiter, halfwin, d, w)

        if (!is.null(wFUN)) w <- wFUN(y, z, w, i, nptperyear, ...)
        if (!missing(ylu)) z %<>% checkfit(ylu)
        # z <- check_ylu(z, ylu)
        I = which(yiter < z)
        if (length(I) > 0 ) yiter[I] <- z[I] # upper envelope
        
        fits[[i]] <- z
    }

    fits %<>% set_names(paste0('ziter', 1:iters))
    ws   %<>% set_names(paste0('witer', 1:iters))

    list(zs = fits, ws = ws)
}

checkfit <- function(y, ylu) {
    I_na = y < ylu[1] | y > ylu[2]
    if (sum(I_na) > 0) {
        y[I_na] = approx(which(!I_na), y[!I_na], which(I_na), "constant", rule = 2)$y
    }
    y
}
