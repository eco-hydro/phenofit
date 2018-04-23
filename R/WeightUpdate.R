# .onUnload <- function (libpath) {
#   library.dynam.unload("rTIMESAT", libpath)
# }

#' bisquare weight update method
#'
#' @references
#' [1]. https://au.mathworks.com/help/curvefit/smoothing-data.html#bq_6ys3-3
#' [2]. Garcia, D., 2010. Robust smoothing of gridded data in one and higher
#' dimensions with missing values. Computational statistics & data analysis,
#' 54(4), pp.1167-1178.
#' @export
bisquare <- function(y, yfit, w, ...){
    if (missing(w))
        w  <- rep(0, length(y))

    re <- abs(y - yfit)
    sc <- 6 * median(re, na.rm = T)

    I <- which(re < sc)
    w[I]  <- (1 - (re[I]/sc)^2)^2
    # w[!I] <- 0
    # constrain growing VI: enlarge the positive bias values, reduce the weights
    #   of negative bias values, as TIMESAT
    # diff = 2 * (yfit(i) - y(i)) / yfitstd;
    # w <- wfact * wfit(i) * exp( - ydiff ^ 2);

    return(w)
}

# https://github.com/kongdd/phenopix/blob/master/R/FitDoubleLogBeck.R
wBECK <- function(){
    # get optimized parameters
    sos <- opt$par[3]
    eos <- opt$par[5]

    m <- lm(c(0, 100) ~ c(sos, eos))
    tr <- coef(m)[2] * t + coef(m)[1]
    tr[tr < 0] <- 0
    tr[tr > 100] <- 100
            
    # estimate weights
    res <- xpred - x
    weights <- 1/((tr * res + 1)^2)
    weights[res > 0 & res <= 0.01] <- 1
    weights[res < 0] <- 4
}

wSELF <- function(y, z, w, ...){w}