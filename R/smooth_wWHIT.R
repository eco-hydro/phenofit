#' Weighted Whittaker smoothing with a second order finite difference penalty
#'
#' This function smoothes signals with a finite difference penalty of order 2.
#' This function is modified from `ptw` package.
#'
#' @param y signal to be smoothed: a vector
#' @param lambda smoothing parameter: larger values lead to more smoothing
#' @param w weights: a vector of same length as y. Default weights are equal to one
#'
#' @return A numeric vector, smoothed signal.
#'
#' @references
#' 1. Eilers, P.H.C. (2004) "Parametric Time Warping", Analytical Chemistry,
#' **76** (2), 404 -- 411. \cr
#' 2. Eilers, P.H.C. (2003) "A perfect smoother", Analytical Chemistry,
#' **75**, 3631 -- 3636.
#'
#' @author Paul Eilers, Jan Gerretzen
#' @examples
#' library(phenofit)
#' data("MOD13A1")
#' dt <- tidy_MOD13.gee(MOD13A1$dt)
#' y <- dt[site == "AT-Neu", ][1:120, y]
#' 
#' plot(y, type = "b")
#' lines(whit2(y, lambda = 2), col = 2)
#' lines(whit2(y, lambda = 10), col = 3)
#' lines(whit2(y, lambda = 100), col = 4)
#' legend("bottomleft", paste("lambda = ", c(2, 10, 15)), col = 2:4, lty = rep(1, 3))
#' @export
whit2 <- function(y, lambda, w = rep(1, ny))
{
    ny <- length(y)
    z <- d <- c <- e <- rep(0, length(y))
    # smooth2(
    .C("smooth2",
         w = as.double(w),
         y = as.double(y),
         z = as.double(z),
         lamb = as.double(lambda),
         mm = as.integer(length(y)),
         d = as.double(d),
         c = as.double(c),
         e = as.double(e))$z
}

#' Weigthed Whittaker Smoother
#'
#' @inheritParams smooth_wHANTS
#' @param lambda whittaker parameter (2-15 is suitable for 16-day VI). Multiple
#' lambda values also are accept, then a list object return.
#' @param second If true, in every iteration, Whittaker will be implemented
#' twice to make sure curve fitting is smooth. If curve has been smoothed
#' enough, it will not care about the second smooth. If no, the second one is
#' just prepared for this situation. If lambda value has been optimized, second
#' smoothing is unnecessary.
#'
#' @inherit smooth_wHANTS return
#' 
#' @references
#' 1. Eilers, P.H.C., 2003. A perfect smoother. Anal. Chem. https://doi.org/10.1021/ac034173t \cr
#' 2. Frasso, G., Eilers, P.H.C., 2015. L- and V-curves for optimal smoothing. Stat.
#'      Modelling 15, 91-111. https://doi.org/10.1177/1471082X14549288
#' 
#' @examples
#' library(phenofit)
#' data("MOD13A1")
#' dt <- tidy_MOD13.gee(MOD13A1$dt)
#' d <- dt[site == "AT-Neu", ]
#' 
#' l <- check_input(d$t, d$y, d$w, nptperyear=23)
#' r_wWHIT <- smooth_wWHIT(l$y, l$w, l$ylu, nptperyear = 23, iters = 2)
#' @export
smooth_wWHIT <- function(y, w, ylu, nptperyear, wFUN = wTSM, iters=1, lambda=15,
    second = FALSE, ...) #, df = NULL
{
    trs <- 0.5
    if (all(is.na(y))) return(y)
    n <- sum(w)

    OUT   <- list()
    yiter <- y
    for (j in seq_along(lambda)){
        lambda_j <- lambda[j]

        fits <- list()
        ws   <- list()
        for (i in 1:iters){
            ws[[i]] <- w
            z <- whit2(yiter, lambda_j, w)
            w <- wFUN(y, z, w, i, nptperyear, ...)

            # If curve has been smoothed enough, it will not care about the
            # second smooth. If no, the second one is just prepared for this
            # situation.
            if (second) z <- whit2(z, lambda_j, w) #genius move

            ## Based on our test, check_ylu and upper envelope will decrease 
            # `wWWHIT`'s performance (20181029). 
            z <- check_ylu(z, ylu) # check ylu

            ylu <- range(z)
            zc <- ylu[1] + (ylu[2] - ylu[1])*trs

            I_fix <- z > yiter & z > zc
            yiter[I_fix] <- z[I_fix] # upper envelope
            
            # browser()
            fits[[i]] <- z
            # wnew <- wFUN(y, z, w, i, nptperyear, ...)
            # yiter <- z# update y with smooth values 
        }
        fits %<>% set_names(paste0('ziter', 1:iters))
        ws   %<>% set_names(paste0('witer', 1:iters))

        OUT[[j]] <- list(zs = fits, ws = ws)
    }
    if (length(lambda) == 1) OUT <- OUT[[1]]
    return(OUT)
}

# # CROSS validation
# if (validation){
#     h   <- fit$dhat
#     df  <- sum(h)
#     r   <- (y - z)/(1 - h)
#     cv  <- sqrt( sum( r^2*w ) /n )
#     gcv <- sqrt( sum( (r/(n-df))^2*w ))
#     LV  <- whit_V(y, z, w) #L curve, D is missing now
#     OUT[[j]] <- c(list(data = as_tibble(c(list(w = w), fits)),
#         df = df, cv = cv, gcv = gcv), LV)
# }else{
# }
