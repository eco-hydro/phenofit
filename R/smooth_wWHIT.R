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
#' [1]. Eilers, P.H.C. (2004) "Parametric Time Warping", Analytical Chemistry,
#' \bold{76} (2), 404 -- 411. \cr
#' [2]. Eilers, P.H.C. (2003) "A perfect smoother", Analytical Chemistry,
#' \bold{75}, 3631 -- 3636.
#'
#' @author Paul Eilers, Jan Gerretzen
#' @examples
#' \dontrun{
#' data(gaschrom)
#' plot(gaschrom[1,], type = "l", ylim = c(0, 100))
#' lines(whit2(gaschrom[1,], lambda = 1e5), col = 2)
#' lines(whit2(gaschrom[1,], lambda = 1e6), col = 3)
#' lines(whit2(gaschrom[1,], lambda = 1e7), col = 4)
#' }
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
#' @inheritParams wHANTS
#' @param lambda whittaker parameter (2-15 is suitable for 16-day VI). Multiple
#' lambda values also are accept, then a list object return.
#' @param second If true, in every iteration, Whittaker will be implemented
#' twice to make sure curve fitting is smooth. If curve has been smoothed
#' enough, it will not care about the second smooth. If no, the second one is
#' just prepared for this situation. If lambda value has been optimized, second
#' smoothing is unnecessary.
#'
#' @references
#' [1]. Eilers, P.H.C., 2003. A perfect smoother. Anal. Chem. https://doi.org/10.1021/ac034173t \cr
#' [2]. Frasso, G., Eilers, P.H.C., 2015. L- and V-curves for optimal smoothing. Stat.
#'      Modelling 15, 91–111. https://doi.org/10.1177/1471082X14549288
#' @export
wWHIT <- function(y, w, ylu, nptperyear, wFUN = wTSM, iters=1, lambda=15,
    second = FALSE, ...) #, df = NULL
{
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

            ## Based on our test, check_fit and upper envelope will decrease 
            # `wWWHIT`'s performance (20181029). 
            
            # z <- check_fit(z, ylu)
            # yiter[yiter < z] <- z[yiter < z] # upper envelope
            fits[[i]] <- z
            # wnew <- wFUN(y, z, w, i, nptperyear, ...)
            # yiter <- z# update y with smooth values 
        }
        fits %<>% set_names(paste0('ziter', 1:iters))
        ws   %<>% set_names(paste0('witer', 1:iters))

        OUT[[j]] <- list(ws = ws, zs = fits)
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
