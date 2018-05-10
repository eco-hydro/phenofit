#' Weighted Whittaker smoothing with a second order finite difference penalty
#'
#' This function smoothes signals with a finite difference penalty of order 2. 
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
#' data(gaschrom)
#' plot(gaschrom[1,], type = "l", ylim = c(0, 100))
#' lines(whit2(gaschrom[1,], lambda = 1e5), col = 2)
#' lines(whit2(gaschrom[1,], lambda = 1e6), col = 3)
#' lines(whit2(gaschrom[1,], lambda = 1e7), col = 4)
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
#' @references
#' [1]. Eilers, P.H.C., 2003. A perfect smoother. Anal. Chem. https://doi.org/10.1021/ac034173t
#' [2]. Frasso, G., Eilers, P.H.C., 2015. L- and V-curves for optimal smoothing. Stat. 
#'      Modelling 15, 91–111. https://doi.org/10.1177/1471082X14549288
#' @export
whitsmw2 <- function(y, w, ylu, nptperyear, wFUN = wTSM, iters=1, lambdas=1000,
    df = NULL, ...)
{
    if (all(is.na(y))) return(y)
    n <- sum(w)

    OUT <- list()
    yiter <- y
    for (j in seq_along(lambdas)){
        lambda <- lambdas[j]

        fits <- list()
        for (i in 1:iters){
            if (i > 1) {
                # print('k1')
                w <- wFUN(y, z, w, iters, nptperyear, ...)
                # w <- phenofit:::wTSM_cpp(y, z, w, iters, nptperyear, 0.5)
            }

            # z_temp <- smooth_HdH(yiter, w, lambda=lambda)$z
            z_temp <- whit2(yiter, lambda, w)
            # If curve has been smoothed enough, it will not care about the
            # second smooth. If no, the second one is just prepared for this
            # situation.
            z <- whit2(z_temp, w, lambda=lambda) #genius move
            z <- check_fit(z, ylu)

            yiter[yiter < z] <- z[yiter < z]
            # if (!validation){
            # }
            fits[[i]] <- z
            # yiter <- z# update y with smooth values
        }
        fits %<>% set_names(paste0('iter', 1:iters))

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
        OUT[[j]] <- list(data = as_tibble(c(list(w = w), fits)))
    }
    if (length(lambdas) == 1) OUT <- OUT[[1]]
    return(OUT)
}
