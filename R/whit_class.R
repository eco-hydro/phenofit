#' Weigthed Whittaker Smoother
#'
#' @references
#' [1]. Eilers, P.H.C., 2003. A perfect smoother. Anal. Chem. https://doi.org/10.1021/ac034173t
#' [2]. Frasso, G., Eilers, P.H.C., 2015. L- and V-curves for optimal smoothing. Stat. 
#'      Modelling 15, 91–111. https://doi.org/10.1177/1471082X14549288
#' @importFrom ptw whit2
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
