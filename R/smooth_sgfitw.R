#' @export
sgolayS <- function(frame, d){
    outer((-(frame-1)/2):((frame-1)/2), 0:d, "^")
}

#' @export
sgfitw <- function(y, w, nptperyear, ylu, wFUN = wTSM, iters = 2,
                   frame = floor(nptperyear/7)*2 + 1, d=2, ...){
    if (all(is.na(y))) return(y)

    S <- sgolayS(frame, d)

    fits <- list()
    for (i in 1:iters){
        z <- sgfitw_rcpp(y, w, S)[, 1]
        w <- wFUN(y, z, w, i, nptperyear, ...)
        z <- check_fit(z, ylu)
        fits[[i]] <- z
    }

    fits %<>% set_names(paste0('iter', 1:iters))
    return(as_tibble(c(list(w = w), fits)))
}

# test about those function using GPP data.
