#' @export
sgolayS <- function(frame, d){
    outer((-(frame-1)/2):((frame-1)/2), 0:d, "^")
}

#' Weighted Savitzky-Golay
#' 
#' @inheritParams wHANTS
#' @param frame Savitzky-Golay windows size
#' @param d polynomial of degree
#' 
#' @export
sgfitw <- function(y, w, nptperyear, ylu, wFUN = wTSM, iters = 2,
                   frame = floor(nptperyear/7)*2 + 1, d=2, ...){
    if (all(is.na(y))) return(y)

    S <- sgolayS(frame, d)
    
    yiter <- y
    fits  <- list()
    for (i in 1:iters){
        z <- sgfitw_rcpp(yiter, w, S)[, 1]
        w <- wFUN(y, z, w, i, nptperyear, ...)

        z <- check_fit(z, ylu)
        yiter[yiter < z] <- z[yiter < z] # upper envelope
        
        fits[[i]] <- z
    }

    fits %<>% set_names(paste0('iter', 1:iters))
    c(list(w = w), fits)
}
