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
    ws    <- list()

    for (i in 1:iters){
        ws[[i]] <- w
        z <- sgfitw_rcpp(yiter, w, S)[, 1]
        wnew <- wFUN(y, z, w, i, nptperyear, ...)

        z <- check_fit(z, ylu)
        yiter[yiter < z] <- z[yiter < z] # upper envelope
        
        fits[[i]] <- z
    }

    fits %<>% set_names(paste0('ziter', 1:iters))
    ws   %<>% set_names(paste0('witer', 1:iters))
    
    list(ws = ws, zs = fits)
}
