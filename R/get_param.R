#' Get parameters from curve fitting result
#'
#' @param fits Multiple methods curve fitting results by \code{curvefits} result.
#' @inheritParams get_GOF
#'
#' @examples
#' library(phenofit)
#' # simulate vegetation time-series
#' fFUN = doubleLog.Beck
#' par  = c(
#'     mn  = 0.1,
#'     mx  = 0.7,
#'     sos = 50,
#'     rsp = 0.1,
#'     eos = 250,
#'     rau = 0.1)
#' t    <- seq(1, 365, 8)
#' tout <- seq(1, 365, 1)
#' y <- fFUN(par, t)
#'
#' methods <- c("AG", "Beck", "Elmore", "Gu", "Zhang") # "Klos" too slow
#' fFITs <- curvefit(y, t, tout, methods)
#' fits  <- list(`2001` = fFITs, `2002` = fFITs)
#' 
#' d_param <- get_param(fits)
#' @export
get_param <- function(fits){
    llply(fits, get_param.fFITs) %>%
        purrr::transpose() %>%
        map(~ldply(., function(x) x, .id = "flag") %>% as_tibble())
}

#' @rdname get_param
#' @export
get_param.fFITs <- function(fFITs){
    map(fFITs$fFIT, "par")
}
