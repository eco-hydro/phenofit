#' getFittings
#'
#' Get curve fitting data.frame
#'
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
#'
#' # multiple years
#' fits <- list(`2001` = fFITs, `2002` = fFITs)
#' pheno <- PhenoExtract(fits, "AG", IsPlot=TRUE)
#' @export
get_fitting <- function(fit){
    llply(fit, get_fitting.fFITs) %>% melt_list("flag")
}

#' @rdname get_fitting
#' 
#' @importFrom purrr map_dfc
#' @export
get_fitting.fFITs <- function(fFITs){
    t <- fFITs$data$t
    I <- match(t, fFITs$tout)

    iters <- length(fFITs$fFIT[[1]]$zs)
    df <- fFITs$fFIT %>% map(function(x){
        d_z <- map_dfc(x$zs, ~.[I]) %>% set_colnames(paste0("ziter", 1:iters))
        # d_w <- map_dfc(x$ws, ~.) %>% set_colnames(paste0("witer", 1:iters))
        cbind(t, d_z) # , d_w
    }) %>% melt_list("meth") %>% as.data.table()

    df <- merge(fFITs$data, df, id = "t")
    df$t %<>% as.Date(date.origin)
    df
}

# tasklist
# --------
# 1. ws not exported, need to add I_out variable
