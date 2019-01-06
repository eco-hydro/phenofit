#' getFittings
#'
#' Get curve fitting data.frame
#'
#' @inheritParams get_GOF
#'
#' @example inst/examples/ex-get_fitting_param_GOF.R
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
