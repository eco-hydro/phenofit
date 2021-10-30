#' getFittings
#'
#' Get curve fitting data.frame
#'
#' @inheritParams get_param
#'
#' @example inst/examples/ex-get_fitting_param_GOF.R
#' @export
get_fitting <- function(x) UseMethod("get_fitting", x)

#' @rdname get_fitting
#' @export
get_fitting.list <- function(x){
    lapply(x, get_fitting.fFITs) %>% melt_list("flag")
}

#' @importFrom purrr map_dfc
#' @rdname get_fitting
#' @export
get_fitting.fFITs <- function(x){
    models = x$model

    t  <- x$data$t
    # fix error: t not in tout
    I  <- match(t, x$tout)
    Ix <- which(!is.na(I))
    I  <- I[Ix]
    t  <- t[Ix]

    iters <- length(models[[1]]$zs)
    df <- models %>% map(function(model){
        d_z <- map(model$zs, ~.[I]) %>% as.data.table() %>%
            set_colnames(paste0("ziter", 1:iters))
        # d_w <- map_dfc(x$ws, ~.) %>% set_colnames(paste0("witer", 1:iters))
        cbind(t, d_z) # , d_w
    }) %>% melt_list("meth") #%>% as.data.table()
    
    df <- merge(x$data[Ix], df, id = "t")
    df$t %<>% as.Date(date.origin)
    df
}

# tasklist
# --------
# 1. ws not exported, need to add I_out variable
