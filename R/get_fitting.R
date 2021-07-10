#' getFittings
#'
#' Get curve fitting data.frame
#'
#' @inheritParams get_GOF
#'
#' @example inst/examples/ex-get_fitting_param_GOF.R
#' @export
get_fitting <- function(fit){
    lapply(fit, get_fitting.fFITs) %>% melt_list("flag")
}

#' @rdname get_fitting
#'
#' @importFrom purrr map_dfc
#' @export
get_fitting.fFITs <- function(fFITs){
    models = fFITs$model

    t  <- fFITs$data$t
    # fix error: t not in tout
    I  <- match(t, fFITs$tout)
    Ix <- which(!is.na(I))
    I  <- I[Ix]
    t  <- t[Ix]

    iters <- length(models[[1]]$zs)
    df <- models %>% map(function(x){
        d_z <- map(x$zs, ~.[I]) %>% as.data.table() %>%
            set_colnames(paste0("ziter", 1:iters))
        # d_w <- map_dfc(x$ws, ~.) %>% set_colnames(paste0("witer", 1:iters))
        cbind(t, d_z) # , d_w
    }) %>% melt_list("meth") #%>% as.data.table()
    
    df <- merge(fFITs$data[Ix], df, id = "t")
    df$t %<>% as.Date(date.origin)
    df
}

# tasklist
# --------
# 1. ws not exported, need to add I_out variable
