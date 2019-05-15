#' @name I_optim
#' @title Interface of unified optimization functions.
#'
#' @description
#' Caution that \pkg{optimx} speed is not so satisfied. So `I_optim`
#' is present.
#'
#' @inheritParams optim_pheno
#' @param FUN Fine curve fitting function for goal function [f_goal()].
#' @param method `method` can be one of `'BFGS','CG','Nelder-Mead',
#' 'L-BFGS-B', 'nlm', 'nlminb', 'ucminf'`. \cr
#' For `I_optimx`, other methods are also supported,
#' e.g. `'spg','Rcgmin','Rvmmin', 'newuoa','bobyqa','nmkb','hjkb'`.
#'
#' @inherit opt_optim return
#'
#' @seealso [stats::optim()], [stats::nlminb()],
#' [stats::nlm()], [optimx::optimx()],
#' [ucminf::ucminf()]
#'
#' @examples
#' library(ggplot2)
#' library(magrittr)
#' library(purrr)
#' 
#' # simulate vegetation time-series
#' fFUN = doubleLog.Beck
#' par = c(
#'   mn  = 0.1,
#'   mx  = 0.7,
#'   sos = 50,
#'   rsp = 0.1,
#'   eos = 250,
#'   rau = 0.1)
#' t    <- seq(1, 365, 8)
#' tout <- seq(1, 365, 1)
#' y <- fFUN(par, t)
#' 
#' # initial parameter
#' par0 <- c(
#'   mn  = 0.15,
#'   mx  = 0.65,
#'   sos = 100,
#'   rsp = 0.12,
#'   eos = 200,
#'   rau = 0.12)
#' 
#' objective <- f_goal # goal function
#' optFUNs   <- c("opt_ucminf", "opt_nlminb", "opt_nlm", "opt_optim") %>% set_names(., .)
#' prior <- as.matrix(par0) %>% t() %>% rbind(., .)
#' 
#' opt1 <- I_optim(prior, fFUN, y, t, tout, c("BFGS", "ucminf", "nlm", "nlminb"))
#' opt2 <- I_optimx(prior, fFUN, y, t, tout, c("BFGS", "ucminf", "nlm", "nlminb"))
#' 
#' # microbenchmark::microbenchmark(
#' #     I_optim(prior, fFUN, y, t, tout, c("BFGS", "ucminf", "nlm", "nlminb")),
#' #     I_optimx(prior, fFUN, y, t, tout, c("BFGS", "ucminf", "nlm", "nlminb")),
#' #     times = 2
#' # )
#' @export
I_optim <- function(prior, FUN, y, t, tout, method = "BFGS", ...)
{
    res <- vector("list", length(method)) %>% set_names(method)

    for (i in seq_along(method)){
        meth <- method[i]
        
        # Get optimization function
        optFUN_name <- ""
        if (meth %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B")){
            optFUN_name <- "opt_optim"
        } else if (meth %in% c("nlm", "nlminb", "ucminf")) {
            optFUN_name <- paste0("opt_", meth)
        } else {
            stop(sprintf("[e] method: %s not supported!", meth))
        }
        optFUN <- get(optFUN_name, mode = "function")

        # browser()
        # optimize parameters
        opt.lst <- alply(prior, 1, optFUN, method = meth,
            objective = f_goal, fun = FUN, y = y, t = t, ...)

        opt.df <- ldply(opt.lst, function(opt){
            with(opt,
                c(par, value = value,
                    fevals = fevals[[1]], niter  = nitns,
                    convcode = convcode)
            )
        }, .id = NULL ) #.id didn't work here
        res[[i]] <- opt.df
    }
    # browser()
    melt_list(res, "method") # quickly return
}

methods <- c(
    'nlm', 'nlminb', 'ucminf',
    "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", # "SANN", "Brent", # (optim)
    'spg','Rcgmin','Rvmmin','newuoa','bobyqa','nmkb','hjkb')

#' @param verbose If `TRUE`, all optimization methods in
#' [optimx::optimx()] are used, and print optimization information
#' of all methods.
#'
#' @rdname I_optim
#' @import optimx
#' @export
I_optimx <- function(prior, FUN, y, t, tout, method, verbose = FALSE, ...){

    # add method column
    opt.df <- alply(prior, 1, optimx, .id = NULL,
        fn = f_goal, fun = FUN, y = y, t = t,
        method = method, ...,
        control = list(maxit = 1000, all.methods = verbose, dowarn = FALSE)
    ) %>% map(~cbind(., method = rownames(.))) %>%
        do.call(rbind, .) %>% set_rownames(NULL)
    opt.df$kkt1 <- NULL
    opt.df$kkt2 <- NULL


    if (verbose){
        opt.df$method <- methods
        opt.df %<>% {.[with(., order(convcode, value, xtimes)), ]}
        print(opt.df) #[, ]

        df <- opt.df[, c("method", "value", "xtimes", "convcode")]#
        # df %<>% reshape2::melt(id.vars = "method", variable.name = "index")
        # df$method %<>% as.factor()
        # ggplot(df, aes(y = method, y = value*1000, fill = method)) +
        #     geom_point() +
        #     scale_y_log10() +
        #     geom_boxplot(outlier.size=2) +
        #     facet_wrap(~index, ncol = 1, scales = "free_y") +
        #     geom_jitter(width = 0.15, size = 1.7, alpha = 1)
        avg <- aggregate(.~method, df, mean)
        avg <- avg[avg$convcode < 0.5, ]
        avg %<>% {.[with(., order(value, xtimes)), ]}
        print(avg)
    }
    return(opt.df)
}
