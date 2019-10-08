methods <- c(
    'nlm', 'nlminb', 'ucminf',
    "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", # "SANN", "Brent", # (optim)
    'spg','Rcgmin','Rvmmin','newuoa','bobyqa','nmkb','hjkb')

#' @name I_optim
#' @title Interface of unified optimization functions.
#'
#' @description
#' \pkg{optimx} speed is not satisfied. So `I_optim` is present.
#'
#' - `I_optim`: Interface of unified optimization functions.
#' - `I_optimx`: deprecated, which is about 10 times slower than `I_optim`.
#'
#' @inheritParams optim_pheno
#' @param FUN Fine curve fitting function for goal function [f_goal()].
#' @param method `method` can be some of `'BFGS','CG','Nelder-Mead',
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
#' @example man/examples/ex-I_optim.R
#'
#' @keywords internal
#' @export
I_optim <- function(prior, FUN, y, t, method = "BFGS", fn = f_goal, ...)
{
    pred = y*0
    if (is.vector(prior)) prior <- t(prior)

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

        # optimize parameters
        # browser()
        opt.lst <- alply(prior, 1, optFUN, method = meth,
            objective = fn, fun = FUN, y = y, t = t, pred = pred, ...)
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

#' @param verbose If `TRUE`, all optimization methods in
#' [optimx::optimx()] are used, and print optimization information
#' of all methods.
#'
#' @rdname I_optim
#' @import optimx
#'
#' @keywords internal
#' @export
I_optimx <- function(prior, FUN, y, t, method, verbose = FALSE, ...){
    pred = y*0
    if (is.vector(prior)) prior <- t(prior)

    # add method column
    opt.lst <- alply(prior, 1, optimx, method = method,
        fn = f_goal, fun = FUN, y = y, t = t, pred = pred, ...,
        control = list(maxit = 1000, all.methods = verbose, dowarn = FALSE)
    )
    opt.df <- map(opt.lst, ~cbind(., method = rownames(.))) %>%
        do.call(rbind, .) %>% set_rownames(NULL)
    opt.df$kkt1 <- NULL
    opt.df$kkt2 <- NULL

    if (verbose){
        opt.df$method <- methods
        opt.df %<>% {.[with(., order(convcode, value, xtimes)), ]}
        print(opt.df)

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
