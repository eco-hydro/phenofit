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
#' @param use.julia whether use julia nlminb optimization?
#'
#' @inherit opt_optim return
#'
#' @seealso [stats::optim()], [stats::nlminb()],
#' [stats::nlm()], [optimx::optimx()],
#' [ucminf::ucminf()]
#'
#' @example R/examples/ex-I_optim.R
#'
#' @keywords internal
#' @export
I_optim <- function(prior, FUN, y, t, method = "BFGS", fn = f_goal, ...,
    use.julia = FALSE)
{
    pred = y*0
    if (is.vector(prior)) prior <- t(prior)
    prior2 <- set_colnames(prior, NULL)

    colnames <- colnames(prior) %>% c("value", "fevals", "niter", "convcode")

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

        lst <- list()
        for (j in 1:nrow(prior)) {
            par0 <- prior2[j, ]

            if (optFUN_name == "opt_nlminb" && use.julia) {
                FUN_name = attr(FUN, "name")
                r   <- opt_nlminb_julia(par0, FUN_name, y, t, ...) # w, ylu
                ans <- c(r$par, value = r$objective,
                    fevals   = r$evaluations[[1]],
                    niter    = r$iterations,
                    convcode = r$convergence)
            } else {
                optFUN <- get(optFUN_name, mode = "function")
                r   <-  optFUN(par0, method = meth,
                    objective = fn, fun = FUN, y = y, t = t, pred = pred, ...)
                ans <- c(r$par, value = r$value,
                    fevals   = r$fevals[[1]], niter  = r$nitns,
                    convcode = r$convcode)
            }
            lst[[j]] <- ans
        }
        res[[i]] <- do.call(rbind, lst) # %>% as.data.frame() %>% set_colnames(colnames)
    }
    res = do.call(rbind, res)
    dimnames(res) <- list(rep(method, each = nrow(prior)), colnames)
    res
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
    opt.lst <- map(1:nrow(prior), function(i) {
        optimx(prior[i, ],
            method = method, fn = f_goal, fun = FUN, y = y, t = t, pred = pred, ...,
            control = list(maxit = 1000, all.methods = verbose, dowarn = FALSE)
        )
    })
    opt.df <- map(opt.lst, ~cbind(., method = rownames(.))) %>%
        do.call(rbind, .) %>% set_rownames(NULL)
    opt.df$kkt1 <- NULL
    opt.df$kkt2 <- NULL

    if (verbose){
        opt.df$method <- methods
        opt.df %<>% {.[with(., order(convcode, value, xtimes)), ]}
        print(opt.df)

        df <- opt.df[, c("method", "value", "xtimes", "convcode")]#
        avg <- aggregate(.~method, df, mean)
        avg <- avg[avg$convcode < 0.5, ]
        avg %<>% {.[with(., order(value, xtimes)), ]}
        print(avg)
    }
    return(opt.df)
}
