#' optim_pheno
#'
#' Interface of optimization functions for double logistics and other parametric
#' curve fitting functions.
#'
#' @inheritParams wHANTS
#' @param prior A vector of initial values for the parameters for which optimal
#' values are to be found. \code{prior} is suggested giving a column name.
#' @param sFUN The name of fine curve fitting functions, can be one of \code{
#' 'FitAG', 'FitDL.Beck', 'FitDL.Elmore', 'FitDL.Gu' and 'FitDL.Klos',
#' 'FitDL.Zhang'}.
#' @param nptperyear Integer, number of images per year, passed to \code{wFUN}.
#' Only \code{\link{wTSM}} needs \code{nptperyear}. If not specified,
#' \code{nptperyear} will be calculated based on \code{t}.
#' @param ylu \code{ymin, ymax}, which is used to force \code{ypred} in the
#' range of \code{ylu}.
#' @param tout Corresponding doy of prediction.
#' @param I_optimFUN Interface of optimization function, can be one of
#' \code{\link{I_optim}} and \code{\link{I_optimx}}.
#' @param method The name of optimization method, passed to \code{I_optimFUN}.
#' @param verbose Whether to display intermediate variables?
#' @param ... other parameters passed to \code{I_optimFUN}
#'
#' @return fFIT object, see \code{\link{fFIT}} for details.
#'
#' @import optimx
#' @export
optim_pheno <- function(
    prior, sFUN,
    y, t, tout,
    I_optimFUN = I_optim, method,
    w, nptperyear, ylu,
    iters = 2, wFUN = wTSM,
    verbose = FALSE, ...)
{
    # To improve the performance of optimization, \code{t} needs to be normalized.
    # t_0    <- t
    tout_0 <- tout
    # tout   <- tout - t[1] + 1
    # t      <- t - t[1] + 1

    # check input parameters
    if (missing(w))          w   <- rep(1, length(y))
    if (missing(nptperyear)) nptperyear <- ceiling(365/mean(as.numeric(diff(t))))
    if (missing(ylu))        ylu <- range(y, na.rm = TRUE)
    A <- diff(ylu) # y amplitude

    FUN <- get(sFUN, mode = "function" )
    # add prior colnames
    parnames <- attr(FUN, 'par')
    colnames(prior) <- parnames
    ############################################################################

    fits <- list()
    ws   <- list() #record weights for every iteration
    ## 1. add weights updating methods here
    par   <- setNames(numeric(length(parnames))*NA, parnames)
    ypred <- rep(NA, length(tout))

    if (verbose){
        boundary <- list(...)[c('lower', 'upper')] %>% do.call(rbind, .) %>%
            set_colnames(parnames) %>% as_tibble()
        print(boundary)
    }

    for (i in 1:iters){
        ws[[i]] <- w
        # pass verbose for optimx optimization methods selection
        opt.df  <- I_optimFUN(prior, FUN, y, t, tout, method = method, w = w, verbose = verbose, ...)
        # colnames(opt.df) <- c(parnames, colnames(opt.df)[(length(parnames) + 1):ncol(opt.df)])
        if (verbose){
            fprintf('Initial parameters:\n')
            print(as_tibble(prior))
            print(opt.df)
        }

        best <- which.min(opt.df$value)
        if (is_empty(best)){
            warning("optimize parameters failed!")
        }else{
            # test for convergence if maximum iterations where reached - restart from best with more iterations
            # If RMSE is small enough, then accept it.
            if (opt.df$convcode[best] != 0) {
                # repeat with more maximum iterations if it did not converge
                par <- opt.df[best, parnames] %>% as.matrix #best par generally
                opt <- I_optimFUN(par, FUN, y, t, tout, method, w = w, verbose = verbose, ...)
            } else { #if (opt.df$convcode[best] == 0)
                opt <- opt.df[best, ]
            }

            # check whether convergence
            if (opt$convcode != 0 & opt$value > 0.1*A) {
                warning("Not convergent!")
            }else{
                par   <- opt[, parnames] %>% unlist()
                # put opt par into prior for the next iteration
                prior <- rbind(prior[1:(nrow(prior)-1), ], par, deparse.level = 0)
                ypred <- FUN(par, tout)
                # too much missing values
                # if (sum(w == 0)/length(w) > 0.5) ypred <- ypred*NA

                # fixed 04 March, 2018;
                w <- tryCatch(
                    # to adapt wTS, set iter = i-1; #20180910
                    wFUN(y, FUN(par, t), w, i, nptperyear, ...),
                    #nptperyear, wfact = 0.5)
                    error = function(e) {
                        message(sprintf('[%s]: %s', sFUN, e$message))
                        return(w) #return original w
                    })
                ypred %<>% check_fit(ylu) #values out of ylu are set as NA
            }
        }
        fits[[i]] <- ypred
    }
    fits %<>% set_names(paste0('iter', 1:iters))
    ws   %<>% set_names(paste0('iter', 1:iters)) # %>% as_tibble()

    # 02. uncertain part also could add here
    # returned object FUN need to be futher optimized
    structure(list(tout = tout_0, zs = fits, ws = ws,
        par  = par, fun = sFUN), class = 'fFIT')
}

# 1. Firstly, using \pkg{optimx} package to select the best performance
#   optimization function.
# 2. Secondly, even though optimx is extraordinary, but it running is not so
#   satisfied, so using self unified optimization function (with prefixion of
#   'p_', opt_FUN means phenology optimization function) to replace it through
#   `optim_p` function. The unified procedure was inspired by `optimx` package.

# 1. nlfr package
# start <- set_names(prior[1, ], attr(FUN, "par"))
# # lower = lower, upper = upper
# f_nlf <- nlfb(start, resfn = function(par, t, y, fun, ...){
#     fun(par, t) - y
# }, trace = FALSE, y = y, t=t, w = w, fun = FUN, ...)
# yfit  <- FUN(f_nlf$coefficients, t)
# plot(t, y, type = "b")
# lines(t, yfit)
# print(f_nlf)
