# 1. Firstly, using optimx package to select the best performance optimization
#   function.
# 2. Secondly, even though optimx is extraordinary, but it running is not so
#   satisfied, so using self unified optimization function (with prefixion of
#   'p_', pfun means phenology optimization function) to replace it through
#   `optim_p` function. The unified procedure was inspired by `optimx` package.
#
# Optimization function for original `phenopix` package. The format of its kenal
#   double logistics function have been modified and unified.
#
#' optim_pheno
#'
#' Interface of optimization functions for double logistics and other parametric
#' curve fitting functions
#' 
#' @inheritParams wHANTS
#' @param prior A vector of initial values for the parameters for which optimal 
#' values are to be found.
#' @param FUN_name Curve fitting function name, can be one of 'FitDL.Zhang', 
#' 'FitAG', 'FitDL.Beck', 'FitDL.Elmore', 'FitDL.Gu' and 'FitDL.Klos'.
#' @param tout corresponding \code{t} of curve fitting result
#' @param optimFun optimization function, can be \code{I_optim} or \code{I_optimx}.
#' @param method String, optimization method, passed to `optimFun`.
#' @param debug boolean
#' @param ... other parameters passed to \code{optimFUN}
#' 
#' @return list(pred, par, fun)
#' @import optimx
#' @export
optim_pheno <- function(prior, FUN_name, y, t, tout, optimFUN = I_optim, method,
    w, ylu, iters = 2, wFUN = wTSM, nptperyear = 46, debug = FALSE, ...)
{
    FUN <- get(FUN_name, mode = "function" )
    # add prior colnames
    parnames <- attr(FUN, 'par')
    colnames(prior) <- parnames
    A <- diff(ylu) # y amplitude

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

    fits <- list()
    ws   <- list() #record weights for every iteration
    ## 1. add weights updating methods here
    par   <- setNames(numeric(length(parnames))*NA, parnames)
    xpred <- rep(NA, length(tout))

    if (debug){
        boundary <- list(...)[c('lower', 'upper')] %>% do.call(rbind, .) %>%
            set_colnames(parnames) %>% as_tibble()
        # upper <- list(...)$upper
        print(boundary)
    }

    for (i in 1:iters){
        ws[[i]] <- w
        # pass debug for optimx optimization methods selection
        opt.df  <- optimFUN(prior, FUN, y, t, tout, method, w = w, debug = debug, ...)
        # colnames(opt.df) <- c(parnames, colnames(opt.df)[(length(parnames) + 1):ncol(opt.df)])
        if (debug){
            fprintf('Initial parameters:\n')
            print(as_tibble(prior))
            print(opt.df)
        }

        best    <- which.min(opt.df$value)
        if (is_empty(best)){
            warning("optimize parameters failed!")
        }else{
            # test for convergence if maximum iterations where reached - restart from best with more iterations
            # If RMSE is small enough, then accept it.
            if (opt.df$convcode[best] != 0) {
                # repeat with more maximum iterations if it did not converge
                par <- opt.df[best, parnames] %>% as.matrix #best par generally
                opt <- optimFUN(par, FUN, y, t, tout, method, w = w, debug = debug, ...)
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
                xpred <- FUN(par, tout)
                # too much missing values
                # if (sum(w == 0)/length(w) > 0.5) xpred <- xpred*NA

                # fixed 04 March, 2018; 
                w <- tryCatch(
                    # to adapt wTS, set iter = i-1; #20180910
                    wFUN(y, FUN(par, t), w, i, nptperyear, ...), 
                    #nptperyear, wfact = 0.5)
                    error = function(e) {
                        message(sprintf('[%s]: %s', FUN_name, e$message))
                        return(w) #return original w
                    })
                xpred %<>% check_fit(ylu) #values out of ylu are set as NA
            }
        }
        fits[[i]] <- xpred
    }
    fits %<>% set_names(paste0('iter', 1:iters))
    ws   %<>% set_names(paste0('iter', 1:iters)) # %>% as_tibble()
    
    # 02. uncertain part also could add here
    # returned object FUN need to be futher optimized
    structure(list(tout = tout, fits = fits, ws = ws,
        par  = par, fun = FUN_name), class = 'phenofit')
}

#' I_optim
#'
#' Interface of self defined optimization functions.
#'
#' @inheritParams optim_pheno
#' @param FUN curve fitting function of `f_goal`
#' @param pfun with prefixion of 'p_', pfun was self modified and unified
#'  optimization function for double logistics and many other functions.
#' @param method was only used for `p_optim` function. Other pfun, e.g. p_nlm, p_nlminb
#'  will ignore it.
#' @export
I_optim <- function(prior, FUN, y, t, tout, pfun = p_optim, method = "BFGS",...){
    opt.lst <- alply(prior, 1, pfun, method = method,
        objective = f_goal, fun = FUN, y = y, t = t, ...)
    opt.df <- ldply(opt.lst, with, .id = NULL,
        expr = c(par, value = value,
            fevals = fevals, gevals = gevals, niter  = nitns,
            convcode = convcode)) #.id didn't work here
    return(opt.df)
}
methods <- c('BFGS','CG','Nelder-Mead','L-BFGS-B','nlm','nlminb',
             'spg','ucminf','Rcgmin','Rvmmin','newuoa','bobyqa','nmkb','hjkb')

#' I_optimx
#'
#' Interface of optimx package which unified optimization functions.
#' Should be caution that optimx speed is not so satisfied. So `I_optim` is also
#' present.
#'
#' @inheritParams optim_pheno
#' @param FUN curve fitting function of `f_goal`
#' @param method Method passed to optimx function, which can be one of:
#'  'BFGS','CG','Nelder-Mead','L-BFGS-B','nlm','nlminb', spg','ucminf',
#'  'Rcgmin','Rvmmin','newuoa','bobyqa','nmkb','hjkb'
#' @param debug If debug is true, all.methods is also true for optimx
#' 
#' @return result same as what optimx function returned
#' @import optimx
#' @export
I_optimx <- function(prior, FUN, y, t, tout, method, debug = FALSE, ...){
    # debug = TRUE
    opt.df <- adply(prior, 1, optimx, .id = NULL,
        fn = f_goal, fun = FUN, y = y, t = t,
        method = method, ...,
        control = list(maxit = 1000, all.methods = debug, dowarn = FALSE)
    )
    if (debug){
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

#' Unified optim function
#' 
#' \itemize{
#'    \item \code{p_nlminb} Optimization using PORT routines, see 
#'                          \code{\link[stats]{nlminb}}.
#'    \item \code{p_ncminf} General-Purpose Unconstrained Non-Linear 
#'                          Optimization, see \code{\link[ucminf]{ucminf}}.
#'    \item \code{p_nlm} Non-Linear Minimization, \code{\link[stats]{nlm}}.
#'    \item \code{p_optim} General-purpose Optimization, see 
#'                         \code{\link[stats]{optim}}.
#' }
#' 
#' @param par0 Initial values for the parameters to be optimized over.
#' @param objective A function to be minimized (or maximized), with first 
#' argument the vector of parameters over which minimization is to take place. 
#' It should return a scalar result.
#' @param method optimization method to be used in \code{p_optim}. See 
#' \code{\link[stats]{optim}}.
#' @param ... other parameters passed to \code{objective}.
#' 
#' @return
#' \describe{
#'    \item{convcode}{ An integer code. 0 indicates successful convergence. 
#'    Various methods may or may not return sufficient information to allow all 
#'    the codes to be specified. An incomplete list of codes includes }
#'    \item{value}{ The value of fn corresponding to par }
#'    \item{par}{ The best parameter found }
#'    \item{nitns}{ the number of iterations }
#'    \item{fevals}{ The number of calls to fn. }
#'    \item{gevals}{ The number of calls to gr. This excludes those calls needed
#'    to compute the Hessian, if requested, and any calls to fn to compute a 
#'    finite-difference approximation to the gradient. }
#' }
#' @rdname optim
#' @export
p_nlminb <- function(par0, objective, ...){
    npar <- length(par0)
    ans <- try(nlminb(start=par0, objective=objective, ...,
      control = list(eval.max=1000, iter.max=1000, trace=0, abs.tol=0)), 
                     silent=TRUE)
    # unify the format of optimization results
    if (class(ans)[1] != "try-error") {
        ans$convcode    <- ans$convergence
        # Translate output to common format and names
        ans$value       <- ans$objective
        ans$objective   <- NULL
        ans$fevals      <- ans$evaluations[[1]]
        ans$gevals      <- ans$evaluations[[2]]
        ans$evaluations <- NULL  # cleanup
        ans$nitns       <- ans$iterations
        ans$iterations  <- NULL
    } else {
        # bad result -- What to do?
        ans <- list(fevals = NA)  # ans not yet defined, so set as list
        ans$convcode <- 9999  # failed in run
        # if (ctrl$trace > 0)
          cat("nlminb function evaluation failure\n")
        ans$value  <- NA#ctrl$badval
        ans$par    <- rep(NA, npar)
        ans$nitns  <- NA  # not used
        ans$gevals <- NA  ## ?? missing 130929
    }
    ans$convergence <- NULL
    return(ans)
}

#' @rdname optim
#' @importFrom ucminf ucminf
#' @export
p_ncminf <- function(par0, objective, ...){
    npar <- length(par0)
    ans <- try(ucminf::ucminf(par = par0, fn = objective,
      control = list(maxeval = 1000), ...), silent = TRUE)
    if (class(ans)[1] != "try-error") {
      ans$convcode <- ans$convergence
      # From ucminf documentation:
      # convergence =
      # 1 Stopped by small gradient (grtol).
      # 2 Stopped by small step (xtol).
      # 3 Stopped by function evaluation limit (maxeval).
      # 4 Stopped by zero step from line search
      # -2 Computation did not start: length(par) = 0.
      # -4 Computation did not start: stepmax is too small.
      # -5 Computation did not start: grtol or xtol <= 0.
      # -6 Computation did not start: maxeval <= 0.
      # -7 Computation did not start: given Hessian not pos. definite.  message: String with reason of
      # termination.
      if (ans$convcode == 1 || ans$convcode == 2 || ans$convcode == 4) {
        ans$convcode <- 0
      } else {
        ans$convcode <- ans$convergence
      }  # Termination criteria are tricky here!  How to determine successful convergence?
      ans$fevals <- ans$info[4]
      ans$gevals <- ans$info[4]  # calls fn and gr together
      ans$info   <- NULL  # to erase conflicting name
      ans$nitns  <- NA
      # if (ctrl$trace > 0)
      # cat("ucminf message:", ans$message, "\n")
    } else {
      # ucminf failed
      # if (ctrl$trace > 0)
      cat("ucminf failed for this problem\n")
      ans <- list(fevals = NA)  # ans not yet defined, so set as list
      ans$value = NA#ctrl$badval
      ans$par      <- rep(NA, npar)
      ans$convcode <- 9999  # failed in run
      ans$gevals   <- NA
      ans$nitns    <- NA
    }
    ans$convergence <- NULL
    return(ans)
}

#' @rdname optim
#' @export
p_nlm <- function(par0, objective, ...){
    npar <- length(par0)
    ans <- try(nlm(start = par0, resfn =objective, ..., iterlim=1000, print.level=0), silent=TRUE)
    # ans <- try(nlm(f=objective, p=par0, ..., iterlim=1000, print.level=0), silent=TRUE)
    # ans <- nlm(f_goal, par0, y = y, t = t, fun = doubleLog.gu)
    if (class(ans)[1] != "try-error") {
        ans$convcode <- ans$code
        if (ans$convcode == 1 || ans$convcode == 2 || ans$convcode == 3)
          ans$convcode <- 0
        if (ans$convcode == 4)
          ans$convcode <- 1
        # Translate output to common format
        ans$value <- ans$minimum
        ans$par <- setNames(ans$estimate, names(par0))
        ans$estimate <- NULL
        ans$minimum <- NULL
        ans$fevals <- NA
        ans$gevals <- NA  # ?? need to fix this somehow in nlm code
        ans$nitns <- ans$iterations
        ans$iterations <- NULL
    } else {
        # if (ctrl$trace > 0)
          # cat("nlm failed for this problem\n")
        ans <- list(fevals = NA)  # ans not yet defined, so set as list
        ans$value = NA#ctrl$badval
        ans$par <- rep(NA, npar)
        ans$convcode <- 9999  # failed in run
        ans$gevals <- NA
        ans$nitns <- NA
    }
    return(ans)
}

# c(p_nlminb, p_ncminf, p_nlm, p_optim)
# methods = c("Nelder-Mead", "BFGS", "L-BFGS-B", "CG", "SANN")

#' @rdname optim
#' @export
p_optim <- function(par0, objective, method = "BFGS", ...){
    npar <- length(par0)
    ans <- try(optim(par=par0, fn=objective, method = method, ...,
        control = list(maxit = 1000, trace = 0)), silent=TRUE)
    # The time is the index=1 element of the system.time for the process, which is a 'try()' of the regular optim() function
    if (class(ans)[1] != "try-error") {
      ans$convcode <- ans$convergence
      ans$convergence <- NULL
      # convergence: An integer code. '0' indicates successful convergence.  if (meth=='SANN') ans$convcode = 1 # always the case for
      # SANN (but it reports 0!)
      ans$fevals <- ans$counts[1]  # save function and gradient count information
      ans$gevals <- ans$counts[2]
      ans$counts <- NULL  # and erase the counts element now data is saved
    } else {
      # bad result -- What to do?
      ans <- list(fevals = NA)  # ans not yet defined, so set as list
      ans$convcode <- 9999  # failed in run
      # if (ctrl$trace > 0)
        cat("optim function evaluation failure\n")
      ans$value = NA#ctrl$badval
      ans$par <- rep(NA, npar)
      ans$fevals <- NA  # save function and gradient count information
      ans$gevals <- NA
    }
    ans$nitns <- NA  # not used
    return(ans)
}
