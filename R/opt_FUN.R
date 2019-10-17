#' @name opt_FUN
#' @title Unified optimization function
#'
#' @description
#' `I_optimx` is rich of functionality, but with a low computing
#' performance. Some basic optimization functions are unified here, with some
#' input and output format. \cr
#'
#' * `opt_ncminf` General-Purpose Unconstrained Non-Linear Optimization, 
#' see [ucminf::ucminf()].
#' * `opt_nlminb` Optimization using PORT routines, see [stats::nlminb()].
#' * `opt_nlm` Non-Linear Minimization, [stats::nlm()].
#' * `opt_optim` General-purpose Optimization, see [stats::optim()].
#'
#' @param par0 Initial values for the parameters to be optimized over.
#' @param objective A function to be minimized (or maximized), with first
#' argument the vector of parameters over which minimization is to take place.
#' It should return a scalar result.
#' @param method optimization method to be used in `p_optim`. See
#' [stats::optim()].
#' @param ... other parameters passed to `objective`.
#'
#' @return
#' * `convcode`: An integer code. 0 indicates successful convergence.
#' Various methods may or may not return sufficient information to allow all
#' the codes to be specified. An incomplete list of codes includes
#' 
#'   * `1`: indicates that the iteration limit `maxit` had been reached.
#'   * `20`: indicates that the initial set of parameters is inadmissible, 
#'    that is, that the function cannot be computed or returns an infinite, 
#'    NULL, or NA value.
#'   * `21`: indicates that an intermediate set of parameters is inadmissible.
#'   * `10`: indicates degeneracy of the Nelder--Mead simplex.
#'   * `51`: indicates a warning from the `"L-BFGS-B"` method; see component 
#'    `message` for further details.
#'   * `52`: indicates an error from the `"L-BFGS-B"` method; see component 
#'    `message` for further details.
#'   * `9999`: error
#' 
#' * `value`: The value of fn corresponding to par 
#' * `par`: The best parameter found 
#' * `nitns`: the number of iterations 
#' * `fevals`: The number of calls to `objective`.
#' 
#' @seealso [optim_pheno()], [I_optim()]
#' @example man/examples/ex-opt_FUN.R
NULL

# ' * `gevals`: The number of calls to gradient function. This excludes
# ' those calls needed to compute the Hessian, if requested, and any calls to
# ' \code{objective} to compute a finite-difference approximation to the
# ' gradient. }
opt_vars <- c("convcode", "value", "par", "nitns", "fevals")

#' @rdname opt_FUN
#' @importFrom ucminf ucminf
#' @export
opt_ucminf <- function(par0, objective, ...){
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
      ans$value <- NA#ctrl$badval
      ans$par      <- rep(NA, npar)
      ans$convcode <- 9999  # failed in run
      ans$gevals   <- NA
      ans$nitns    <- NA
    }
    # ans$convergence <- NULL
    ans[opt_vars]
}


#' @rdname opt_FUN
#' @export
opt_nlm <- function(par0, objective, ...){
    npar <- length(par0)

    # ans <- try(nlm(start = par0, resfn =objective, ..., iterlim=1000, print.level=0), silent=TRUE)
    ans <- try(nlm(f=objective, p=par0, ..., iterlim=1000, print.level=0), silent=TRUE)
      # , silent=TRUE)
    # ans <- nlm(f_goal, par0, y = y, t = t, fun = doubleLog.Gu)
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
        ans$gradient <- NULL
        ans$code <- NULL
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
    ans[opt_vars]
}

# c(p_nlminb, p_ncminf, p_nlm, p_optim)
# methods = c("Nelder-Mead", "BFGS", "L-BFGS-B", "CG", "SANN")

#' @rdname opt_FUN
#' @export
opt_optim <- function(par0, objective, method = "BFGS", ...){
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
    ans[opt_vars]
}
