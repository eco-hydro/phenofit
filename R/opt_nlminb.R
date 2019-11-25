#' @rdname opt_FUN
#' @export
opt_nlminb <- function(par0, objective, ...){
    npar <- length(par0)
    ans  <- nlminb(start=par0, objective=objective, ...,
      control = list(eval.max=1000, iter.max=1000, trace=0, abs.tol=0))
    # unify the format of optimization results
    # if (class(ans)[1] != "try-error") {
        ans$convcode    <- ans$convergence
        # Translate output to common format and names
        ans$value       <- ans$objective
        ans$objective   <- NULL
        ans$fevals      <- ans$evaluations[[1]]
        ans$gevals      <- ans$evaluations[[2]]
        ans$evaluations <- NULL  # cleanup
        ans$nitns       <- ans$iterations
        ans$iterations  <- NULL
    # } else {
    #     # bad result -- What to do?
    #     ans <- list(fevals = NA)  # ans not yet defined, so set as list
    #     ans$convcode <- 9999  # failed in run
    #     # if (ctrl$trace > 0)
    #       cat("nlminb function evaluation failure\n")
    #     ans$value  <- NA#ctrl$badval
    #     ans$par    <- rep(NA, npar)
    #     ans$nitns  <- NA  # not used
    #     ans$gevals <- NA  ## ?? missing 130929
    # }
    ans[opt_vars]
}

#' Optimization using PORT routines
#'
#' Unconstrained and box-constrained optimization using PORT routines.
#'
#' @inheritParams optim_pheno
#' @inheritParams stats::nlminb
#' @param par0 Initial values for the parameters to be optimized over.
#' @param fitMeth Curve fitting methods, one of `c("doubleLog_Beck",
#' "doubleLog_Elmore", "doubleLog_AG", "doubleLog_Zhang")`
#' @param ... ignored parameters
#'
#' @return A list object of
#' - `par`: The optimal parameters
#' - `convergence`:
#'   - 0: convergent;
#'   - 1: Non-convergent
#' - iterations
#' - evaluations: list(function, gradient)
#' - objective
#'
#' @seealso [stats::nlminb()]
#' @example man/examples/ex-nlminb_julia.R
#'
#' @export
opt_nlminb_julia <- function(
    par0,
    fitMeth = "doubleLog_Beck",
    y, t,
    w = NULL, ylu = NULL,
    lower = NULL, upper = NULL, ...)
{
    pred <- y*0
    methods = c("doubleLog_Beck", "doubleLog_Elmore", "doubleLog_AG", "doubleLog_Zhang", "doubleLog_Gu")
    fitMeth = methods[pmatch(fitMeth, methods)]

    julia_assign("par0" , par0)
    julia_assign("t"    , t)
    julia_assign("y"    , y)
    julia_assign("ypred", pred)
    julia_assign("w"    , w)
    # julia_assign("ylu"  , NULL) #check_double(ylu)
    julia_assign("lower", check_double(lower))
    julia_assign("upper", check_double(upper))

    # in the same order of goal function in julia
    # ylu,
    cmd <- sprintf("optim_nlminb(par0, goal!, %s!, y, t, ypred, w,
        lower = lower, upper = upper, verbose = true, eval_max = 1000, iter_max = 1000)", fitMeth)
    julia_eval(cmd)
}

check_double <- function(x) {
    if(!is.null(x)) as.double(x) else NULL
}

#' julia_init
#' 
#' @keywords internal
#' @import JuliaCall
#' @export
julia_init <- function(){
    julia <- julia_setup()
    julia_library("nlminb")
}
