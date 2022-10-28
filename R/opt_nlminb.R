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
