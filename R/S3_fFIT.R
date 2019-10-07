#' @name fFIT
#' @title S3 class of fine curve fitting object.
#' 
#' @description
#' `fFIT` is returned by [optim_pheno()].
#' 
#' @format
#' * `tout`: Corresponding doy of prediction
#' * `zs`: curve fitting values of every iteration
#' * `ws`: weight of every iteration
#' * `par`: Optimized parameter of fine curve fitting method
#' * `fun`: The name of fine curve fitting function.
#' 
#' @keywords internal
NULL

#' @export
print.fFIT <- function(x, ...){
    FUN <- get(x$fun, mode = 'function')
    cat(sprintf(" formula:\t%s\n", attr(FUN, "formula") ))
    cat("pars:\n")
    print(x$par)
}
