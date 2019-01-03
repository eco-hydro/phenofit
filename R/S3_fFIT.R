#' @name fFIT
#' @title S3 class of fine curve fitting object.
#' 
#' @description
#' \code{fFIT} is returned by \code{\link{optim_pheno}}.
#' 
#' @format
#' \describe{
#'   \item{tout}{Corresponding doy of prediction}
#'   \item{zs}{curve fitting values of every iteration}
#'   \item{ws}{weight of every iteration}
#'   \item{par}{Optimized parameter of fine curve fitting method}
#'   \item{fun}{The name of fine curve fitting function.}
#' }
NULL

#' @export
print.fFIT <- function(x, ...){
    FUN <- get(x$fun, mode = 'function')
    cat(sprintf(" formula:\t%s\n", attr(FUN, "formula") ))
    cat("pars:\n")
    print(x$par)
}
