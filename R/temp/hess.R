#' @export
hessian.default <- 
function (func, x, method = "Richardson", method.args = list(), 
          ...) 
{
    # if (1 != length(func(x, ...))) 
        # stop("Richardson method for hessian assumes a scalar valued function.")
    if (method == "complex") {
        args <- list(eps = 1e-04, d = 0.1, zero.tol = sqrt(.Machine$double.eps/7e-07), 
                     r = 4, v = 2)
        args[names(method.args)] <- method.args
        fn <- function(x, ...) {
            grad(func = func, x = x, method = "complex", side = NULL, 
                 method.args = list(eps = .Machine$double.eps), 
                 ...)
        }
        return(jacobian(func = fn, x = x, method = "Richardson", 
                        side = NULL, method.args = args, ...))
    }
    else if (method != "Richardson") 
        stop("method not implemented.")
    args <- list(eps = 1e-04, d = 0.1, zero.tol = sqrt(.Machine$double.eps/7e-07), 
                 r = 4, v = 2, show.details = FALSE)
    args[names(method.args)] <- method.args
    D <- genD(func, x, method = method, method.args = args, ...)$D
    if (1 != nrow(D)) 
        stop("BUG! should not get here.")
    H <- diag(NA, length(x))
    u <- length(x)
    for (i in 1:length(x)) for (j in 1:i) {
        u <- u + 1
        H[i, j] <- D[, u]
    }
    H <- H + t(H)
    diag(H) <- diag(H)/2
    H
}

environment(hessian.default) <- environment(numDeriv:::hessian.default)