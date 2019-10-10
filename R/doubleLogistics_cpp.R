doubleLogMain <- function(.NAME, par, t, pred) {
    miss_pred = missing(pred)
    if (miss_pred) pred = t*0

    invisible(.Call(.NAME, par, t, pred))
    if (miss_pred) return(pred)
}

#' Double logistics in Rcpp
#'
#' @inheritParams doubleLog.Beck
#' @param pred Numeric Vector, predicted values
#' 
#' @seealso [doubleLog.Beck()]
#'
#' @keywords internal
#' @export
logistic <- function(par, t, pred) {
    doubleLogMain(`_phenofit_clogistic`, par, t, pred)
}

#' @rdname logistic
#' @export
doubleLog_Zhang <- function(par, t, pred) {
    doubleLogMain(`_phenofit_cdoubleLog_Zhang`, par, t, pred)
}

#' @rdname logistic
#' @export
doubleLog_AG <- function(par, t, pred) {
    doubleLogMain(`_phenofit_cdoubleLog_AG`, par, t, pred)
}

#' @rdname logistic
#' @export
doubleLog_Beck <- function(par, t, pred) {
    doubleLogMain(`_phenofit_cdoubleLog_Beck`, par, t, pred)
}

#' @rdname logistic
#' @export
doubleLog_Elmore <- function(par, t, pred) {
    doubleLogMain(`_phenofit_cdoubleLog_Elmore`, par, t, pred)
}

#' @rdname logistic
#' @export
doubleLog_Gu <- function(par, t, pred) {
    doubleLogMain(`_phenofit_cdoubleLog_Gu`, par, t, pred)
}

#' @rdname logistic
#' @export
doubleLog_Klos <- function(par, t, pred) {
    doubleLogMain(`_phenofit_cdoubleLog_Klos`, par, t, pred)
}

#' @export
f_goal2 <- function(par, fun, y, t, ypred, w = NULL, ylu = NULL, ...) {
    .Call(`_phenofit_f_goal_cpp`, par, fun, y, t, ypred, w, ylu)
}

# set par and names for double Logistics functions
funcs = lsf.str(pattern = "^doubleLog_")

for (func in funcs) {
    funr = gsub("_", ".", func)
    eval(parse(text = sprintf("attr(%s, 'name') <- '%s'", func, func)))
    eval(parse(text = sprintf("attr(%s, 'name') <- '%s'", funr, funr)))
    eval(parse(text = sprintf("attr(%s, 'par')  <- attr(%s, 'par')", func, funr)))
    eval(parse(text = sprintf("attr(%s, 'formula')  <- attr(%s, 'formula')", func, funr)))
    eval(parse(text = sprintf("attr(%s, 'gradient') <- attr(%s, 'gradient')", func, funr)))
    eval(parse(text = sprintf("attr(%s, 'hessian')  <- attr(%s, 'hessian')", func, funr)))
}
