# ' devirs of double logistics functions (analytical resolution)
# '
# ' Analytical derivative can also have some weird result. `numDeriv`'s result
# ' was accurate enough.
# '
# ' Enclosure env of `ans` was this funcion's execution env. But environment of
# ' C language based deriv function was wired.
# '
# ' grad is the first order derivative, and hess is the second order derivative.
# ' `f` means function, generated from `deriv` function
# ' `e` means expression, generated from `D` or `D` function
# ' 
# ' gradf_t   : gradient function f(par, t) (from the aspect of t)
# ' hessf_t   : hessian  function f(par, t) (from the aspect of t)
# ' grade_t   : gradient expression function f(par, t), get values through `eval`, (from the aspect of t)
# ' hesse_t   : gradient expression function f(par, t), get values through `eval`, (from the aspect of t)
# '
# ' grad_fpar : gradient function f(par, t), get values through `eval`, (from the aspect of t)
# ' hess_fpar : hessian  function f(par, t), get values through `eval`, (from the aspect of t)
# ' @param FUN Curve fitting function.
# ' @rdname curvefit_deriv
# ' @examples
# ' \dontrun{
# ' # grade_t' only cost 3/4 time of 'gradf_t'
# ' # hesse_t' used more 1/2 times of 'hessf_t'
# ' microbenchmark(
# '     gradf_t(FUN)(par,t),
# '     grade_t(FUN)(par,t))
# ' microbenchmark(
# '  hessf_t(FUN)(par,t),
# '  hesse_t(FUN)(par,t))
# ' }
# ' @export
grad_ft <- function(FUN){
    formula  <- attr(FUN, 'formula')
    parnames <- c("t", attr(FUN, 'par'))
    if (length(formula) == 2){
        # piece wise function
        f1 <- deriv(formula[1], "t", func = TRUE, function.arg = parnames)
        f2 <- deriv(formula[2], "t", func = TRUE, function.arg = parnames)
        # environment(f1) <- environment() #assign enclosure env
        # environment(f2) <- environment()
        # ans <- function(t) f1(t) * (t <= t0) + f2(t) * (t > t0)
        ans <- function(t, t0,  ...){
            values <- NA
            grad <- array(c(
                attr(f1(t[t <= t0], t0, ...), 'gradient'),
                attr(f2(t[t  > t0], t0, ...), 'gradient')), dim = c(length(t), 1))
            attr(values, 'gradient') <- grad
            return(values)
        }
    }else if (length(formula) == 1){
        ans <- deriv(formula[1], "t", func = TRUE, function.arg = parnames)
    }
    # environment(ans) <- environment()
    return(ans)
}

# ' @rdname curvefit_deriv
# ' @export
hess_ft <- function(FUN){
    formula <- attr(FUN, 'formula')
    parnames <- c("t", attr(FUN, 'par'))
    if (length(formula) == 2){
        f1 <- deriv3(formula[1], "t", func = TRUE, function.arg = parnames)
        f2 <- deriv3(formula[2], "t", func = TRUE, function.arg = parnames)
        # environment(f1) <- environment() #assign enclosure env
        # environment(f2) <- environment()
        # ans <- function(t) f1(t) * (t <= t0) + f2(t) * (t > t0)
        ans <- function(t, t0, ...){
            values <- NA
            hess <- array(c(
                attr(f1(t[t <= t0], t0, ...), 'hessian'),
                attr(f2(t[t  > t0], t0, ...), 'hessian')), dim = c(length(t), 1, 1))
            attr(values, 'hessian') <- hess
            return(values)
        }
    }else if (length(formula) == 1){
        ans <- deriv3(formula[1], "t", func = TRUE, function.arg = parnames)
    }
    # environment(ans) <- environment()
    return(ans)
}
# hess_ft = . %>% {deriv3(attr(., 'formula'), 't', func = TRUE)}

grad_fpar = . %>% {
    deriv (attr(., 'formula'), attr(., 'par'), func = TRUE)}

hess_fpar = . %>% {
    deriv3(attr(., 'formula'), attr(., 'par'), func = TRUE)}

# ' Higher derivatives
# '
# ' @param expr expression to be derivated.
# ' @param name Character, which variable to calculate derivative.
# ' @param order Integer, derivative order. 
# ' @return expression was returned.
# ' @rdname curvefit_deriv
DD <- function(expr, name, order = 1) {
   if(order < 1) stop("'order' must be >= 1")
   if(order == 1) D(expr, name)
   else DD(D(expr, name), name, order - 1)
}


# ' DERIVS from the aspect of par
# '
# ' For par aspect DERIVS, it had better to use deriv return function.
# ' For t aspect DERIVS, it had better to use D, then eval D expression
# '
# ' @rdname curvefit_deriv
gradf_t <- function(FUN) {
    # f <- deriv(attr(FUN, 'formula'), 't', func = TRUE) #attr(FUN, 'par')
    f <- grad_ft(FUN) #binding env, look variables in enclosure env
    # f is weird, it's envinment is global env

    function(par, t){
        ans <- do.call(f, c(list(t = t), as.list(par)))
        return(attr(ans, 'gradient'))
    }
}
# fun <- doubleLog.Zhang
# par <- setNames(1:7, attr(fun,"par"))
# fun_grad <- attr(fun, "grad")

# using `do.call` to avoid the difficult of environment finding
# ' @rdname curvefit_deriv
hessf_t <- function(FUN) {
    f <- hess_ft(FUN)
    # f <- deriv3(attr(FUN, 'formula'), 't', func = TRUE) #attr(FUN, 'par')
    function(par, t){
        ans <- do.call(f, c(list(t = t), as.list(par)))
        return(attr(ans, 'hessian'))
    }
}


## 02 ---------- DERIVS from the aspect of par --------------
# ' grad_par
# '
# ' @return return the gradient function of double logistics
# ' @rdname curvefit_deriv
grad_par <- function(FUN) {
    f <- grad_fpar(FUN)
    environment(f) <- environment()
    # print(environment(f))
    # whether eval and expression will be much suitable?

    function(par, t){
        e <- list2env(list(t = t), envir = environment(f))
        # 1. In this way, list will copy into a new env, cost spare time
        ans <- do.call(f, as.list(par)) #find t in enclosure env
        return(attr(ans, 'gradient'))
    }
}

hess_par <- function(FUN) {
    f <- hess_fpar(FUN) #attr(FUN, 'par')
    environment(f) <- environment()
    function(par, t){
        e <- list2env(list(t = t), envir = environment(f))
        ans <- do.call(f, as.list(par)) #find t in enclosure env
        return(attr(ans, 'hessian'))
    }
}

#' add hessian to this double logistic functions
# microbenchmark::microbenchmark({
    # DLs <- lapply(c('doubleLog.Beck', 'doubleLog.Elmore', 'doubleLog.Gu', 'doubleLog.Klos', 'doubleLog.Zhang', 'doubleAG'),
    #               function (FUN){
    #                 fun <- match.fun(FUN)
    #                 attr(fun, 'gradient') <- gradf_t(fun)   # gradient
    #                 attr(fun, 'hessian')  <- hessf_t(fun) # hessian
    #                 assign(FUN, fun, envir = globalenv())
    #                 fun
    #               })
# })
