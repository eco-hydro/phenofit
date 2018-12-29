context("Analytical Derivs")

expect_silent({
    # attach gradient and hessian analytical function to curve fitting functions
    .dls <- lapply(
        c("doubleLog.Beck", "doubleLog.Elmore", "doubleLog.Gu",
          "doubleLog.Klos", "doubleLog.Zhang", "doubleAG"),
        function (FUN){
            # FUN <- deparse(substitute(fun))
            fun <- get(FUN)
            attr(fun, 'gradient') <- gradf_t(fun) # gradient
            attr(fun, 'hessian')  <- hessf_t(fun) # hessian
            # print(environment(fun))
            # print(fun)
            # print(FUN)
            # assign(FUN, fun, envir = environment(fun)) #environment("namespace:phenofit"))#
            fun
        })
})
