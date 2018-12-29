## ' curvefit_optimx; still under debug
## '
## ' With the help of `optimx` package, try to find which optimization function
## ' is best.
## ' @param optimFUN `I_optim` or `I_optimx`
## ' @param meth c('BFGS','CG','Nelder-Mead','L-BFGS-B','nlm','nlminb',
## ' spg','ucminf','Rcgmin','Rvmmin','newuoa','bobyqa','nmkb','hjkb')
## ' @param pfun c(p_nlminb, p_ncminf, p_nlm, p_optim)
## '
## ' @export
curvefit_optimx <- function(x, t = index(x), tout = t,
    optimFUN = I_optimx,
    methods = c('spline', 'Beck', 'Elmore', 'AG'),
    meth, pfun, ...)
{
    if (all(is.na(x))) return(NULL)
    ##2. curve fitting
    if (length(methods) == 1 && methods == 'all')
        methods <- c('spline', 'AG', 'Beck', 'Elmore', 'Gu', 'Klos', 'Zhang')

    # failed: 'BFGS', 'Nelder-Mead', 'L-BFGS-B'
    # meth = 'L-BFGS-B'
    # ok: 'L-BFGS-B'; failed: 'BFGS', 'Nelder-Mead'
    params <- list(x, t, tout, optimFUN = optimFUN, pfun = pfun, method = meth, ...)
    # fit.beck   <- FitDL.Beck   #even Nelder-Mead was faster and convergent, but nlminb was better
    # ok: BFGS; failed: 'L-BFGS-B'
    # if ('spline' %in% methods) fit.spline <- splinefit(x, t, tout)
    if ('Beck'   %in% methods) fit.Beck   <- do.call(FitDL.Beck, params)        #nlminb
    if ('Elmore' %in% methods) fit.Elmore <- do.call(FitDL.Elmore, params)      #nlminb

    # best: BFGS, but its speed lower than other function, i.e. nlm
    if ('Gu'     %in% methods) fit.Gu     <- do.call(FitDL.Gu, params)          #nlm, ucminf
    if ('Zhang'  %in% methods) fit.Zhang  <- do.call(FitDL.Zhang, params)       #nlm
    if ('AG'     %in% methods) fit.AG     <- do.call(FitAG, params)             #nlm
    if ('Klos'   %in% methods) fit.Klos   <- do.call(FitDL.Klos, params)        #BFGS, Nelder-Mead, L-BFGS-B

    # test for optimx methods
    # fit   <- FitDL.Zhang  (x, t, tout, optimFUN = optimx_fun, debug = T, method = 'BFGS')
    names <- ls(pattern = "fit\\.") %>% set_names(., .)
    fits  <- lapply(names, get, envir = environment()) %>%
        set_names(toupper(gsub("fit\\.", "", names))) #remove `fit.` and update names
    return(fits)
}
