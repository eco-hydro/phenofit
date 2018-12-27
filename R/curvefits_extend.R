phenonames <- c('TRS2.SOS', 'TRS2.EOS', 'TRS5.SOS', 'TRS5.EOS', 'TRS6.SOS', 'TRS6.EOS',
    'DER.SOS', 'DER.POP', 'DER.EOS',
    'UD', 'SD', 'DD','RD',
    'GreenUp', 'Maturity', 'Senescence', 'Dormancy')



#' curve fit vegetation index (VI) time-series
#'
#' curve fit VI using  methods
#'
#' @param y Vegetation time-series index, numeric vector
#' @param t The corresponding doy of x
#' @param tout The output interpolated time.
#' @param meth optimization method
#' @param methods Fine curve fitting methods, can be one or more of \code{c('AG', 
#' 'Beck', 'Elmore', 'Gu', 'Klos', 'Zhang')}. 
#' @param ... other parameters passed to curve fitting function.
#' 
#' @note 'Klos' and 'Gu' have many parameters. It will be slow and not stable.
#' 
#' @return fits
#' @seealso \code{\link{FitAG}}, \code{\link{FitDL.Beck}}, 
#' \code{\link{FitDL.Elmore}}, \code{\link{FitDL.Gu}}, 
#' \code{\link{FitDL.Klos}}, \code{\link{FitDL.Zhang}}
#' 
#' @export
curvefit <- function(y, t = index(y), tout = t, meth = 'BFGS',
    methods = c('AG', 'Beck', 'Elmore', 'Gu', 'Klos', 'Zhang'), ...)
{
    if (all(is.na(y))) return(NULL)
    if (length(methods) == 1 && methods == 'all')
        methods <- c('AG', 'Beck', 'Elmore', 'Gu', 'Klos', 'Zhang')

    params <- list(y, t, tout, optimFUN = I_optim, method = meth, ...)

    # if ('spline' %in% methods) fit.spline <- splinefit(y, t, tout)

    if ('AG'     %in% methods) fit.AG     <- do.call(FitAG,       c(params, pfun = p_nlminb))     #nlm
    if ('Beck'   %in% methods) fit.Beck   <- do.call(FitDL.Beck,  c(params, pfun = p_nlminb))  #nlminb
    if ('Elmore' %in% methods) fit.Elmore <- do.call(FitDL.Elmore,c(params, pfun = p_nlminb))  #nlminb

    # best: BFGS, but its speed lower than other function, i.e. nlm
    if ('Gu'     %in% methods) fit.Gu     <- do.call(FitDL.Gu,    c(params, pfun = p_nlminb))  #nlm, ucminf
    if ('Klos'   %in% methods) fit.Klos   <- do.call(FitDL.Klos,  c(params, pfun = p_optim))   #BFGS, Nelder-Mead, L-BFGS-B
    if ('Zhang'  %in% methods) fit.Zhang  <- do.call(FitDL.Zhang, c(params, pfun = p_nlminb))  #nlm

    names <- ls(pattern = "fit\\.") %>% set_names(., .)
    fits  <- lapply(names, get, envir = environment()) %>%
        set_names(toupper(gsub("fit\\.", "", names))) #remove `fit.` and update names
    return(fits)
}


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
    methods = c('spline', 'beck', 'elmore', 'AG'),
    meth, pfun, ...)
{
    if (all(is.na(x))) return(NULL)
    ##2. curve fitting
    if (length(methods) == 1 && methods == 'all')
        methods <- c('spline', 'beck', 'elmore', 'klos', 'AG', 'Gu', 'zhang')

    # failed: 'BFGS', 'Nelder-Mead', 'L-BFGS-B'
    # meth = 'L-BFGS-B'
    # ok: 'L-BFGS-B'; failed: 'BFGS', 'Nelder-Mead'
    params <- list(x, t, tout, optimFUN = optimFUN, pfun = pfun, method = meth, ...)
    # fit.beck   <- FitDL.Beck   #even Nelder-Mead was faster and convergent, but nlminb was better
    # ok: BFGS; failed: 'L-BFGS-B'
    # if ('spline' %in% methods) fit.spline <- splinefit(x, t, tout)
    if ('beck'   %in% methods) fit.beck   <- do.call(FitDL.Beck, params)        #nlminb
    if ('elmore' %in% methods) fit.elmore <- do.call(FitDL.Elmore, params)      #nlminb

    # best: BFGS, but its speed lower than other function, i.e. nlm
    if ('Gu'     %in% methods) fit.Gu     <- do.call(FitDL.Gu, params)          #nlm, ucminf
    if ('zhang'  %in% methods) fit.zhang  <- do.call(FitDL.Zhang, params)       #nlm
    if ('AG'     %in% methods) fit.AG     <- do.call(FitAG, params)             #nlm
    if ('klos'   %in% methods) fit.klos   <- do.call(FitDL.Klos, params)        #BFGS, Nelder-Mead, L-BFGS-B

    # test for optimx methods
    # fit   <- FitDL.Zhang  (x, t, tout, optimFUN = optimx_fun, debug = T, method = 'BFGS')
    names <- ls(pattern = "fit\\.") %>% set_names(., .)
    fits  <- lapply(names, get, envir = environment()) %>%
        set_names(toupper(gsub("fit\\.", "", names))) #remove `fit.` and update names
    return(fits)
}



#' Get parameters from curve fitting result
#'
#' @param fit Curve fitting result by \code{curvefits} result. 
#' @export
getparam <- function(fit){
    llply(fit$fits, function(x){
        ldply(x, . %>% .$par, .id = "flag") %>% as_tibble()
    })
}

#' @param fits Multiple methods curve fitting results by \code{curvefits} result. 
#' @rdname getparam
#' @export
getparams <- function(fits){
    pars <- map(fits, getparams) %>% purrr::transpose() %>%
        map(~melt_list(.x, "site") %>% as_tibble())
    return(pars)
}


##  curvefits2
## 
##  @param t A date vector
##  @param y A numberic vector, same length as t
##  @param methods A character vector, c('spline', 'beck', 'elmore', 'klos', 'AG', 'zhang')
##  @param ... other parameters pass to season or curvefit_site
##  @export
curvefits2 <- function(t, y, w, nptperyear = 46,
                          wFUN = wTSM, iters = 2,
                          lambda, south = FALSE,
                          IsPlot = FALSE,
                          Aymin_less = 0.6, ymax_min = 0.1,
                          methods = c('AG', 'zhang', 'beck', 'elmore', 'Gu'),
                          debug = FALSE, ...)
{
    # 1. Check input data and initial parameters for phenofit
    INPUT <- check_input(t, y, w, maxgap = nptperyear / 4, alpha = 0.02)

    # 2. The detailed information of those parameters can be seen in `season`.
    brks  <- season(INPUT, lambda, nptperyear, iters = 3, wFUN = wFUN,
                    IsPlot = IsPlot, south = south,
                    Aymin_less = Aymin_less, ymax_min = ymax_min,
                    max_MaxPeaksperyear=2.5, max_MinPeaksperyear=3.5, ...)
    ## 3.1. curve fitting
    fit <- curvefits(INPUT, brks, lambda =lambda,
                         methods = c("AG", "zhang", "beck", "elmore", 'Gu'), #,"klos"
                         nptperyear = nptperyear, debug = F, wFUN = wTSM, ...)
    ## 3.2 Get GOF information
    stat  <- ldply(fit$fits, function(fits_meth){
        ldply(fits_meth, statistic.phenofit, .id = "flag")
    }, .id = "meth")

    # 4. extract phenology
    # pheno: list(p_date, p_doy)
    p <- lapply(fit$fits, ExtractPheno)
    pheno  <- map(p, tidyFitPheno, origin = t[1]) %>% purrr::transpose()

    fit$INPUT   <- INPUT
    fit$seasons <- brks
    fit$stat    <- stat
    fit$pheno   <- pheno
    return(fit)
}