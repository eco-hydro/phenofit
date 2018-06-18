#' Fine Curve fitting
#'
#' Fine Curve fitting for INPUT time-series. Procedures of initial weight, growing
#' season dividing and curve fitting are separated.
#'
#' @param extend_month For every growing season, previous and afterwards
#' `extend_month` are added to fitting daily time-series. In order to
#' @param minPercValid If valid percentage is less than \code{minPercValid}, the
#' fits are set to \code{list()}.
#'
#' @examples
#' INPUT <- check_input(d$date, d$EVI_500m, d$w, trim = T, maxgap = nptperyear / 4, alpha = 0.02)
#' # The detailed information of those parameters can be seen in `season`.
#' brks  <- season(INPUT, lambda, nptperyear, iters = 3, wFUN = wFUN, IsPlot = F,
#'                 south = south,
#'                 Aymin_less = 0.7,
#'                 max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5, ...) #, ...
#'
#' fit <- curvefit_site(INPUT, brks, lambda =lambda, IsPlot = T,
#'                      methods = c("AG", "zhang", "beck", "elmore", 'Gu'), #,"klos"
#'                      nptperyear = nptperyear, debug = F, wFUN = wTSM,
#'                      ymax_min = ymax_min,
#'                      south = d$lat[1] < 0)
#' plot_phenofit(fit, d) # plot to check the curve fitting
#' @export
curvefits <- function(INPUT, brks, nptperyear = 23,
                      wFUN = wTSM, iters = 2, wmin = 0.1,
                      south = FALSE,
                      extend_month = 2, minT = 0,
                      methods = c('AG', 'zhang', 'beck', 'elmore', 'Gu'),
                      qc, minPercValid = 0.3,
                      debug = FALSE, ...)
{
    t    <- INPUT$t
    n    <- length(t)
    doys <- as.numeric(difftime(t, t[1], units = "day") + 1)#days from origin

    # Tn for background module
    Tn     <- INPUT$Tn #if has no Tn, NULL will be return
    has_Tn <- ifelse(is_empty(Tn), F, T)
    y0     <- INPUT$y0 #original y
    if (is.null(y0)) y0 <- INPUT$y

    # title(x$site[1])
    if (all(is.na(INPUT$y))) return(NULL)
    # also constrained in `optim_pheno` function
    # if (sum(INPUT$w == 0)/length(INPUT$w) > 0.5) return(NULL) #much rigorous than all is.na
    w   <- INPUT$w
    I_w <- match(brks$whit$t, t) %>% rm_empty()
    w[I_w] <- brks$whit$w

    di <- brks$di
    if (is.null(di)){
        getDateId <- function(dates) match(dates, t) #%>% rm_empty()
        di <- data.table( beg  = getDateId(brks$dt$beg),
                          peak = getDateId(brks$dt$peak),
                          end  = getDateId(brks$dt$end)) %>% na.omit()
    }
    # possible snow or cloud, replaced with whittaker smooth.
    # I_fix <- which(w == wmin)
    # INPUT$y[I_fix] <- brks$whit %>% {.[[ncol(.)]][I_fix]}
    # w[I_fix]       <- 0.2 # exert the function of whitaker smoother

    # plot(y, type = "b"); grid()
    # lines(brks$whit$iter3, col = "blue")
    # lines(INPUT$y        , col = "red")

    ## 1. Curve fitting
    fits <- list()
    for (i in 1:nrow(di)){ #
        runningId(i)
        I    <- di$beg[i]:di$end[i]

        # extend curve fitting period
        period <- floor(nptperyear/12*extend_month)
        I_beg2 <- max(1, di$beg[i] - period)
        I_end2 <- min(n, di$end[i] + period)
        I_extend <- I_beg2:I_end2

        ti   <- doys[I_extend]
        yi   <- INPUT$y[I_extend]
        wi   <- w[I_extend]
        # original weights, put in w0 incurvefitting is unwisdom, but for plot
        w0   <- qc[I_extend] #INPUT$w
        # add background module here, 20180513
        if (has_Tn){
            Tni        <- Tn[I_extend]
            back_value <- backval(yi, ti, wi, Tni, minT, nptperyear)
            if (!is.na(back_value)){
                I_back     <- yi < back_value
                yi[I_back] <- back_value
                wi[I_back] <- 0.5
            }
        }
        beginI = ifelse(i == 1, 1, 2) # make sure no overlap
        tout <- doys[I] %>% {.[beginI]:last(.)} # make sure return the same length result.

        fit  <- curvefit(yi, ti, tout, nptperyear = nptperyear,
                         w = wi, w0 = w0, ylu = INPUT$ylu, iters = iters,
                         methods = methods, meth = 'BFGS', wFUN = wFUN, ...)
        # add original input data here, global calculation can comment this line
        data <- data.table(y = y0[I], t = doys[I], w = qc[I]) #INPUT$w[I]
        for (j in seq_along(fit)) fit[[j]]$data <- data
        # x <- fit$ELMORE
        # plot(y~t, x$data, type = "b"); grid()
        # lines(x$tout, x$fits$iter2)

        #if too much missing values
        if (sum(wi > pmax(wmin, 0.2))/length(wi) < minPercValid){
            fit %<>% map(function(x){
                x$fits %<>% map(~.x*NA) # list()
                return(x)
            })
        }
        fits[[i]] <- fit
    }
    # L1:curve fitting method, L2:yearly flag
    fits %<>% set_names(brks$dt$flag) %>% purrr::transpose()
    return(list(tout = t[first(di$beg):last(di$end)],  #dates for OUTPUT curve fitting VI
                fits = fits))
}


#' curvefits2
#'
#' @param t A date vector
#' @param y A numberic vector, same length as t
#' @param methods A character vector, c('spline', 'beck', 'elmore', 'klos', 'AG', 'zhang')
#' @param ... other parameters pass to season or curvefit_site
#' @export
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


phenonames <- c('TRS2.SOS', 'TRS2.EOS', 'TRS5.SOS', 'TRS5.EOS', 'TRS6.SOS', 'TRS6.EOS',
    'DER.SOS', 'DER.POP', 'DER.EOS',
    'UD', 'SD', 'DD','RD',
    'GreenUp', 'Maturity', 'Senescence', 'Dormancy')


#'
#' curve fit vegetation index (VI) time-series
#'
#' curve fit VI using 'spline', 'beck', 'elmore', 'klos', 'AG', 'Gu', 'zhang' methods
#'
#' @param y Vegetation time-series index, numeric vector
#' @param t The corresponding doy of x
#' @param tout The output interpolated time.
#'
#' @export
curvefit <- function(y, t = index(y), tout = t, meth = 'BFGS',
    methods = c('spline', 'beck', 'elmore', 'klos', 'AG', 'zhang'), ...)
{
    if (all(is.na(y))) return(NULL)
    if (length(methods) == 1 && methods == 'all')
        methods <- c('spline', 'beck', 'elmore', 'klos', 'AG', 'Gu', 'zhang')

    params <- list(y, t, tout, optimFUN = I_optim, method = meth, ...)

    if ('spline' %in% methods) fit.spline <- splinefit(y, t, tout)

    if ('beck'   %in% methods) fit.beck   <- do.call(FitDL.Beck,  c(params, pfun = p_nlminb))  #nlminb
    if ('elmore' %in% methods) fit.elmore <- do.call(FitDL.Elmore,c(params, pfun = p_nlminb))  #nlminb

    # best: BFGS, but its speed lower than other function, i.e. nlm
    if ('Gu'     %in% methods) fit.Gu     <- do.call(FitDL.Gu,    c(params, pfun = p_nlminb))  #nlm, ucminf
    if ('zhang'  %in% methods) fit.zhang  <- do.call(FitDL.Zhang, c(params, pfun = p_nlminb))  #nlm
    if ('AG'     %in% methods) fit.AG     <- do.call(FitAG,       c(params, pfun = p_nlminb))     #nlm
    if ('klos'   %in% methods) fit.klos   <- do.call(FitDL.Klos,  c(params, pfun = p_optim))   #BFGS, Nelder-Mead, L-BFGS-B

    names <- ls(pattern = "fit\\.") %>% set_names(., .)
    fits  <- lapply(names, get, envir = environment()) %>%
        set_names(toupper(gsub("fit\\.", "", names))) #remove `fit.` and update names
    return(fits)
}


#' curvefit_optimx
#'
#' With the help of `optimx` package, try to find which optimization function
#' is best.
#' @param optimFUN `I_optim` or `I_optimx`
#' @param meth c('BFGS','CG','Nelder-Mead','L-BFGS-B','nlm','nlminb',
#' spg','ucminf','Rcgmin','Rvmmin','newuoa','bobyqa','nmkb','hjkb')
#' @param pfun c(p_nlminb, p_ncminf, p_nlm, p_optim)
#'
#' @export
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
    if ('spline' %in% methods) fit.spline <- splinefit(x, t, tout)
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


#'
#' Get parameters from curve fitting result
#'
#' @param fit Object returned by curvefits
#' @export
getparam <- function(fit){
    llply(fit$fits, function(x){
        ldply(x, . %>% .$par, .id = "flag") %>% as_tibble()
    })
}

#'
#' Get parameters from multiple curve fitting results
#'
#' @param fits List Object of curvefits returned
#' @export
getparams <- function(fits){
    pars <- map(fits, getparams) %>% purrr::transpose() %>%
        map(~melt_list(.x, "site") %>% as_tibble())
    return(pars)
}
