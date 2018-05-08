#' curvefits
#'
#' Curve fitting for INPUT time-series. Procedures of initial weight, growing
#' season dividing and curve fitting are separated.
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
curvefits <- function(INPUT, brks, nptperyear = 46,
                          wFUN = wTSM, iters = 2,
                          lambda, south = FALSE,
                          IsPlot = FALSE,
                          methods = c('AG', 'zhang', 'beck', 'elmore', 'Gu'),
                          debug = FALSE, ...)
{
    t    <- INPUT$t
    doys <- as.numeric(difftime(t, t[1], units = "day") + 1)#days from origin

    # title(x$site[1])
    if (all(is.na(INPUT$y))) return(NULL)
    # also constrained in `optim_pheno` function
    # if (sum(INPUT$w == 0)/length(INPUT$w) > 0.5) return(NULL) #much rigorous than all is.na

    w  <- brks$whit$w
    di <- brks$di

    # possible snow or cloud, replaced with whittaker smooth.
    I_fix <- which(w == 0)
    INPUT$y[I_fix] <- brks$whit %>% {.[[ncol(.)]][I_fix]}

    # plot(y, type = "b"); grid()
    # lines(brks$whit$iter3, col = "blue")
    # lines(INPUT$y        , col = "red")

    w[I_fix]       <- 0.2 #exert the function of whitaker smoother

    if (debug){
        fits <- stat <- pheno<- NULL
    }else{
        ## 1. Curve fitting
        fits <- list()
        for (i in 1:nrow(di)){ #
            runningId(i)
            I    <- di$beg[i]:di$end[i]
            ti   <- doys[I]
            yi   <- INPUT$y[I]
            wi   <- w[I]
            tout <- doys[I] %>% {first(.):last(.)} # make sure return the same length result.

            fit <- curvefit(yi, ti, tout, nptperyear = nptperyear,
                            w = wi, ylu = INPUT$ylu, iters = iters,
                            methods = methods, meth = 'BFGS', wFUN = wFUN, ...)
            #if too much missing values
            if (sum(wi > 0.2)/length(wi) < 0.3){
                fit %<>% map(function(x){
                    x$fits %<>% map(~.x*NA)
                    x$pred %<>% multiply_by(NA)
                    return(x)
                })
            }
            fits[[i]] <- fit
            # p <- getFits_pheno(fits[[i]], IsPlot = IsPlot)
        }
        # L1:curve fitting method, L2:yearly flag
        fits %<>% set_names(brks$dt$flag) %>% purrr::transpose()

        ## 2. GOF, try to give fitting informations for every curve fitting
        stat  <- ldply(fits, function(fits_meth){
            ldply(fits_meth, statistic.phenofit, .id = "flag")
        }, .id = "meth")

        # 3. phenology
        p <- lapply(fits, getFits_pheno)
        # pheno: list(p_date, p_doy)
        pheno  <- map(p, tidyFits_pheno, origin = t[1]) %>% purrr::transpose()
    }
    return(list(INPUT = tibble(t, y = INPUT$y),
                seasons = brks,
                tout = t[first(di$beg):last(di$end)],  #dates for OUTPUT curve fitting VI
                fits = fits, stat = stat,
                pheno = pheno))
}

#' 
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
    # Check input data and initial parameters for phenofit
    INPUT <- check_input(t, y, w, trim = T, maxgap = nptperyear / 4, alpha = 0.02)

    # The detailed information of those parameters can be seen in `season`. 
    brks  <- season(INPUT, lambda, nptperyear, iters = 3, wFUN = wFUN, IsPlot = T,
                    south = south,
                    Aymin_less = Aymin_less, ymax_min = ymax_min, 
                    max_MaxPeaksperyear=2.5, max_MinPeaksperyear=3.5, ...)

    fit <- curvefits(INPUT, brks, lambda =lambda, IsPlot = T, 
                         methods = c("AG", "zhang", "beck", "elmore", 'Gu'), #,"klos"
                         nptperyear = nptperyear, debug = F, wFUN = wTSM, ...)
    return(fit)
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
