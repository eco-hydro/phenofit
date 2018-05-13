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
                      minT = 0,
                      methods = c('AG', 'zhang', 'beck', 'elmore', 'Gu'),
                      debug = FALSE, ...)
{
    t    <- INPUT$t
    n    <- length(t)
    doys <- as.numeric(difftime(t, t[1], units = "day") + 1)#days from origin

    # Tn for background module
    Tn <- INPUT$Tn #if has no Tn, NULL will be return
    has_Tn <- ifelse(is_empty(Tn), F, T)

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

            # extend curve fitting period
            period <- floor(nptperyear/12*3)
            I_beg2 <- max(1, di$beg[i] - period)
            I_end2 <- min(n, di$end[i] + period)
            I_extend <- I_beg2:I_end2

            ti   <- doys[I_extend]
            yi   <- INPUT$y[I_extend]
            wi   <- w[I_extend]

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
        }
        # L1:curve fitting method, L2:yearly flag
        fits %<>% set_names(brks$dt$flag) %>% purrr::transpose()
    }
    return(list(tout = t[first(di$beg):last(di$end)],  #dates for OUTPUT curve fitting VI
                fits = fits))
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
    # 1. Check input data and initial parameters for phenofit
    INPUT <- check_input(t, y, w, trim = T, maxgap = nptperyear / 4, alpha = 0.02)

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

#' ExtractPheno
#'
#' Get yearly vegetation phenological metrics of a curve fitting method
#'
#' @param outdir save calculated phenology in oudir txts
#' @export
ExtractPheno <- function(fits, TRS = c(0.1, 0.2, 0.5, 0.6), IsPlot = FALSE){
    names <- names(fits)
    pheno_list <- list()
    methods    <- c(paste0("TRS", TRS*10),"DER","GU", "ZHANG")
    TRS_last   <- last(TRS) # only display last threshold figure

    if (IsPlot)
        op <- par(mfrow = c(length(fits), 5),
            oma = c(1, 2, 3, 1), mar = rep(0, 4), yaxt = "n", xaxt = "n")
    for (i in seq_along(fits)){
        fit    <- fits[[i]]
        ypred  <- last(fit$fits)
        all_na <- all(is.na(ypred))
        # 1. show curve fitting RMSE
        # need to fix here, about status variable. 31 Jan, 2018
        if (IsPlot && !all_na){
            PhenoPlot(fit$tout, ypred)
            lines(fit$data$t, fit$data$x, type = "b")
            # legend(min(fit$data$t) - 40, max(fit$pred),
            #     sprintf("rmse = %.2e\nnash = %.2f\nR=%.2f, p = %.3f",
            #         status$rmse, status$nash, status$r, status$pvalue),
            #     adj = c(0, 0), bty='n', text.col = "red")
            mtext(names[i], side = 2)
        }
        if (i == 1 && IsPlot) mtext("Fitting")

        p_TRS <- map(TRS, function(trs) {
            PhenoTrs(fit, approach = "White", trs = trs, IsPlot = FALSE)
        })

        if (IsPlot && !all_na) {
            trs6 <- PhenoTrs(fit, approach="White", trs=TRS_last, IsPlot = IsPlot)
            if (i == 1) mtext(sprintf('TRS%d', TRS_last*10))
        }

        der   <- PhenoDeriv(fit, IsPlot);    if (i == 1 && IsPlot) mtext("DER")
        gu    <- PhenoGu(fit, IsPlot)[1:4];  if (i == 1 && IsPlot) mtext("GU")
        zhang <- PhenoKl(fit, IsPlot);       if (i == 1 && IsPlot) mtext("ZHANG")
        pheno_list[[i]] <- c(p_TRS, list(der, gu, zhang)) %>% set_names(methods)
    }
    pheno_list %<>% set_names(names)
    return(pheno_list)
}
