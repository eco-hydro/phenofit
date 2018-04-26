curvefit_site <- function(t, y, w, nptperyear = 46,
                          wFUN = wTSM, iters = 2,
                          lambda, south = FALSE,
                          IsPlot = FALSE,
                          methods = c('AG', 'zhang', 'beck', 'elmore', 'Gu'),
                          debug = FALSE, ...)
{
    doys       <- as.numeric(difftime(t, t[1], units = "day") + 1)#days from origin

    INPUT <- check_input(t, y, w, trim = T, maxgap = nptperyear / 4)
    brks  <- season(INPUT, lambda, nptperyear, iters = 3, wFUN = wTSM, IsPlot = F,
                    south = south,
                    max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5, ...) #, ...
    # title(x$site[1])
    # fit <- curvefit_site(x$date, x$GPP_NT, lambda =1e4,
    #                      methods = c("AG", "zhang", "beck", "elmore", 'Gu'), #,"klos"
    #                      nptperyear = 365, debug = F, wFUN = wTSM,
    #                      south = x$lat[1] < 0)
    if (all(is.na(INPUT$y))) return(NULL)
    # also constrained in `optim_pheno` function
    # if (sum(INPUT$w == 0)/length(INPUT$w) > 0.5) return(NULL) #much rigorous than all is.na

    # season(INPUT, lambda, nptperyear = 46, iters = 2, wFUN = wTSM, IsPlot = TRUE)
    # brks <- season(INPUT, lambda, nptperyear, south = south, iters = 3, wFUN = wTSM, IsPlot = FALSE)
    w  <- brks$whit$w
    di <- brks$di

    if (debug){
        fits <- stat <- pheno<- NULL
    }else{
        ## 1. Curve fitting
        fits <- list()
        for (i in 1:nrow(di)){ #
            runningId(i)
            I    <- di$beg[i]:di$end[i]
            # Try to remove NA values at head and tail. Fialed, Na value may not
            # at head or tail.
            # I_nona <- which(!is.na(y[I])) %>% {I[first(.):last(.)]}
            I_nona <- I#checking year brk is enough
            ti     <- doys[I_nona]

            # check year brokens
            I_brkyear <- which(diff(ti) >= 365)
            nbrk      <- length(I_brkyear)

            if (nbrk > 0 & nbrk <= 2){
                if (nbrk == 1) {
                    I_1 <- I_nona[1:I_brkyear]
                    I_2 <- I_nona[(I_brkyear + 1):length(I_nona)]

                    lst <- list(I_1, I_2)
                } else if (nbrk == 2) {
                    I_1 <- I_nona[1:I_brkyear[1]]
                    I_2 <- I_nona[(I_brkyear[1]+1):I_brkyear[2]]
                    I_3 <- I_nona[(I_brkyear[2]+1):length(I_nona)]

                    lst <- list(I_1, I_2, I_3)
                }
                #select the longest segment
                I_nona <- lst[[which.max(sapply(lst, length))]]
                ti   <- doys[I_nona]
            }

            yi   <- INPUT$y[I_nona]
            wi   <- w[I_nona]
            tout <- doys[I] %>% {first(.):last(.)} # make sure return the same length result.

            fit <- curvefit(yi, ti, tout, nptperyear = nptperyear,
                            w = wi, ylu = INPUT$ylu, iters = iters,
                            methods = methods, meth = 'BFGS', wFUN = wFUN, ...)
            #if too much missing values
            if (sum(is.na(yi))/length(I_nona) > 0.5){
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
    return(list(INPUT = tibble(t, y),
                seasons = brks,
                tout = t[first(di$beg):last(di$end)],  #dates for OUTPUT curve fitting VI
                fits = fits, stat = stat,
                pheno = pheno))
}