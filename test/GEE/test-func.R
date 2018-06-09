#' @param FUN Curve fitting functions, call be one of `sgfitw`, `whitsmw2` and `wHANTS`.
season <- function(INPUT, nptperyear = 46, south = FALSE,
                   FUN = whitsmw2, wFUN = wTSM, iters = 2,
                   lambda, nf  = 2, frame = floor(nptperyear/7) * 2 + 1,
                   minpeakdistance = nptperyear/6,
                   ymax_min = 0.5,
                   rymin_less = 0.6, #ymin < rymin_less * A
                   threshold_max = 0.2, threshold_min = 0.05,
                   IsPlot = FALSE,
                   max_MaxPeaksperyear = 2, max_MinPeaksperyear = 3,
                   wmin = 0.1,
                   plotdat = INPUT, print = FALSE,
                   ...)
{
    t   <- INPUT$t
    y   <- INPUT$y

    if (all(is.na(y))) return(NULL)
    n     <- length(y)
    npt   <- sum(INPUT$w > wmin)
    nyear <- ceiling(npt/nptperyear)

    frame <- floor(nptperyear/7) * 2 + 1 #13, reference by TSM SG filter
    if (missing(lambda)) lambda <- max(nyear*frame, 15)

    ## 1. weighted curve fitting help to divide growing season
    iloop <- 1
    while (iloop <= 3){
        # sgfitw(INPUT$y, INPUT$w, nptperyear, INPUT$ylu, wFUN, iters, frame, d=2)
        # whitsmw(y, w, ylu, wFUN, iters = 1, lambda = 100, ..., d = 2, missval)
        param <- c(INPUT, nptperyear = nptperyear,
            wFUN = wFUN, wmin = wmin, iters = iters,
            lambda = lambda,  # param for whittaker
            nf     = nf,      # param for HANTS,
            frame  = frame    # param for
        )
        yfits <- do.call(FUN, param)
        ypred <- last(yfits) #as.numeric(runmed(ypred, frame))
        alpha <- 0.01

        # default is three year data input, median will be much better
        ylu_min <- aggregate(ypred, list(year = year(t)), min)$x %>% median()
        ylu_max <- aggregate(ypred, list(year = year(t)), max)$x %>% median()
        ylu <- c(pmax(ylu_min, INPUT$ylu[1]), #quantile(ypred, alpha/2)
                 pmin(ylu_max, INPUT$ylu[2]))
        # ylu   <- quantile(ypred, c(alpha/2, 1 - alpha), na.rm = T)
        A         <- diff(ylu)
        ## Plateau peak will lead to failed to find local extreme values.
        #  To avoid fluctuating in peak of growing season or flat peak growing
        #  season, like fluxsite: ZM-Mon
        # max_slp <- 2*A/nptperyear
        # pek_slp <- abs(coefficients(lm(ypred[I]~I))[[2]])
        I <- which(ypred > (0.8*A + ylu[1]))
        if (length(I)/length(y) > 0.3){
            ypred[I] <- median(ypred[I])
        }
        # local minimum values
        # peak values is small for minimum values, so can't use threshold_min here
        peaks <- findpeaks(-ypred,
                           threshold_max = threshold_max*A,
                           threshold_min = threshold_min*0,
                           minpeakdistance = minpeakdistance, zero = "-", nups = 0)
        pos_min   <- peaks$X
        pos_min[, 1] %<>% multiply_by(-1)
        npeak_MinPerYear <- length(peaks$gregexpr)/nyear#max peaks
        # local maximum values,
        peaks <- findpeaks(ypred, zero = "+",
                           threshold_max = threshold_max*A,
                           threshold_min = threshold_min*A,
                           minpeakdistance = minpeakdistance,
                           minpeakheight = max(0.2*A + ylu[1], ymax_min), nups = 1)
        pos_max <- peaks$X
        npeak_MaxPerYear <- length(peaks$gregexpr)/nyear#max peaks

        if (print)
            cat(sprintf('iloop = %d: lambda = %.1f, npeak_MinPerYear = %.2f, npeak_MaxPerYear = %.2f\n',
                iloop, lambda, npeak_MinPerYear, npeak_MaxPerYear))
        # maxpeaksperyear <- 2
        if (npeak_MaxPerYear > max_MaxPeaksperyear | npeak_MinPerYear > max_MinPeaksperyear){
            lambda <- lambda*2
        }else if (npeak_MaxPerYear < 1  | npeak_MinPerYear < 1){
            lambda <- lambda/2
        }else{
            break
        }
        iloop <- iloop + 1
    }

    if (is.null(pos_max) || is.null(pos_min)){
        warning("Can't find a complete growing season before trim!")
        return(NULL)
    }
    # plot(ypred, type = "b"); grid()

    # 1.1 the local minimum value should small than rymin_less * A
    pos_min <- pos_min[(val - ylu[1]) <= rymin_less *A, ]
    pos_min[, type := -1]
    pos_max[, type :=  1]
    pos <- rbind(pos_min, pos_max)
    pos <- pos[order(pos), ] # c("val", "pos", "left", "right", "type")

    # 1.2 remove both points (date or value of min and max too close)
    # I_del <- union(I_date, I_date+1)
    # I_del <- union(I_date + 1, I_val + 1)
    for (i in 1:2){
        if (i == 1){
            # remove dates too close first, then rm continue max or min value
            I_del  <- which(diff(pos$pos) < (nptperyear/12*1)) + 1 # for date
        } else if(i == 2){
            # remove values too close second, then rm continue max or min value again.
            I_del  <- which(abs(diff(pos$val)) < 0.1*A) + 1 #for value
        }
        if (length(I_del) > 0) pos <- pos[-I_del, ]
        # 1.3 remove replicated
        pos$flag <- cumsum(c(1, diff(pos$type) != 0))
        pos      <- ddply(pos, .(flag), rm_duplicate, y = ypred, threshold = threshold_min*A)[, 2:6]
    }
    pos$t    <- t[pos$pos]
    # print(nrow(pos))
    if (nrow(pos) == 0){
        warning("Can't find a complete growing season before!")
        return(NULL)
    }
    ############################################################################
    ## 5. check head and tail break point, and reform breaks
    locals <- pos[, c("pos", "type")]
    ns     <- nrow(locals)
    # check the head and tail minimum values
    minlen <- nptperyear/3 #distance from peak point
    if (last(pos$type) == 1 && (n - nth(pos$pos, -2)) > minlen &&
        abs(last(ypred) - nth(pos$val, -2)) < 0.15*A )
        locals %<>% rbind.data.frame(., data.frame(pos = n, type = -1))
    if (pos$type[1] == 1 && pos$pos[2] > minlen && abs(ypred[1] - pos$val[2]) < 0.15*A)
        locals %<>% rbind.data.frame(data.frame(pos = 1, type = -1), .)

    # a complete growing season, from minimum to minimum
    I      <- which(locals$type == -1)
    locals <- locals[I[1]:I[length(I)], ]

    s  <- locals$pos; ns <- length(s)
    if (ns < 3) {
        warning("Can't find a complete growing season!")
        return(NULL)
    }
    locals %<>% mutate(val = ypred[pos], t = t[pos])

    pos_max <- subset(locals, type == 1)
    pos_min <- subset(locals, type == -1)
    ## 6. divide into multiple growing seasons
    di <- data.table(beg  = s[seq(1, ns-1, 2)],
                     peak = s[seq(2, ns, 2)],
                     end  = s[seq(3, ns, 2)])
    di %<>% fix_di(t = t) #fix whole year data missing

    dt <- map(di, ~t[.x]) %>% as.data.table() %>%
        .[, `:=`( y_beg  = ypred[di$beg],
                  y_peak = ypred[di$peak],
                  y_end  = ypred[di$end],
                  len    = as.integer(difftime(end, beg, units = "days") + 1),
                  year   = year(peak) )]
    fix_dt(dt) # c++ address operation
    # get the growing season year, not only the calendar year
    if (south) dt[, year := year + as.integer(peak >= ymd(sprintf('%d0701', year))) - 1L]

    dt[, `:=`(season = 1:.N, flag   = sprintf("%d_%d", year, 1:.N)), .(year)]

    ## 7. plot
    if (IsPlot){
        # 7.1 PLOT CURVE FITTING TIME-SERIES
        # need to plot outside, because y, w have been changed.
        plotdata(plotdat, nptperyear)
        colors <- c("red", "blue", "green")
        for (i in 1:(length(yfits) - 1)){
            lines(INPUT$t, yfits[[i+1]], col = colors[i], lwd = 2)
        }
        if (!is.null(INPUT$ylu)) abline(h=INPUT$ylu, col="red", lty=2) # show ylims

        # 7.2 plot break points
        I_max <- di$peak
        I_min <- union(di$beg, di$end)
        points(t[I_max], ypred[I_max], pch=20, cex = 1.5, col="red")
        points(t[I_min], ypred[I_min], pch=20, cex = 1.5, col="blue")
    }
    return(list(whit = bind_cols(data_frame(t, y), yfits),
                pos = pos, dt = dt, di = di))
}



curvefits <- function(INPUT, brks, nptperyear = 23,
                      wFUN = wTSM, iters = 2, wmin = 0.1,
                      south = FALSE,
                      extend_month = 2, minT = 0,
                      methods = c('AG', 'zhang', 'beck', 'elmore', 'Gu'),
                      qc,
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
    w  <- brks$whit$w
    di <- brks$di
    if (is.null(di)){
        origin <- t[1]
        di <- data.table(
            beg = brks$dt$beg - origin + 1,
            peak = brks$dt$peak - origin + 1,
            end  = brks$dt$end  - origin + 1)
    }

    # possible snow or cloud, replaced with whittaker smooth.
    I_fix <- which(w == wmin)
    # INPUT$y[I_fix] <- brks$whit %>% {.[[ncol(.)]][I_fix]}
    # w[I_fix]       <- 0.2 # exert the function of whitaker smoother
    # plot(y, type = "b"); grid()
    # lines(brks$whit$iter3, col = "blue")
    # lines(INPUT$y        , col = "red")

    if (debug){
        fits <- stat <- pheno<- NULL
    }else{
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

            #if too much missing values
            if (sum(wi > pmax(wmin, 0.2))/length(wi) < 0.25){
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


curvefits <- function(INPUT, brks, nptperyear = 23,
                      wFUN = wTSM, iters = 2, wmin = 0.1,
                      south = FALSE,
                      extend_month = 2, minT = 0,
                      methods = c('AG', 'zhang', 'beck', 'elmore', 'Gu'),
                      qc,
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
    w  <- brks$whit$w
    di <- brks$di
    if (is.null(di)){
        origin <- t[1]
        di <- data.table(
            beg = brks$dt$beg - origin + 1,
            peak = brks$dt$peak - origin + 1,
            end  = brks$dt$end  - origin + 1)
    }

    # possible snow or cloud, replaced with whittaker smooth.
    I_fix <- which(w == wmin)
    # INPUT$y[I_fix] <- brks$whit %>% {.[[ncol(.)]][I_fix]}
    # w[I_fix]       <- 0.2 # exert the function of whitaker smoother
    # plot(y, type = "b"); grid()
    # lines(brks$whit$iter3, col = "blue")
    # lines(INPUT$y        , col = "red")

    if (debug){
        fits <- stat <- pheno<- NULL
    }else{
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

            #if too much missing values
            if (sum(wi > pmax(wmin, 0.2))/length(wi) < 0.25){
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
