season <- function(INPUT, lambda, nptperyear = 46, south = FALSE,
                   iters = 2, wFUN = wTSM, IsPlot = TRUE,
                   minpeakdistance = nptperyear/6,
                   ymax_min = 0.5,
                   rymin_less = 0.6, #ymin < rymin_less * A
                   threshold_max = 0.2, threshold_min = 0.05,
                   # TRS = 0.05, meth = c('whit', 'sg'), ...,
                   max_MaxPeaksperyear = 2, max_MinPeaksperyear = 3,
                   wmin = 0.1,
                   plotdat = INPUT,
                   ...)
{
    t   <- INPUT$t
    y   <- INPUT$y

    if (all(is.na(y))) return(NULL)
    # npt   <- length(y)
    npt   <- sum(INPUT$w > wmin)
    nyear <- ceiling(npt/nptperyear)
    # if (nyear <= 3) nyear <- ceiling(nyear)

    frame <- floor(nptperyear/7) * 2 + 1 #13, reference by TSM SG filter
    if (missing(lambda)) lambda <- max(nyear*frame, 15)

    ## 3. weighted whittaker curve fitting
    # wfun <- wTSM#wTSM, bisquare
    iloop <- 1
    while (iloop <= 3){
        # if (meth[1] == 'sg'){
        #     yfits <- sgfitw(INPUT$y, INPUT$w, nptperyear, INPUT$ylu, wFUN, iters, frame, d=2)
        # }else if(meth[1] == 'whit'){
        # whitsmw(y, w, ylu, wFUN, iters = 1, lambda = 100, ..., d = 2, missval)
        yfits <- whitsmw2(INPUT$y, INPUT$w, INPUT$ylu, nptperyear, wFUN, iters, lambda)$data

        ## 4. find local extreme values
        ypred <- yfits[, ncol(yfits), drop = T]
        # ypred <- as.numeric(runmed(ypred, frame))
        alpha <- 0.01

        # default is three year data input, median will be much better
        ylu_min <- aggregate(ypred, list(year = year(t)), min)$x %>% median()
        ylu_max <- aggregate(ypred, list(year = year(t)), max)$x %>% median()
        ylu <- c(pmax(ylu_min, INPUT$ylu[1]), #quantile(ypred, alpha/2)
                 pmin(ylu_max, INPUT$ylu[2]))
        # ylu   <- quantile(ypred, c(alpha/2, 1 - alpha), na.rm = T)
        # ylu[1] <- pmax(ylu[1], INPUT$ylu[1])
        # ylu[2] <- pmax(ylu[2], INPUT$ylu[2])
        A         <- diff(ylu)
        # For avoid fluctuating in peak of growing season or flat peak growing season,
        # like fluxsite: ZM-Mon
        # max_slp <- 2*A/nptperyear
        # pek_slp <- abs(coefficients(lm(ypred[I]~I))[[2]])
        I <- which(ypred > (0.8*A + ylu[1]))
        if (length(I)/length(y) > 0.3){
            ypred[I] <- median(ypred[I])
        }
        # local minimum values
        # threshold for local extreme values
        # peak values is small for minimum values, so can't use threshold here
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
    pos_min %<>% subset((val - ylu[1]) <= rymin_less *A)
    # pos column: c("val", "pos", "left", "right", "type")
    pos <- rbind(add_column(pos_max, type = 1), add_column(pos_min, type = -1))
    pos <- pos[order(pos$pos), ]

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
    ns <- nrow(locals)
    # check the head and tail minimum values
    minlen <- nptperyear/3 #distance from peak point
    n <- length(y)
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
    di <- data_frame(beg  = s[seq(1, ns-1, 2)]+1,
                     peak = s[seq(2, ns, 2)],
                     end  = s[seq(3, ns, 2)])
    di %<>% fix_di(t = t) #fix whole year data missing

    dt <- map_df(di, ~t[.x]) %>%
        mutate(len = difftime(end, beg, units = "days") + 1, year = year(peak)) %>%
        bind_cols(mval = ypred[di$peak], .)
    # get the growing season year, not only the calendar year
    if (south){
        dt %<>% mutate(year = year + as.integer(peak > ymd(sprintf('%d0701', year))) - 1L)
    }
    dt %<>% group_by(year) %>% dplyr::mutate(season = 1:n(), flag = sprintf("%d_%d", year, season))

    ## 7. plot
    #  7.1 PLOT CURVE FITTING TIME-SERIES
    if (IsPlot){
        # need to plot outside, because y, w have been changed.
        plotdata(plotdat, nptperyear)
        colors <- c("red", "blue", "green")
        for (i in 1:(ncol(yfits) - 1)){
            lines(INPUT$t, yfits[, i+1, drop = T], col = colors[i], lwd = 2)
        }
        if (!is.null(INPUT$ylu)) abline(h=INPUT$ylu, col="red", lty=2) # show ylims
    }
    # 7.2 plot break points
    I_max <- di$peak
    I_min <- union(di$beg, di$end)
    if (IsPlot){
        points(t[I_max], ypred[I_max], pch=20, cex = 1.5, col="red")
        points(t[I_min], ypred[I_min], pch=20, cex = 1.5, col="blue")
    }

    return(list(whit = bind_cols(data_frame(t, y), yfits),
                pos = pos, dt = dt, di = di))
}
