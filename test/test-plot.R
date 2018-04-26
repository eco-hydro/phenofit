season <- function(INPUT, lambda, nptperyear = 46, south = FALSE,
                   iters = 2, wFUN = wTSM, IsPlot = TRUE,
                   minpeakdistance = nptperyear/6, ymax_min = 0.5,
                   TRS = 0.05, meth = c('whit', 'sg'), ...,
                   max_MaxPeaksperyear = 2, max_MinPeaksperyear = 3)
{
    t   <- INPUT$t
    y   <- INPUT$y
    ylu <- INPUT$ylu

    if (all(is.na(y))) return(NULL)
    # npt   <- length(y)
    npt   <- sum(INPUT$w > 0)
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
        # }else{
        #     stop('Invalid method input! Should be "sg" or "whit".')
        # }

        ## 4. find local extreme values
        ypred <- yfits[, ncol(yfits), drop = T]
        # ypred <- as.numeric(runmed(ypred, frame))
        alpha <- 0.01
        ylu   <- quantile(ypred, c(alpha/2, 1 - alpha), na.rm = T)
        A     <- diff(ylu)
        threshold <- TRS*A

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
                           threshold_max = 0.2*A,
                           minpeakdistance = minpeakdistance, zero = "-", nups = 0)
        pos_min   <- peaks$X
        pos_min[, 1] %<>% multiply_by(-1)
        npeak_MinPerYear <- length(peaks$gregexpr)/nyear#max peaks
        # local maximum values,
        peaks <- findpeaks(ypred, zero = "+",
                           threshold_max = 0.2*A,
                           threshold_min = 0*A,
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

    # plot curve fitting time-series
    if (IsPlot){
        plotdata(INPUT, nptperyear)
        colors <- c("red", "blue", "green")
        for (i in 1:(ncol(yfits) - 1)){
            lines(INPUT$t, yfits[, i+1, drop = T], col = colors[i], lwd = 2)
        }
    }
    # plot(ypred, type = "b"); grid()

    if (is.null(pos_max) || is.null(pos_min)){
        stop("Can't find a complete growing season before trim!")
    }
    # 1.1 the local minimum value should small than 0.4*A
    pos_min %<>% subset((val - ylu[1]) <= 0.7*A)
    # add column type: max is 1; min is -1.
    # pos column: c("val", "pos", "left", "right", "type")
    pos <- rbind(add_column(pos_max, type = 1), add_column(pos_min, type = -1))
    pos <- pos[order(pos$pos), ]

    # 1.2 remove both points (date or value of min and max too close)
    I_date  <- which(diff(pos$pos) < (nptperyear/12*1)) #15
    I_val   <- which(abs(diff(pos$val)) < 0.1*A) #0.15

    # I_del <- union(I_date, I_date+1)
    I_del <- union(I_date + 1, I_val + 1)
    if (length(I_del) > 0) pos <- pos[-I_del, ]
    pos$flag <- cumsum(c(1, diff(pos$type) != 0))

    # 1.3 remove replicated
    pos   <- ddply(pos, .(flag), rm_duplicate, y = ypred, threshold = threshold)[, 2:6]
    pos$t <- t[pos$pos]

    ############################################################################
    ## 5. check head and tail break point, and reform breaks
    locals <- pos[, c("pos", "type")]
    ns <- nrow(locals)
    # check the head and tail minimum values
    minlen <- nptperyear/3 #distance from peak point
    if (last(pos$type) == 1 && (npt - nth(pos$pos, -2)) > minlen &&
        abs(last(ypred) - nth(pos$val, -2)) < 0.15*A )
        locals %<>% rbind.data.frame(., data.frame(pos = npt, type = -1))
    if (pos$type[1] == 1 && pos$pos[2] > minlen && abs(ypred[1] - pos$val[2]) < 0.15*A)
        locals %<>% rbind.data.frame(data.frame(pos = 1, type = -1), .)

    # a complete growing season, from minimum to minimum
    I      <- which(locals$type == -1)
    locals <- locals[I[1]:I[length(I)], ]

    s  <- locals$pos
    ns <- length(s)
    if (ns < 3) stop("Can't find a complete growing season!")
    locals %<>% mutate(val = ypred[pos], t = t[pos])

    pos_max <- subset(locals, type == 1)
    pos_min <- subset(locals, type == -1)
    ## 6. divide into multiple growing seasons
    di <- data_frame(beg  = s[seq(1, ns-1, 2)]+1,
                     peak = s[seq(2, ns, 2)],
                     end  = s[seq(3, ns, 2)])
    dt <- map_df(di, ~t[.x]) %>%
        mutate(len = difftime(end, beg, units = "days") + 1, year = year(peak)) %>%
        bind_cols(mval = ypred[di$peak], .)
    # get the growing season year, not only the calendar year
    if (south){
        dt %<>% mutate(year = year + as.integer(peak > ymd(sprintf('%d0701', year))) - 1L)
    }
    dt %<>% group_by(year) %>% dplyr::mutate(season = 1:n(), flag = sprintf("%d_%d", year, season))

    ## 7. plot
    if (IsPlot){
        points(pos_max$t, pos_max$val, pch=20, cex = 1.5, col="red")
        points(pos_min$t, pos_min$val, pch=20, cex = 1.5, col="blue")
    }
    # Then curve fitting VI index in every segment, according to local minimum values
    # If begin at maximum value, then add the 1th point as min value. Or just begin
    # from the original the first minimum value.
    return(list(whit = bind_cols(data_frame(t, y), yfits),
                pos = pos, dt = dt, di = di))
}
