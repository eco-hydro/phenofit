#' season
#'
#' First smooth VI timeseries by weighted whittaker, then use findpeak function
#' to get the local maximum and local minimum values. Two local minimum defined
#' a growing season. If two local minimum(maximum) are too closed, then only the
#' smaller(biger) is left.
#'
#' Then according to season pos, based to local maximum position divide yearly
#' growing season. lambda need to set carefully.
#'
#' @param INPUT A list object with the elements of 't', 'y', 'w', 'Tn' (option) 
#' and 'ylu', returned by \code{check_input}.
#' @param nptperyear Integer, points per year.
#' @param south Boolean. In south hemisphere, growing year is 1 July to the 
#' following year 31 June; In north hemisphere, growing year is 1 Jan to 31 Dec.
#' 
#' @param FUN Coarse curve fitting function, can be one of `sgfitw`, `whitsmw2` 
#' and `wHANTS`.
#' @param wFUN weights updating function, can be one of 'wTSM', 'wChen' and 
#' 'wBisquare'.
#' @param iters How many times curve fitting is implemented.
#' @param wmin Double, minimum weigth (i.e. weight of snow, ice and cloud).
#' @param lambda the parameter of \code{whitsmw2}
#' @param nf the parameter of \code{wHANTS}, number of frequencies to be 
#' considered above the zero frequency
#' @param frame the parameter of \code{sgfitw}, moving window size. Suggested by 
#' TIMESAT, default frame = floor(nptperyear/7)*2 + 1.
#' @param minpeakdistance The minimum distance (in indices) peaks have to have
#' to be counted. If the distance of two maximum extreme value less than
#' `minpeakdistance`, only the real maximum value will be left.
#' @param threshold_min Threshold is defined as the difference of peak value with
#' trough value. There are two threshold (left and right). The minimum threshold
#' should be greater than threshold_min.
#' @param threshold_max Similar as `threshold_min`, The maximum threshold should
#' be greater than `threshold_max`.
#' @param ypeak_min ypeak >= ypeak_min 
#' @param rytrough_max ytrough <= rytrough_max*A, A is the amplitude of y.
#' @param MaxPeaksPerYear This parameter is used to adjust lambda in iterations.
#' If PeaksPerYear > MaxPeaksPerYear, then lambda = lambda*2.
#' @param MaxTroughsPerYear This parameter is used to adjust lambda in iterations.
#' If TroughsPerYear > MaxTroughsPerYear, then lambda = lambda*2.
#' @param IsPlot Boolean
#' @param plotdat Is IsPlot = true, plotdata is used to plot original input, 
#' known that y and w in \code{INPUT} have been changed.
#' @param print Whether to print progress information
#' @param ... Other parameters passed to findpeaks
#'
#' @export
#' @return A list object with the elements of 'fit' and 'dt'.
#' list(dt, di)
#' $whit
# # A tibble: 574 x 5
#    t              y     w iter1 iter2
#    <date>     <dbl> <dbl> <dbl> <dbl>
#  1 2002-07-04  9.43 0.149 10.1  10.5
#   $dt
# # A tibble: 11 x 7
# # Groups:   year [11]
#    beg        peak       end        len     year season flag
#    <date>     <date>     <date>     <time> <dbl>  <int> <chr>
#  1 2003-01-09 2003-07-12 2004-02-18 406     2003      1 2003_1
#   $di
# # A tibble: 11 x 3
#      beg  peak   end
#    <dbl> <dbl> <dbl>
#  1  25.0  48.0  76.0
# if more than one continuous maximum(minimum) values, only kept the bigger
# (smaller) one
season <- function(INPUT, nptperyear = 46, south = FALSE,
                   FUN = whitsmw2, wFUN = wTSM, iters = 2, wmin = 0.1,
                   lambda, nf  = 2, frame = floor(nptperyear/7) * 2 + 1,
                   minpeakdistance = nptperyear/6,
                   threshold_max = 0.2, threshold_min = 0.05,
                   ypeak_min   = 0.1, rytrough_max = 0.6,
                   MaxPeaksPerYear = 2, MaxTroughsPerYear = 3,
                   IsPlot  = FALSE, plotdat = INPUT, print = FALSE,
                   ...)
{
    t   <- INPUT$t
    y   <- INPUT$y

    if (all(is.na(y))) return(NULL)
    n     <- length(y)
    npt   <- sum(INPUT$w > wmin)
    nyear <- ceiling(npt/nptperyear)
    ylu0  <- INPUT$ylu

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
        zs    <- yfits$zs
        ypred <- last(zs) #as.numeric(runmed(ypred, frame))
        alpha <- 0.01

        # default is three year data input, median will be much better
        ylu_min <- aggregate(ypred, list(year = year(t)), min)$x %>% median()
        ylu_max <- aggregate(ypred, list(year = year(t)), max)$x %>% median()
        ylu <- c(pmax(ylu_min, INPUT$ylu[1]), #quantile(ypred, alpha/2)
                 pmin(ylu_max, INPUT$ylu[2]))
        A         <- diff(ylu)

        INPUT$ylu <- ylu
        # ylu   <- quantile(ypred, c(alpha/2, 1 - alpha), na.rm = T)
        ## Plateau peak will lead to failed to find local extreme values.
        #  To avoid fluctuating in peak of growing season or flat peak growing
        #  season, like fluxsite: ZM-Mon
        # max_slp <- 2*A/nptperyear
        # pek_slp <- abs(coefficients(lm(ypred[I]~I))[[2]])
        #
        # I <- which(ypred > (0.8*A + ylu[1]))
        # if (length(I)/length(y) > 0.3){
        #     ypred[I] <- median(ypred[I])
        # }
        #
        # local minimum values
        # peak values is small for minimum values, so can't use threshold_min here
        peaks <- findpeaks(-ypred,
                           threshold_max = threshold_max*A,
                           threshold_min = threshold_min*0,
                           minpeakdistance = minpeakdistance, zero = "-", nups = 0)
        pos_min   <- peaks$X
        pos_min[, 1] %<>% multiply_by(-1)
        ntrough_PerYear <- length(peaks$gregexpr)/nyear#max peaks
        # local maximum values,
        peaks <- findpeaks(ypred, zero = "+",
                           threshold_max = threshold_max*A,
                           threshold_min = threshold_min*A,
                           minpeakdistance = minpeakdistance,
                           minpeakheight = max(0.2*A + ylu[1], ypeak_min), nups = 1)
        pos_max <- peaks$X
        npeak_PerYear <- length(peaks$gregexpr)/nyear#max peaks

        if (print)
            cat(sprintf('iloop = %d: lambda = %.1f, ntrough_PerYear = %.2f, npeak_PerYear = %.2f\n',
                iloop, lambda, ntrough_PerYear, npeak_PerYear))
        # maxpeaksperyear <- 2
        if (npeak_PerYear > MaxPeaksPerYear | ntrough_PerYear > MaxTroughsPerYear){
            lambda <- lambda*2
        }else if (npeak_PerYear < 1  | ntrough_PerYear < 1){
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

    # 1.1 the local minimum value should small than rytrough_max * A
    pos_min <- pos_min[(val - ylu[1]) <= rytrough_max *A, ]
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
    # get the growing season year, not only the calendar year
    if (south) dt[, year := year + as.integer(peak >= ymd(sprintf('%d0701', year))) - 1L]

    dt[, `:=`(season = 1:.N, flag   = sprintf("%d_%d", year, 1:.N)), .(year)]

    ## 7. plot
    if (IsPlot){
        # 7.1 PLOT CURVE FITTING TIME-SERIES
        # need to plot outside, because y, w have been changed.
        plotdata(plotdat, nptperyear)
        colors <- c("red", "blue", "green")
        for (i in 1:(length(zs))){
            lines(INPUT$t, zs[[i]], col = colors[i], lwd = 2)
        }
        if (is.null(INPUT$ylu)) abline(h=INPUT$ylu, col="red", lty=2) # show ylims

        # 7.2 plot break points
        I_max <- di$peak
        I_min <- union(di$beg, di$end)
        points(t[I_max], ypred[I_max], pch=20, cex = 1.5, col="red")
        points(t[I_min], ypred[I_min], pch=20, cex = 1.5, col="blue")
    }
    return(list(whit = as.data.table(c(list(t = t, y = y), yfits$ws, yfits$zs)),
                pos = pos, dt = dt, di = di))
}

# rm duplicated max or min values
rm_duplicate <- function(d, y, threshold){
    d <- d[, 1:5]
    if (nrow(d) > 1){
        type <- d$type[1] #1(-1) represent max(min)
        # if range amplitude less than TRS, get median
        if (diff(range(d$val)) < threshold){
            I <- floor(median(d$pos))
            data_frame(val = y[I], pos = I,
                       left = min(d$left),
                       right = max(d$right), type = type)
        }else{
            # else, get the local extreme value
            fun <- ifelse(type == 1, which.max, which.min)
            d[fun(d$val), ]
        }
    }else{
        d
    }
}

# fix across multi-year breaks points, when whole year data are missing
# 
# This function is only for fluxsites data
fix_di <- function(di, t){
    for (i in 1:nrow(di)){
        I    <- di$beg[i]:di$end[i]
        # Try to remove NA values at head and tail. Fialed, Na value may not
        # at head or tail.
        # I_nona <- which(!is.na(y[I])) %>% {I[first(.):last(.)]}
        I_nona <- I#checking year brk is enough
        ti     <- t[I_nona]

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
            di$beg[i] <- first(I_nona)
            di$end[i] <- last(I_nona)
        }
    }
    di#quickly return
}