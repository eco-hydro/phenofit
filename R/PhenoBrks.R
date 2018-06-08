# source('R/mainfunc/PhenoBrks.R', encoding = "utf-8")
#' @export
plotdata <- function(INPUT, nptperyear, wmin = 0.1, ...){
    t <- INPUT$t
    y <- INPUT$y
    w <- INPUT$w
    
    npt <- length(y)
    # show grid lines
    par(mgp = c(1.5, 0.5, 0)) #oma = c(1, 2, 3, 1)
    # at <- t[seq(1, npt, nptperyear)]
    # fmt <- ifelse(yday(at[1]) == 1, "%Y", "%Y/%m/%d")
    # axis(side=1, at = at, labels = format(at, fmt))
    
    wf <- 4 - findInterval(w, c(-Inf, wmin, 0.5, 1), left.open = T)
    
    colors <- c("grey60", "#00BFC4", "#F8766D", "#C77CFF", "blue", "red", "black")
    pch    <- c(19, 15, 4)
    
    plot(t, y, type = "l", ...)
    Ids <- unique(wf)
    for (i in 1:3){
        I = wf == i
        add <- ifelse(i == 1, F, T)
        points(t[I], y[I], pch = pch[i], col = colors[i], cex = 0.8)
    }
    # ylab = expression(paste('GPP ( gC ', m^-2, d^-1, ')'))
    abline(v = t[seq(1, length(y), nptperyear)], col = "grey60", lty = 3)
    grid(nx = NA, NULL)
    ylu <- INPUT$ylu
    if (!is.null(ylu)) abline(h=ylu, col="red", lty=2) # show ylims
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

#' season
#'
#' First smooth VI timeseries by weighted whittaker, then use findpeak function
#' to get the local maximum and local minimum values. Two local minimum defined
#' a growing season. If two local minimum(maximum) was adjacent, then only the
#' smaller(biger) is left.
#'
#' Then according to season pos, based to local maximum position divide yearly
#' growing season. lambda need to set carefully.
#'
#' divide growing season circle just like TIMESAT season function
#' @param INPUT returned list(t, y, w, ylu) object from `check_input`
#' @param ... Other parameters passed to findpeaks, e.g. 'threshold',
#' 'minpeakheight', 'minpeakdistance'
#' @param TRS multiply A (diff(ylu)) becomes the threshold parameter of
#' `findpeak` function
#' @param rymin_less y_tough < rymin_less * A
#' @param plotdat Is IsPlot = true, plotdata need to be passed. 
#' Since, y, w have been changed.
#' @param FUN Curve fitting functions, call be one of `sgfitw`, `whitsmw2` and `wHANTS`.
#' @export
#' @return list(dt, di)
#' $whit
# # A tibble: 574 x 5
#    t              y     w iter1 iter2
#    <date>     <dbl> <dbl> <dbl> <dbl>
#  1 2002-07-04  9.43 0.149 10.1  10.5
#  2 2002-07-12  9.06 0.137  9.81 10.2
#  3 2002-07-20 10.2  1.00   9.56  9.86
#   $dt
# # A tibble: 11 x 7
# # Groups:   year [11]
#    beg        peak       end        len     year season flag
#    <date>     <date>     <date>     <time> <dbl>  <int> <chr>
#  1 2003-01-09 2003-07-12 2004-02-18 406     2003      1 2003_1
#  2 2004-02-26 2004-07-11 2005-02-10 351     2004      1 2004_1
#  3 2005-02-18 2005-07-12 2006-01-09 326     2005      1 2005_1
#   $di
# # A tibble: 11 x 3
#      beg  peak   end
#    <dbl> <dbl> <dbl>
#  1  25.0  48.0  76.0
#  2  77.0  94.0 121
#  3 122   140   163
# if more than one continuous maximum(minimum) values, only kept the bigger
# (smaller) one
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

#' findpeaks
#'
#' Find peaks (maxima) in a time series. This function is modified from
#' \code{pracma::findpeaks}.
#'
#' @param nups minimum number of increasing steps before a peak is reached
#' @param ndowns minimum number of decreasing steps after the peak
#' @param zero can be `+', `-', or `0'; how to interprete succeeding steps
#' of the same value: increasing, decreasing, or special
#' @param peakpat define a peak as a regular pattern, such as the default
#' pattern ``[+]{1,}[-]{1,}''; if a pattern is provided, the parameters
#' \code{nups} and \code{ndowns} are not taken into account
#' @param minpeakheight The minimum (absolute) height a peak has to have
#' to be recognized as such
#' @param minpeakdistance The minimum distance (in indices) peaks have to have
#' to be counted. If the distance of two maximum extreme value less than
#' `minpeakdistance`, only the real maximum value will be left.
#' @param threshold_min Threshold is defined as the difference of peak value with
#' trough value. There are two threshold (left and right). The minimum threshold
#' should be greater than threshold_min.
#' @param threshold_max Similar as `threshold_min`, The maximum threshold should
#' be greater than `threshold_max`.
#' @param npeaks  the number of peaks to return. If \code{sortstr} = true, the 
#' largest npeaks maximum values will be returned; If \code{sortstr} = false, 
#' just the first npeaks are returned in the order of index.
#' @param sortstr Boolean, Should the peaks be returned sorted in decreasing oreder of 
#' their maximum value?
#'
#' @examples
#' x <- seq(0, 1, len = 1024)
#' pos <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81)
#' hgt <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
#' wdt <- c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005)
#' pSignal <- numeric(length(x))
#' for (i in seq(along=pos)) {
#'     pSignal <- pSignal + hgt[i]/(1 + abs((x - pos[i])/wdt[i]))^4
#' }
#'
#' plot(pSignal, type="l", col="navy"); grid()
#' x <- findpeaks(pSignal, npeaks=3, threshold_min=4, sortstr=TRUE)
#' points(x[, 2], x[, 1], pch=20, col="maroon")
#'
#' @export
findpeaks <- function (x, IsDiff = TRUE, nups = 1, ndowns = nups, zero = "0", peakpat = NULL,
                       minpeakheight = -Inf, minpeakdistance = 1,
                       threshold_min = 0, threshold_max = 0,
                       npeaks = 0, sortstr = FALSE)
{
    stopifnot(is.vector(x, mode = "numeric") ||
                  is.vector(x, mode = "logical") || length(is.na(x)) == 0)

    if (minpeakdistance < 1)
        warning("Handling 'minpeakdistance < 1' is logically not possible.")
    if (!zero %in% c("0", "+", "-"))
        stop("Argument 'zero' can only be '0', '+', or '-'.")

    # extend the use of findpeaks:
    # If want to find extreme values, `IsDiff` should be true;
    # If just want to find the continue negative or positive values, just set
    # `IsDiff` as false.
    if (IsDiff){
        xc <- sign(diff(x))
    }else{
        xc <- x
    }
    xc <- paste(as.character(sign(xc)), collapse = "")
    xc <- gsub("1", "+", gsub("-1", "-", xc))
    if (zero != "0")      xc      <- gsub("0", zero, xc)
    if (is.null(peakpat)) peakpat <- sprintf("[+]{%d,}[-]{%d,}", nups, ndowns)

    rc <- gregexpr(peakpat, xc)[[1]]

    if (rc[1] < 0) return(NULL)
    x1 <- rc
    x2 <- rc + attr(rc, "match.length") - 1
    attributes(x1) <- NULL
    attributes(x2) <- NULL
    n <- length(x1)
    xv <- xp <- numeric(n)
    for (i in 1:n) {
        # for duplicated extreme values, get the median
        vals <- x[x1[i]:x2[i]]
        maxI <- which(vals == max(vals, na.rm = T))
        xp[i] <- floor(median(maxI)) + x1[i] - 1
        xv[i] <- x[xp[i]]
    }
    inds <- which(xv >= minpeakheight &
                      xv - pmin(x[x1], x[x2]) >= threshold_max &
                      xv - pmax(x[x1], x[x2]) >= threshold_min)
    X <- cbind(xv[inds], xp[inds], x1[inds], x2[inds])

    if (length(X) == 0) return(NULL)
    # remove near point where dist < minpeakdistance
    rm_near <- function(x){
        I <- which(diff(x[, 2]) < minpeakdistance)[1]
        if (is.na(I) | length(x) == 0){
            return(x)
        }else{
            I_del <- I + as.integer(x[I, 1] > x[I+1, 1]) #remove the min
            x <- x[-I_del, , drop = F]
            if (nrow(x) <= 1) return(x);
            rm_near(x)
        }
    }

    X <- X[order(X[, 2]), ,drop = F] # update 20180122; order according to index
    X <- rm_near(X) # sort index is necessary before `rm_near`

    if (sortstr) { # || minpeakdistance > 1
        sl <- sort.list(X[, 1], na.last = NA, decreasing = TRUE)
        X <- X[sl, , drop = FALSE]
    }

    if (npeaks > 0 && npeaks < nrow(X)) {
        X <- X[1:npeaks, , drop = FALSE]
    }

    X <- setNames(data.table(X), c("val", "pos", "left", "right"))
    return(list(gregexpr = rc, X = X))
}