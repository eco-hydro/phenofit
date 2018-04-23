# source('R/mainfunc/PhenoBrks.R', encoding = "utf-8")
# ylab = expression(paste('GPP ( gC ', m^-2, d^-1, ')'))
#' @export
plotdata <- function(INPUT, nptperyear, ...){
    t <- INPUT$t
    y <- INPUT$y
    ylu <- INPUT$ylu

    npt <- length(y)
    # show grid lines
    par(mgp = c(1.5, 0.5, 0)) #oma = c(1, 2, 3, 1)
    plot(t, y, xaxt = "s", xaxs = "i", type = "b",
         xlab = 'date', ...)
    # at <- t[seq(1, npt, nptperyear)]
    # fmt <- ifelse(yday(at[1]) == 1, "%Y", "%Y/%m/%d")

    # axis(side=1, at = at, labels = format(at, fmt))
    grid(nx = NA, NULL)
    abline(v=t[seq(1, npt, nptperyear)], col="grey60", lty=3)
    if (!missing(ylu)) abline(h=ylu, col="red", lty=2) # show ylims
}


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
season <- function(INPUT, lambda, nptperyear = 46, south = FALSE,
                   iters = 2, wFUN = wTSM, IsPlot = TRUE,
                   minpeakdistance = nptperyear/6, 
                   ymax_min = 0.5,
                   Aymin_less = 0.6, #ymin < Aymin_less * A
                   TRS = 0.05, meth = c('whit', 'sg'), ...,
                   max_MaxPeaksperyear = 2, max_MinPeaksperyear = 3)
{
    t   <- INPUT$t
    y   <- INPUT$y

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
        ylu[1] <- pmax(ylu[1], INPUT$ylu[1])
        ylu[2] <- pmax(ylu[2], INPUT$ylu[2])
        
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
    # 1.1 the local minimum value should small than Aymin_less * A
    pos_min %<>% subset((val - ylu[1]) <= Aymin_less *A)
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

#' findpeaks
#' @export
findpeaks <- function (x, nups = 1, ndowns = nups, zero = "0", peakpat = NULL,
                       minpeakheight = -Inf, minpeakdistance = 1,
                       threshold_min = 0, threshold_max = 0,
                       npeaks = 0, sortstr = FALSE)
{
    stopifnot(is.vector(x, mode = "numeric") || length(is.na(x)) == 0)
    if (!zero %in% c("0", "+", "-"))
        stop("Argument 'zero' can only be '0', '+', or '-'.")
    xc <- paste(as.character(sign(diff(x))), collapse = "")
    xc <- gsub("1", "+", gsub("-1", "-", xc))
    if (zero != "0")
        xc <- gsub("0", zero, xc)
    if (is.null(peakpat)) {
        peakpat <- sprintf("[+]{%d,}[-]{%d,}", nups, ndowns)
    }
    rc <- gregexpr(peakpat, xc)[[1]]

    if (rc[1] < 0)
        return(NULL)
    x1 <- rc
    x2 <- rc + attr(rc, "match.length")
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
    X <- X[order(X[, 2]), ,drop = F] # update 20180122; order according to index

    if (minpeakdistance < 1)
        warning("Handling 'minpeakdistance < 1' is logically not possible.")
    if (sortstr) { # || minpeakdistance > 1
        sl <- sort.list(X[, 1], na.last = NA, decreasing = TRUE)
        X <- X[sl, , drop = FALSE]
    }
    if (length(X) == 0) return(c())

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

    X <- rm_near(X)
    if (npeaks > 0 && npeaks < nrow(X)) {
        X <- X[1:npeaks, , drop = FALSE]
    }

    X <- setNames(as_tibble(X), c("val", "pos", "left", "right"))
    return(list(gregexpr = rc, X = X))
}
