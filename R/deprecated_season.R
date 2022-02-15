#' @title Growing season division
#' @name season
#'
#' @description
#' Divide growing seasons according to rough fitting (`rFUN`) result .
#'
#' For `season`, rough fitting is applied for whole.
#' For `season_mov` rough fitting is applied in every year, during which
#' `maxExtendMonth` is extended.
#'
#' @details
#' Before growing season division, `INPUT` should be added a year in head
#' and tail first by `add_HeadTail`.
#'
#' Finally, use [findpeaks()] to get local maximum and local minimum values.
#' Two local minimum define a growing season.
#' If two local minimum(maximum) are too closed, then only the smaller(biger)
#' is left.
#'
#' @param INPUT A list object with the elements of `t`, `y`, `w`,
#' `Tn` (optional) and `ylu`, returned by [check_input()].
#' @param rFUN Rough curve fitting function, character, can be one of 
#' c("smooth_wSG", "smooth_wWHIT", "smooth_wHANTS"), which are corresponding to 
#' [smooth_wSG()], [smooth_wWHIT()] and [smooth_wHANTS()].
#' @param wFUN weights updating function, can be one of [wTSM()],
#' [wChen()], [wBisquare()] and [wSELF()].
#' @param iters How many times curve fitting is implemented.
#' @param wmin Double, minimum weigth (i.e. weight of snow, ice and cloud).
#' @param lambda The smoothing parameter of [smooth_wWHIT()]. For
#' [season_mov()], if lambda is `NULL`, [init_lambda()]
#' will be used. Generally, it was set as 10000, 15, and 5 for daily, 8-day
#' and 16-day inputs respectively.
#' @param nf The parameter of [smooth_wHANTS()], number of frequencies to be
#' considered above the zero frequency.
#' @param frame The parameter of [smooth_wSG()], moving window size. Suggested by
#' TIMESAT, default `frame = floor(nptperyear/7)*2 + 1`.
#' @param minpeakdistance Numberic, in the unit of points (default as
#' `nptperyear/12`). The minimum distance of two peaks. If the distance of two
#' maximum extreme value less than `minpeakdistance`, only the real maximum
#' value will be left.
#' @param ypeak_min `y_peak >= ypeak_min`
#' @param r_min Threshold is defined as the difference of peak value with
#' trough value. There are two threshold (left and right). The minimum threshold
#' should be greater than r_min.
#' @param r_max Similar as `r_min`, The maximum threshold should
#' be greater than `r_max`.
#' @param rtrough_max `ytrough <= rtrough_max*A`, A is the amplitude of y.
#' @param MaxPeaksPerYear This parameter is used to adjust lambda in iterations.
#' If PeaksPerYear > MaxPeaksPerYear, then lambda = lambda*2.
#' @param MaxTroughsPerYear This parameter is used to adjust lambda in iterations.
#' If TroughsPerYear > MaxTroughsPerYear, then lambda = lambda*2.
#' @param calendarYear If true, only one static calendar growing season will be
#' returned.
#' @param adj.param Adjust rough curve fitting function parameters automatically,
#' if too many or to less peak and trough values.
#' @param rm.closed boolean. Whether check the two closest peaks (or troughs).
#' @param is.continuous boolean. Whether the input is continuous? This parameter
#' is for fluxsite site-year data.
#' @param .check_season not used (only for debug)
#' @param ... ignored.
#' @param verbose whether to print `options_season` into console?
#'
#' @return
#' - `fit`: A data.table of Rough fitting result, with the columns of
#' (`t`, `y`, `witer1`, ..., `witerN`, `ziter1`, ..., `ziterN`).
#'
#' - `dt`: A data.table of growing season division information, with the columns
#' of (`beg`, `peak`, `end`, `y_beg`, `y_peak`, `y_end`, `len`, `year`,
#' `season`, `flag`).
#' @seealso [findpeaks()].
#'
#' @example R/examples/ex-season.R
#' @export
season <- function(
    INPUT, rFUN, wFUN, iters = 2, wmin = 0.1,
    lambda, nf  = 3, frame = floor(INPUT$nptperyear/5)*2 + 1,

    minpeakdistance = NULL,
    ypeak_min = 0.1,
    r_max = 0.2, r_min = 0.05,
    rtrough_max = 0.6,
    MaxPeaksPerYear = 2, MaxTroughsPerYear = 3,
    calendarYear = FALSE,
    adj.param = TRUE,
    rm.closed = TRUE,
    is.continuous = TRUE,
    .check_season = TRUE,
    verbose = FALSE,
    ...)
{
    # message("`season` has been deprecated after v0.3.0!")
    if (missing(wFUN)) wFUN = .options$season$wFUN_rough
    if (missing(rFUN)) rFUN = .options$season$rFUN
    rFUN = check_function(rFUN)
    wFUN = check_function(wFUN)
    
    nptperyear <- INPUT$nptperyear
    south      <- INPUT$south
    t          <- INPUT$t
    y          <- INPUT$y

    nlen       <- length(t)

    # 1. How many years data
    if (is.null(frame)) frame <- floor(nptperyear/5) * 2 + 1
    if (is.null(minpeakdistance)) minpeakdistance <- nptperyear/6

    if (all(is.na(y))) return(NULL)
    c(nyear, years) %<-% guess_nyear(INPUT)

    ylu0  <- INPUT$ylu
    ylu   <- ylu0
    A0    <- diff(ylu0)

    # frame <- floor(nptperyear/7) * 2 + 1 #13, reference by TSM SG filter
    if (missing(lambda)) lambda <- max(nyear*frame, 15)

    ## 1. weighted curve fitting help to divide growing season
    iloop <- 1
    while (iloop <= 3){
        # sgfitw(INPUT$y, INPUT$w, nptperyear, INPUT$ylu, wFUN, iters, frame, d=2)
        # whitsmw(y, w, ylu, wFUN, iters = 1, lambda = 100, ..., d = 2, missval)
        param <- c(INPUT,
            wFUN = wFUN, wmin = wmin, iters = iters,
            lambda = lambda,  # param for whittaker
            nf     = nf,      # param for wHANTS,
            frame  = frame    # param for wSG
        )
        # str(param, 1)
        yfits <- do.call(rFUN, param)
        ypred <- last(yfits$zs) #as.numeric(runmed(ypred, frame))
        alpha <- 0.01

        # default is three year data input, median will be much better
        # This module was primarily designed for `season_mov`. It also works for
        # group large than 3-year. 2-year median will be underestimated.
        if (nyear >= 2.5){ # considering NA values, nyear of 3-year will be smaller.
            ylu_min <- aggregate(ypred, list(year = year(t)), min)$x %>% median()
            ylu_max <- aggregate(ypred, list(year = year(t)), max)$x %>% median()

            # If multiple years are completely missing, ylu_min possiable equals ylu_max
            if (ylu_max - ylu_min > A0*0.2){
                ylu <- c(pmax(ylu_min, INPUT$ylu[1]), #quantile(ypred, alpha/2)
                         pmin(ylu_max, INPUT$ylu[2]))
            }
        }
        INPUT$ylu <- ylu
        A         <- diff(ylu)
        # minPeakHeight <- pmax(ypeak_min, A*0.1 + ylu[1])
        nups <- default_nups(nptperyear)
        info_peak = findpeaks_season(ypred, r_max, r_min,
            minpeakdistance = 0, minpeakheight = ypeak_min,
            nups = nups, nyear = nyear)
        npeak_PerYear   <- info_peak$npeak_PerYear
        ntrough_PerYear <- info_peak$ntrough_PerYear

        if (verbose)
            cat(sprintf('iloop = %d: lambda = %.1f, ntrough_PerYear = %.2f, npeak_PerYear = %.2f\n',
                iloop, lambda, ntrough_PerYear, npeak_PerYear))
        ## This module will automatically update lambda, nf and wHANTS
        if (adj.param) {
            # auto_adjust <- function(nptperyear, npeak_PerYear, ntrough_PerYear, MaxPeaksPerYear, MaxTroughsPerYear, lambda, nf, frame) {
              delta_frame <- ceiling(nptperyear/12)
               # adjust frame in the step of `month`
              if (is.null(npeak_PerYear) || is.null(ntrough_PerYear)) browser()

              if (npeak_PerYear > MaxPeaksPerYear | ntrough_PerYear > MaxTroughsPerYear) {
                  lambda <- lambda*2
                  nf     <- max(1, nf - 1)
                  frame  <- min(frame + delta_frame, nptperyear*2)
              } else if (npeak_PerYear < 0.8  | ntrough_PerYear < 0.8) {
                  lambda <- lambda/2
                  nf     <- min(5, nf + 1)
                  frame  <- max(frame - delta_frame, delta_frame)
              } else {
                  break
              }
            # }
            iloop <- iloop + 1
        } else {
            break
        }
    }
    pos_max = info_peak$pos_max
    pos_min = info_peak$pos_min

    dt <- di <- NULL
    # rough curve fitting time-series
    rfit <- as.data.table(c(list(t = t, y = y), yfits$ws, yfits$zs))
    res  <- list(fit = rfit, dt = dt) # , pos = pos, di = di
    if (is.null(pos_max) || is.null(pos_min)){
        warning("Can't find a complete growing season before trim!")
        return(res)
    }
    # plot(ypred, type = "b"); grid()
    # 1.1 the local minimum value should small than rtrough_max*A

    ## 6. divide into multiple growing seasons
    # update 20180913, commented at 201911221
    # solve the problem of growing season too long, (e.g. US-Cop).
    # I  <- which(dt$y_peak >= ypeak_min)
    # dt <- dt[I, ]
    if (calendarYear) {
        # need to remove incomplete year
        dt <- season_calendar(years, south) # [2:(nyear-1)]
    } else {
        if (rm.closed) {
            pos_min <- pos_min[(val - ylu[1]) <= rtrough_max*A, ] # `y_trough <= rtrough_max*A + ylu[1]`
            pos <- rbind(pos_min, pos_max)[order(pos), ]
            # rm peak value if peak value smaller than the nearest trough values
            I   <- !with(pos, (c(diff(val) > 0, FALSE) & c(diff(type) == -2, FALSE)) |
                    (c(FALSE, diff(val) < 0) & c(FALSE, diff(type) == 2)))
            pos <- pos[I, ]

            pos   <- removeClosedExtreme(pos, ypred, A, y_min = r_min*A)
            pos$t <- t[pos$pos]
            if (nrow(pos) < 2){ # at least two points, begin and end
                warning("Can't find a complete growing season before!"); return(res)
            }
            di  <- check_GS_HeadTail(pos, ypred, minlen = nptperyear/3, A)
        } else {
            # pos <- rbind(pos_min, pos_max)[order(pos), ]
            # di  <- check_GS_HeadTail(pos, ypred, minlen = nptperyear / 3)
            di = pos_max[, .(beg = left, peak = pos, end = right)]
        }
        # fix whole year data missing in FLUXNET data, di: beg, peak and end
        dt <- di2dt(di, t, ypred)
        if (!is.continuous) dt %<>% fixYearBroken(t, ypred)

        dt = dt[len > 45 & len < 650, ] # mask too long and short gs
        if (.check_season) {
            dt %<>% check_season_dt()
            if (!is.continuous) dt %<>% fixYearBroken(t, ypred)
        }
    }

    if (is_empty(dt)) return(NULL)
    # get the growing season year, not only the calendar year
    if (south) dt[, year := year + as.integer(peak >= ymd(sprintf('%d0701', year))) - 1L]
    dt[, `:=`(season = as.numeric(1:.N), flag = sprintf("%d_%d", year, 1:.N)), .(year)]
    res <- list(fit = rfit, dt = dt)
    return(res)
}
# if (IsPlot) plot_season(INPUT, res, plotdat, INPUT)
