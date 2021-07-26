#' find_season
#'
#' @param rfit data.frame with the columns of t and `ziter...`, the first column
#' should be `t`, and the last should be `ziter...`.
#'
#' @keywords internal
#' @export
find_season.peaks <- function(
    rfit,
    info_peak,
    options = list(),
    # minpeakdistance = NULL,
    # ypeak_min = 0.1,
    # r_max = 0.2, r_min = 0.05,
    # rtrough_max = 0.6,
    # MaxPeaksPerYear = 2, MaxTroughsPerYear = 3,
    # calendarYear = FALSE,
    # # adj.param = TRUE,
    # rm.closed = TRUE,
    # is.continuous = TRUE,
    # .check_season = TRUE,
    # verbose = FALSE,
    ...)
{
    old <- .options$season
    on.exit(options$season <- old)

    .options$season %<>% modifyList(options) %>%
        modifyList(list(...))
    list2env(.options$season, envir = environment())

    ypred = rfit %>% last() # ypred has no NA
    t = rfit$t
    A = get_A(ypred, na.rm = FALSE)

    if (is.null(minpeakdistance)) minpeakdistance <- info_peak$nptperyear / 6
    # frame <- floor(nptperyear/7) * 2 + 1 #13, reference by TSM SG filter
    pos_max = info_peak$pos_max
    pos_min = info_peak$pos_min

    # rough curve fitting time-series
    if (is.null(pos_max) || is.null(pos_min)){
        warning("Can't find a complete growing season before trim!")
        return(NULL)
    }

    ## 6. divide into multiple growing seasons
    # update 20180913, commented at 201911221
    # solve the problem of growing season too long, (e.g. US-Cop).
    # I  <- which(dt$y_peak >= ypeak_min)
    # dt <- dt[I, ]
    if (calendarYear) {
        # need to remove incomplete year
        dt <- season_calendar(info_peak$year, info_peak$south) # [2:(nyear-1)]
    } else {
        # peaks and toughs put together, and eliminate feaks
        if (rm.closed) {
            # 1.1 the local minimum value should small than rtrough_max*A
            pos_min <- pos_min[(val - min(ypred)) <= rtrough_max * A, ] # `y_trough <= rtrough_max*A + ylu[1]`
            pos <- rbind(pos_min, pos_max)[order(pos), ]
            # rm peak value if peak value smaller than the nearest trough values
            I   <- !with(pos, (c(diff(val) > 0, FALSE) & c(diff(type) == -2, FALSE)) |
                    (c(FALSE, diff(val) < 0) & c(FALSE, diff(type) == 2)))
            pos <- pos[I, ]

            pos   <- removeClosedExtreme(pos, ypred, A, y_min = r_min*A, minpeakdistance)
            pos$t <- t[pos$pos]
            if (nrow(pos) < 2){ # at least two points, begin and end
                warning("Can't find a complete growing season before!"); return(NULL)
            }
            di <- check_GS_HeadTail(pos, ypred, minlen = info_peak$nptperyear / 3, A)
        } else {
            # pos <- rbind(pos_min, pos_max)[order(pos), ]
            # di  <- check_GS_HeadTail(pos, ypred, minlen = nptperyear / 3)
            di = pos_max[, .(beg = left, peak = pos, end = right)]
        }
        # browser()

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
    if (.options$south) dt[, year := year + as.integer(peak >= ymd(sprintf("%d0701", year))) - 1L]
    dt[, `:=`(season = as.numeric(1:.N), flag = sprintf("%d_%d", year, 1:.N)), .(year)]
    dt
}
# if (IsPlot) plot_season(INPUT, res, plotdat, INPUT)

#' @rdname find_season.peaks
#' @export
find_season.default <- function(
    ypred, t = seq_along(ypred),
    nptperyear = NULL, south = NULL,
    options = list(),
    # r_max = 0.2, r_min = 0.05, minpeakdistance = 30,
    # ypeak_min = 0.1, rtrough_max = 0.6,
    # nptperyear = NULL,
    # calendarYear = FALSE, is.continuous = TRUE, rm.closed = TRUE,
    # .check_season = TRUE,
    ...)
{
    old <- .options$season
    on.exit(options$season <- old)

    .options$season %<>% modifyList(options) %>%
        modifyList(list(...))
    opt = .options$season

    if (is.null(nptperyear)) nptperyear = .options$nptperyear
    if (is.null(south)) south = .options$south

    nups <- default_nups(nptperyear)
    info_peak = findpeaks_season(ypred, opt$r_max, opt$r_min,
        minpeakdistance = opt$minpeakdistance, minpeakheight = opt$ypeak_min,
        nups = nups, nyear = NULL)

    years = NULL
    if (is.Date(t)) {
        years = table(year(t)) %>% {.[. > nptperyear*0.2]} %>% names() %>% as.numeric()
    }

    info_peak %<>% c(listk(year = years, nptperyear, south))
    d_fit = data.table(t, ypred)
    find_season.peaks(d_fit, info_peak, .options$season)
}

# ... could also pass to `options`

#' @keywords internal
#' @importFrom utils str
#' @export
opt_season <- function(INPUT,
                       # rFUN, wFUN, lambda,
                       options = NULL, verbose = FALSE, ...) {
    # grasp all parameters
    # params = as.list(environment()) %>% {.[seq(1, length(.)-1)]} # rm the last option
    # .options为全局变量，内部修改存在潜在bug
    .options$season %<>% modifyList(options)

    # OPTIONS_default <- getOption("phenofit.season")
    # options = modifyList(OPTIONS_default, options)
    dots <- list(...)
    ind <- match(names(dots), names(.options$season)) %>% which.notna()
    if (length(ind) > 0) {
        .options$season %<>% modifyList(dots[ind])
        dots <- dots[-ind]
    }
    if (verbose) print(str(.options$season))
    # all of those parameters are put in options

    # separate rough_fitting and find_season
    c(d_fit, info_peak) %<-% rough_fitting(INPUT)
    d_season = find_season.peaks(d_fit, info_peak)
    listk(fit = d_fit, dt = d_season)
    # do.call(season, c(listk(INPUT), .options$season, dots))
}
