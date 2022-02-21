# ' @param IsPlot.vc Whether to plot V-curve optimized time-series.
# ' @param IsPlot.OnlyBad If true, only plot partial figures whose NSE < 0.3.

#' @name season_mov
#' @title Moving growing season division
#' 
#' @inheritParams season
#'
#' @param years.run Numeric vector. Which years to run? If not specified, it is
#' all years.
#' @param options see the following section `options for season` for details.
#' 
#' @param ... others parameter to [set_options()]
#' 
#' @section options for season:
#' #### (a) Parameters for rough fitting
#' - `rFUN`              : character (default `smooth_wWHIT`), the name of rough
#'      curve fitting function, can be one of `c("smooth_wSG", "smooth_wWHIT",
#'      "smooth_wHANTS")`, which are corresponding to [smooth_wSG()],
#'      [smooth_wWHIT()] and [smooth_wHANTS()].
#'
#' - `wFUN`              : character (default `wTSM`), the name of weights
#'      updating functions, can be one of c("wTSM", "wChen", "wBisquare",
#'      "wSELF"). See [wTSM()], [wChen()], [wBisquare()] and [wSELF()] for
#'      details.
#'
#' - `iters`             : integer (default 2), the number of rough fitting
#'   iterations.
#'
#' - `wmin`              : double, the minimum weight of bad points (i.e. snow,
#'   ice and cloud).
#'
#' - `verbose`           : logical (default `FALSE`). If `TRUE`,
#'   `options$season` will be printed on the console.
#'
#' - `lambda`            : double (default NULL), the smoothing parameter of
#'      [smooth_wWHIT()].
#'      + If `lambda = NULL`, V-curve theory will be employed to find the
#'      optimal `lambda`. See [lambda_vcurve()] for details.
#' - `frame`             : integer (default NULL), the parameter of
#'      [smooth_wSG()], moving window size.
#'      + If `frame = NULL`, `frame` will be reset as `floor(nptperyear/5)*2 +
#'      1` (refered by TIMESAT).
#' - `nf`                : integer (default 4), the number of frequencies in
#'   [smooth_wHANTS()].
#'
#' - `maxExtendMonth`: integer (default 12), previous and subsequent
#'   `maxExtendMonth` (in month) data were added to the current year for rough
#'   fitting.
#'
#' - `nextend`           : integer (default NULL), same as `maxExtendMonth`, but
#'    in points.
#'    + If `nextend` provided, `maxExtendMonth` will be ignored.
#'    + If `nextend = NULL`, `nextend` will be reset as
#'      `ceiling(maxExtendMonth/12*nptperyear)`
#'
#' #### (b) Parameters for growing season division
#' - `minpeakdistance`   : double (default NULL), the minimum distance of two
#'   peaks (in points). If the distance of two maximum extreme value less than
#'   `minpeakdistance`, only the maximum one will be kept.
#'   + If `minpeakdistance = NULL`, it will be reset as `nptperyear/6`.
#'
#' - `r_max`             : double (default 0.2; in (0, 1)). `r_max` and `r_min`
#'   are used to eliminate fake peaks and troughs.
#'   + The real peaks should satisfy:
#'      1. \eqn{max(h_{peak, L}, h_{peak, R}) > r_{max} A}
#'      2. \eqn{min(h_{peak, L}, h_{peak, R}) > r_{min} A,} where \eqn{h_{peak,
#'      L}, h_{peak, R}} are height difference from the peak to the left- and
#'      right-hand troughs.
#'   + The troughs should satisfy:
#'      1. \eqn{max(h_{trough, L}, h_{trough, R}) > r_{max} A,} where
#'      \eqn{h_{trough, L}, h_{trough, R}} are height difference from the trough
#'      to the left- and right-hand peaks.
#' - `r_min`             : double (default 0.05; in (0, 1)), see above `r_max`
#'   for details. `r_min` < `r_max`.
#'
#' - `rtrough_max`       : double (default 0.6, in (0, 1)), \eqn{y_{peak} <=
#'   rtrough_max * A + ylu[1]}.
#' - `ypeak_min`         : double 0.1 (in VI unit), \eqn{y_{peak} >= ypeak_min}.
#'
#' - `.check_season`     : logical (default `TRUE`). check the growing season
#'   length according to `len_min` and `len_max`. If `FALSE`, `len_min` and
#'   `len_max` will lose their effect.
#' - `len_min`           : integer (default 45), the minimum length (in days) of
#' growing season
#' - `len_max`           : integer (default 650), the minimum length (in days)
#' of growing season
#'
#' - `adj.param`         : logical. If `TRUE` (default), if there are too many
#'   or too less peaks and troughs, `phenofit` will automatically adjust rough
#'   curve fitting function parameters. See `MaxPeaksPerYear` and
#'   `MaxTroughsPerYear` for details.
#'
#' - `MaxPeaksPerYear` (optional)   : integer (default 2), the max number of
#'   peaks per year. If `PeaksPerYear` > `MaxPeaksPerYear`, then `lambda =
#'   lambda*2`.
#' - `MaxTroughsPerYear` (optional) : integer (default 3), the max number of
#'   troughs per year. If `TroughsPerYear` > `MaxTroughsPerYear`, then `lambda =
#'   lambda*2`.
#'
#' - `calendarYear`      : logical (default `FALSE`). If `TRUE`, the start and
#'   end of a calendar year will be regarded as growing season division (North
#'   Hemisphere is from 01 Jan to 31 Dec; South Hemisphere is from 01 Jul to 30
#'   Jun).
#'
#' - `rm.closed`         : logical (default `TRUE`). If `TRUE`, closed peaks (or troughs)
#'   will be further tidied. Only the maximum
#'
#' - `is.continuous` (not used): logical (default `TRUE`). This parameter is for
#'   `fluxnet2015` fluxsite data, where the input might be not continuous.
#' 
#' @references
#' 1. Kong, D., Zhang, Y., Wang, D., Chen, J., & Gu, X. (2020). Photoperiod
#'    Explains the Asynchronization Between Vegetation Carbon Phenology and
#'    Vegetation Greenness Phenology. Journal of Geophysical Research:
#'    Biogeosciences, 125(8), e2020JG005636.
#'    https://doi.org/10.1029/2020JG005636
#' 2. Kong, D., Zhang, Y., Gu, X., & Wang, D. (2019). A robust method for
#'    reconstructing global MODIS EVI time series on the Google Earth Engine.
#'    ISPRS Journal of Photogrammetry and Remote Sensing, 155, 13-24.
#'
#' @seealso [season()]
#' @example R/examples/ex-season_mov.R
#' 
#' @importFrom lubridate leap_year
#' @export
season_mov <- function(INPUT,
    # rFUN, wFUN,
    # lambda = NULL, 
    options = list(),
    # nf  = 3, frame = floor(INPUT$nptperyear/5)*2 + 1,
    # iters = 2, wmin = 0.1,
    # calendarYear = FALSE,
    # r_min = 0.05,
    # rtrough_max = 0.6,
    # .check_season = TRUE,
    # maxExtendMonth = 12,
    # len_min = 45, len_max = 650,
    # verbose = FALSE,
    ...,
    years.run = NULL)
{
    ## side effect on global parameter
    # opt_old = .options$season
    # on.exit(.options$season <- opt_old)
    set_options(season = options, ...)
    opt = .options$season
    # if (opt$verbose) print(str(opt))

    lambda = opt$lambda
    has_lambda = !(is.null(lambda) || is.na(lambda))
    # r_min = r_min * 0 # 20191128: major update, set `r_min = 0`
    dots <- list(...)

    nptperyear <- INPUT$nptperyear
    south      <- INPUT$south
    t          <- INPUT$t
    nlen       <- length(t)

    # 1. How many years data
    date_year <- year(t) + ((month(t) >= 7)-1)*south # ecology date
    info  <- table(date_year) # rm years with limited obs
    years <- info[info > nptperyear*0.2] %>% {as.numeric(names(.))}
    nyear <- length(years)

    nextend = .options$season$nextend
    if (is.null(nextend)) nextend <- ceiling(.options$season$maxExtendMonth/12 * nptperyear)
    # width_ylu <- nptperyear*0 # already 3y group, moving window for ylu unnecessary
    width_ylu <- nptperyear*2 # This is quite important, to make time-series continuous.

    # years0 = years[-c(1, nyear)] # original years before `add_HeadTail`
    years.run = if (is.null(years.run)) years else intersect(years.run, years)

    brks  <- list()
    vcs   <- vector("list", length(years.run)) %>% set_names(years.run)
    for (year_i in years.run) {
        i = which(year_i == years.run)
        if (.options$season$verbose) fprintf("  [season_mov] running %d ... \n", i)

        # I <- which(date_year %in% years[(i-ny_extend):(i+ny_extend)]) # 3y index
        I   <- which(date_year %in% years.run[i])
        ylu <- get_ylu(INPUT$y, date_year, INPUT$w, width = width_ylu, I, Imedian = TRUE,
            opt$wmin)
        ylu %<>% merge_ylu(INPUT$ylu, .) # curvefits.R

        # extend curve fitting period, for continuity.
        I <- seq( max(1, first(I) - nextend), min(last(I) + nextend, nlen) )

        input <- lapply(INPUT[c("t", "y", "w")], `[`, I) # y, t, w
        input <- c(input, list(ylu = ylu, nptperyear=nptperyear, south=south))

        if (opt$rFUN == "smooth_wWHIT" && !has_lambda) {
            vc <- v_curve(input,
                lg_lambdas = seq(0, 5, by = 0.005), 
                wFUN = opt$wFUN, iters = opt$iters, plot = FALSE
            )
            lambda = vc$lambda; vcs[[i]] <- vc
        }
        params_i = c(list(INPUT = input, lambda = lambda), dots)
        # assign("params_i", params_i, envir = .GlobalEnv)
        brk <- do.call(season_input, params_i)

        if (is_empty(brk$dt)){
            # if have no brks, try to decrease r_max
            params_i$r_max <- max(params_i$r_max-0.1, 0.05)
            brk <- do.call(season_input, params_i)
        }
        if (!is_empty(brk$dt))
            brk$dt %<>% subset(year == year_i) %>% mutate(lambda = lambda)

        # improve for head and tail, v0.3.4
        if (i == 1) {
            rfit = brk$fit[date_year[I] <= year_i, ]
        } else if (i == length(years.run)) {
            rfit = brk$fit[date_year[I] >= year_i, ]
        } else {
            rfit = brk$fit[date_year[I] == year_i, ]
        }
        ans = list(fit = rfit, dt = brk$dt)
        if (.options$debug) ans %<>% c(list(fit.raw = brk$fit), .)
        brks[[i]] <- ans
    }
    brks  = set_names(brks, years.run) %>% rm_empty() %>% purrr::transpose()
    brks$fit %<>% do.call(rbind, .)

    if (opt$calendarYear) {
        # BUG: need to remove incomplete year
        brks$dt <- season_calendar(years.run, south)
    } else {
        dt  <- do.call(rbind, brks$dt)
        if (is_empty(dt)) {
            warning( 'No growing season found!'); return(NULL)
        }
        if (opt$.check_season) brks$dt %<>% check_season_list()
    }
    brks$GOF <- stat_season(INPUT, brks$fit)
    # compatible for previous versions
    if (is.character(opt$rFUN) && opt$rFUN == "smooth_wWHIT" && !has_lambda)
        brks$optim <- vcs
    return(brks)
}

# IsPlot = FALSE, show.legend = TRUE, plotdat = INPUT,  title = ""
# plot_season(INPUT, brks, plotdat, IsPlot.OnlyBad = FALSE, show.legend = show.legend)

#' @keywords internal
#' @rdname rcpp_season_filter
#' @export
check_season_dt <- function(dt) {
    # rtrough_max, r_min, len_min = 45, len_max = 650
    len_min = .options$season$len_min
    len_max = .options$season$len_max

    dt <- dt[len > len_min & len < len_max, ] # mask too long and short gs
    rcpp_season_filter(dt, .options$season$rm.closed,
        rtrough_max = .options$season$rtrough_max,
        r_min = .options$season$r_min)
    dt[y_peak >= .options$season$ypeak_min & (len > len_min & len < len_max), ] %>%
        .[peak <= end & beg < end] # a temporary resolution
}

#' @param dt data.table of growing season division info
#' @param lst_dt list of `dt`. Every year is corresponding to a `dt`.
#'
#' @inheritParams season_mov
#' @rdname rcpp_season_filter
#' @export
check_season_list <- function(lst_dt) {
    if (is.data.frame(lst_dt)) lst_dt = list(lst_dt)
    lst_dt %<>% rm_empty()
    dt2 = map(lst_dt, ~ check_season_dt(.x)) %>%
        do.call(rbind, .)
    check_season_dt(dt2)
}

season_calendar <- function(years, south = FALSE){
    date_begin <- ifelse(south, "0701", "0101") %>% paste0(years, .) %>% ymd()
    date_end   <- {if (south) paste0(years+1, "0630") else paste0(years, "1231") } %>% ymd()

    dt <- data.table(
        beg = date_begin, end = date_end,
        year= years,
        len = as.numeric(difftime(date_end, date_begin, units = "days")) + 1) %>%
        cbind(season = 1, flag = sprintf("%d_1", years))
    dt
}

#' statistics
#'
#' @param d_fit A data.frame with the columns of `t`, `y`, `witer...` and `ziter...`.
#'
#' @keywords internal
#' @rdname season
stat_season <- function(INPUT, d_fit){
    nseason <- ifelse(is.data.frame(d_fit), nrow(d_fit), NA)

    # d_fit <- brks$fit %>% .[,.SD,.SDcols=c(1, ncol(.))] %>% set_colnames(c("t", "ypred"))
    ypred = d_fit %>% dplyr::select(.,starts_with("ziter")) %>% last()
    d_fit = data.table(t = d_fit$t, ypred)

    d_org <- if (!is.data.table(INPUT)) {
        ans = as.data.table(INPUT[c("t", "y0", "w")])
        if(is.null(INPUT$y0)) ans$y0 = INPUT$y
        ans
    } else INPUT
    d <- merge(d_org, d_fit, by = "t")

    stat <- with(d, GOF(y0, ypred, w, include.cv = TRUE, include.r = TRUE))# %>% as.list()
    stat['nseason'] <- nseason

    # str_title <- sprintf("[%s] IGBP = %s, %s, lat = %.2f", sitename, IGBP_name, stat_str, lat)
    # str_title <- paste(titlestr, stat_txt)
    # NSE <- stat$NSE
    # cv  <- stat$cv
    stat
}

# triplicate HANTS test, 2018-09-19
# Not perfect at all for regions with multiple growing season.
# No one method can cope with all the situation.
# {
#     nextend <- length(I)
#     I_beg <- max(1, first(I) - nextend)
#     I_end <- min(last(I) + nextend, nlen)

#     yi <- INPUT$y[I]
#     yhead <- I_beg:(first(I) - 1) %>% { . - .[1] + 1} %>% yi[.]
#     ytail <- (last(I)+1):I_end %>% { . - .[1] + 1} %>% yi[.]
#     yi <- c(yhead, yi, ytail)
#     I <- I_beg:I_end
#     input$y <- yi
# }
