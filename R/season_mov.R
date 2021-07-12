# ' @param IsPlot.vc Whether to plot V-curve optimized time-series.
# ' @param IsPlot.OnlyBad If true, only plot partial figures whose NSE < 0.3.

#' @name season_mov
#' @title Moving growing season division
#'
#' @inheritParams season
#'
#' @param years.run Numeric vector. Which years to run? If not specified, it is
#' all years.
#' @param options see details
#' @param ... others to [season()]
#'
#' @section options:
#' - `len_min`, `len_max`: minimum and maximum length (in the unit of days)
#' of growing season
#'
#' - `.lambda_vcurve`: Boolean. If the Whittaker's parameter lambda not provided,
#' whether to optimize lambda by V-curve theory? This parameter only works when
#' `lambda` not provided.
#'
#' - `maxExtendMonth`: Previous and subsequent `maxExtendMonth` data were added
#' for every year curve fitting.
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
#' @example R/examples/ex-season.R
#'
#' @importFrom lubridate leap_year
#' @export
season_mov <- function(INPUT,
    # rFUN, wFUN,
    # lambda = NULL, .lambda_vcurve = FALSE,
    options = list(r_min = 0),
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
    # params = as.list(match.call())
    # params = c(as.list(environment()), list(...))
    .options$season %<>% modifyList(options) %>%
        modifyList(list(...))
    .options$season$wFUN %<>% check_function()
    .options$season$rFUN %<>% check_function()
    opt = .options$season

    if (opt$verbose) print(str(opt))

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
        i = which(year_i == years)
        if (.options$verbose_season_mov) fprintf("  [season_mov] running %d ... \n", i)

        ## TODO:
        # `nextend` might be not enough, 3y moving could be an option
        # I <- which(date_year %in% years[(i-ny_extend):(i+ny_extend)]) # 3y index
        I   <- which(date_year %in% years[i])
        ylu <- get_ylu(INPUT$y, date_year, INPUT$w, width = width_ylu, I, Imedian = TRUE,
            opt$wmin)
        ylu %<>% merge_ylu(INPUT$ylu, .) # curvefits.R

        # extend curve fitting period, for continuity.
        I <- seq( max(1, first(I) - nextend), min(last(I) + nextend, nlen) )

        input <- lapply(INPUT[c("t", "y", "w")], `[`, I) # y, t, w
        input <- c(input, list(ylu = ylu, nptperyear=nptperyear, south=south))

        if (!has_lambda) {
            vc = guess_lambda(input) # IsPlot.vc
            lambda = vc$lambda; vcs[[i]] <- vc
        }

        params_i = c(list(INPUT = input, lambda = lambda), dots)
        brk <- do.call(opt_season, params_i)

        if (is_empty(brk$dt)){
            # if have no brks, try to decrease r_max
            params_i$r_max <- max(params_i$r_max-0.1, 0.05)
            brk <- do.call(opt_season, params_i)
        }
        if (!is_empty(brk$dt))
            brk$dt %<>% subset(year == year_i) %>% mutate(lambda = lambda)

        brks[[i]] <- list(
            fit = brk$fit[date_year[I] == year_i, ],
            fit.raw = brk$fit,
            dt = brk$dt)
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
        if (opt$.check_season) brks$dt %<>% cheak_season_list()
    }
    brks$GOF <- stat_season(INPUT, brks$fit)

    if (!has_lambda && opt$.lambda_vcurve) brks$optim <- vcs
    return(brks)
}

# IsPlot.vc = FALSE, IsPlot.OnlyBad = FALSE,
# IsPlot = FALSE, show.legend = TRUE, plotdat = INPUT,  title = ""
# plot_season(INPUT, brks, plotdat, IsPlot.OnlyBad = FALSE, show.legend = show.legend)

guess_lambda <- function(input, IsPlot.vc = FALSE, ...) {
    opt <- .options$season
    if (opt$.lambda_vcurve) {
        y <- input$y %>% rm_empty() # should be NA values now
        # update 20181029, add v_curve lambda optimiazaiton in season_mov
        vc <- v_curve(input,
                      lg_lambdas = seq(0, 5, by = 0.005), d = 2,
                      wFUN = opt$wFUN, iters = opt$iters,
                      IsPlot = IsPlot.vc)
    } else {
        vc <- NULL
        vc$lambda <- init_lambda(input$y) #* 2
    }
    vc
}

#' @rdname check_season
#' @export
check_season_dt <- function(dt) {
    # rtrough_max, r_min, len_min = 45, len_max = 650
    len_min = .options$season$len_min
    len_max = .options$season$len_max

    dt <- dt[len > len_min & len < len_max, ] # mask too long and short gs
    check_season(dt, rtrough_max = .options$season$rtrough_max,
                 r_min = .options$season$r_min)
    dt[y_peak != -9999.0 & (len > len_min & len < len_max), ]
}

#' @param dt data.table of growing season dividing info
#' @param lst_dt list of `dt`. Every year is corresponding to a `dt`.
#'
#' @inheritParams season_mov
#' @rdname check_season
#' @export
cheak_season_list <- function(lst_dt) {
    if (is.data.frame(lst_dt)) lst_dt = list(lst_dt)
    lst_dt %<>% rm_empty()
    res <- list()
    for(i in seq_along(lst_dt)) {
        dt = lst_dt[[i]]
        res[[i]] <- check_season_dt(dt)
    }

    dt2 = do.call(rbind, res)
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

    d_org <- if (!is.data.table(INPUT)) as.data.table(INPUT[c("t", "y0", "w")]) else INPUT

    d <- merge(d_org, d_fit, by = "t")

    stat <- with(d, GOF(y0, ypred, w, include.cv = TRUE, include.r = TRUE))# %>% as.list()
    stat['nseason'] <- nseason

    # str_title <- sprintf("[%s] IGBP = %s, %s, lat = %.2f", sitename, IGBP_name, stat_str, lat)
    # str_title <- paste(titlestr, stat_txt)
    # NSE <- stat$NSE
    # cv  <- stat$cv
    stat
}
