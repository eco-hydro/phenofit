#' @param IsOptim_lambda Whether to optimize Whittaker's parameter lambda by
#' V-curve theory?
#' @param maxExtendMonth Previous and subsequent `maxExtendMonth` data were added
#' for every year curve fitting.
#' @param titlestr string for title
#' @param IsPlot.vc Whether to plot V-curve optimized time-series.
#' @param IsPlot.OnlyBad If true, only plot partial figures whose NSE < 0.3.
#' @param years.run Numeric vector. Which years to run? If not specified, it is
#' all years.
#' @param len_min,len_max the minimum and maximum length (in the unit of days)
#' of growing season
#' 
#' @importFrom lubridate leap_year
#' @rdname season
#' @export
season_mov <- function(INPUT, rFUN, wFUN, iters = 2, wmin = 0.1,
    IsOptim_lambda = FALSE,
    lambda = NULL, nf  = 3, frame = floor(INPUT$nptperyear/5)*2 + 1,
    maxExtendMonth = 12,
    calendarYear = FALSE,
    r_min = 0.05,
    rtrough_max = 0.6,
    ...,
    len_min = 45, len_max = 650,
    .check_season = TRUE,
    years.run = NULL,
    IsPlot = FALSE, IsPlot.vc = FALSE, IsPlot.OnlyBad = FALSE,
    plotdat = INPUT,  titlestr = "")
{
    if (missing(wFUN)) wFUN = get(.options$wFUN_rough)
    if (missing(rFUN)) rFUN = .options$rFUN
    rFUN = check_function(rFUN)
    wFUN = check_function(wFUN)

    nptperyear <- INPUT$nptperyear
    south      <- INPUT$south
    t          <- INPUT$t
    nlen       <- length(t)

    # 1. How many years data
    date_year <- year(t) + ((month(t) >= 7)-1)*south # ecology date
    info  <- table(date_year) # rm years with limited obs
    years <- info[info > nptperyear*0.2] %>% {as.numeric(names(.))}
    nyear <- length(years)

    # 20191128: major update, set `r_min = 0`
    params <- list(
        rFUN = rFUN, wFUN = wFUN, iters = iters, wmin = wmin,
        nf  = nf, frame = frame,
        IsPlot = FALSE, plotdat = plotdat,
        .check_season = .check_season,
        rtrough_max = rtrough_max, r_min = r_min*0, ...)

    has_lambda = !(is.null(lambda) || is.na(lambda))
    brks  <- list()
    vcs   <- vector("list", nyear-2) %>% set_names(years[2:(nyear-1)])

    width_ylu <- nptperyear*0 # already 3y group, moving window for ylu unnecessary

    nextend   <- ceiling(maxExtendMonth/12*nptperyear)
    width_ylu <- nptperyear*2 # This is quite important, to make time-series continuous.

    # modified options：run all years
    # years0 = years[-c(1, nyear)] # original years before `add_HeadTail`
    if (is.null(years.run)) {
        years.run = years
    } else years.run = intersect(years.run, years)

    for (year_i in years.run) {
        i = which(year_i == years)
        if (.options$verbose_season_mov) fprintf("  [season_mov] running %d ... \n", i)

        # I <- which(date_year %in% years[(i-ny_extend):(i+ny_extend)]) # 3y index
        I   <- which(date_year %in% years[i]) # 3y index
        # `nextend` is not enough

        ylu <- get_ylu (INPUT$y, date_year, INPUT$w, width = width_ylu, I, Imedian = TRUE, wmin)
        ylu <- merge_ylu(INPUT$ylu, ylu) # curvefits.R

        # extend curve fitting period, for continuity.
        I <- seq( max(1, first(I) - nextend), min(last(I) + nextend, nlen) )

        input <- lapply(INPUT[c("t", "y", "w")], `[`, I) # y, t, w
        input <- c(input, list(ylu = ylu, nptperyear=nptperyear, south=south))

        if (!has_lambda) {
            vc = guess_lambda(input, wFUN, iters, IsOptim_lambda, IsPlot.vc)
            lambda = vc$lambda; vcs[[i]] <- vc
        }

        params_i = c(list(INPUT = input, lambda = lambda), params)
        brk    <- do.call(season, params_i)

        if (!is.null(brk$dt)){
            brk$dt %<>% subset(year == year_i)
            brk$dt$lambda <- lambda
        }

        if (is.null(brk$dt) || nrow(brk$dt) == 0){
            # if have no brks, try to decrease r_max
            params_i$r_max <- max(params_i$r_max-0.1, 0.05)
            brk <- do.call(season, params_i)
            # we need `rfit` time-series, so can't skip NULL brks.
            if (!is.null(brk$dt)){
                brk$dt %<>% subset(year == year_i)
                brk$dt$lambda <- lambda
            }
        }
        brks[[i]] <- list(whit = brk$whit[date_year[I] == year_i, ], dt   = brk$dt)
    }
    brks  = set_names(brks, years.run) %>% rm_empty() %>% purrr::transpose()
    brks$whit %<>% do.call(rbind, .)

    if (calendarYear) {
        # BUG: need to remove incomplete year
        brks$dt <- season_calendar(years.run, south)
    } else {
        dt  <- do.call(rbind, brks$dt)
        if (is.null(dt)) {
            warning( 'No growing season found!'); return(NULL)
        }
        if (.check_season) {
            brks$dt <- cheak_season_list(brks$dt, rtrough_max, r_min, len_min, len_max)
        }
    }
    brks$GOF <- stat_season(INPUT, brks)

    ## VISUALIZATION
    if (IsPlot) plot_season(INPUT, brks, plotdat, ylu = INPUT$ylu, IsPlot.OnlyBad)
    if (!has_lambda && IsOptim_lambda) brks$optim <- vcs
    return(brks)
}

guess_lambda <- function(input, wFUN = wTSM, iters = 2, IsOptim_lambda = FALSE, IsPlot.vc = FALSE, ...) {
    if (IsOptim_lambda) {
        y <- input$y %>% rm_empty() # should be NA values now
        # update 20181029, add v_curve lambda optimiazaiton in season_mov
        vc <- v_curve(input,
            lg_lambdas = seq(0, 3, by = 0.005), d = 2,
            wFUN = wFUN, iters = iters, IsPlot = IsPlot.vc)
    } else {
        vc <- NULL
        vc$lambda <- init_lambda(input$y) #* 2
    }
    vc
}

#' @rdname check_season
#' @export
check_season_dt <- function(dt, rtrough_max, r_min,
    len_min = 45, len_max = 650)
{
    dt <- dt[len > len_min & len < len_max, ] # mask too long and short gs
    check_season(dt, rtrough_max = rtrough_max, r_min = r_min)
    dt[y_peak != -9999.0 & (len > len_min & len < len_max), ]
}

#' @param dt data.table of growing season dividing info
#' @param lst_dt list of `dt`. Every year is corresponding to a `dt`.
#'
#' @inheritParams season_mov
#' @rdname check_season
#' @export
cheak_season_list <- function(lst_dt, rtrough_max, r_min,
    len_min = 45, len_max = 650)
{
    if (is.data.frame(lst_dt)) lst_dt = list(lst_dt)
    lst_dt %<>% rm_empty()
    res <- list()
    for(i in seq_along(lst_dt)) {
        dt = lst_dt[[i]]
        res[[i]] <- check_season_dt(dt, rtrough_max, r_min, len_min, len_max)
    }
    dt2 = do.call(rbind, res)
    check_season_dt(dt2, rtrough_max, r_min, len_min, len_max)
}

season_calendar <- function(years, south = FALSE){
    date_begin <- ifelse(south, "0701", "0101") %>% paste0(years, .) %>% ymd()
    date_end   <- {if (south) paste0(years+1, "0630") else paste0(years, "1231") } %>% ymd()

    dt <- data.table(
        beg = date_begin,
        end = date_end,
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
#' @param brks A list object returned by `season` or `season_mov`.
#'
#' @keywords internal
#' @rdname season
stat_season <- function(INPUT, brks){
    d_org <- as.data.table(INPUT[c("t", "y0", "w")])
    d_fit <- brks$whit %>% .[,.SD,.SDcols=c(1, ncol(.))] %>% set_colnames(c("t", "ypred"))

    d <- merge(d_org, d_fit, by = "t")

    stat <- with(d, GOF(y0, ypred, w, include.cv = TRUE, include.r = TRUE))# %>% as.list()
    nseason <- ifelse(is.data.frame(brks$dt), nrow(brks$dt), NA)
    stat['nseason'] <- nseason

    # str_title <- sprintf("[%s] IGBP = %s, %s, lat = %.2f", sitename, IGBP_name, stat_str, lat)
    # str_title <- paste(titlestr, stat_txt)
    # NSE <- stat$NSE
    # cv  <- stat$cv
    stat
}
