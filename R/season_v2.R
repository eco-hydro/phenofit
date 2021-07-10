div_season <- function(ypred, t = seq_along(ypred),
                       r_max = 0.2, r_min = 0.05, minpeakdistance = 30,
                       ypeak_min = 0.1, rtrough_max = 0.6,
                       nptperyear = NULL,
                       calendarYear = FALSE, is.continuous = TRUE, rm.closed = TRUE,
                       .check_season = TRUE,
                       ...) {
    ylu <- range(ypred)
    A <- diff(range(ypred))

    if (is.null(nptperyear) && is.Date(t)) {
        nptperyear <- table(year(t)) %>% median()
    }
    nyear <- length(ypred) / nptperyear

    nups <- ifelse(nptperyear >= 100, 2, 1)
    info_peak <- findpeaks_season(ypred, r_max * A, r_min * A, minpeakdistance,
        ypeak_min, nyear,
        nups = nups
    )

    npeak_PerYear <- info_peak$npeak_PerYear
    ntrough_PerYear <- info_peak$ntrough_PerYear

    pos_max <- info_peak$pos_max
    pos_min <- info_peak$pos_min
    dt0 <- rbind(pos_min, pos_max) # [order(pos), ]

    dt <- di <- NULL
    ## rough curve fitting time-series
    # rfit <- as.data.table(c(list(t = t, y = y), yfits$ws, yfits$zs))
    # res  <- list(fit = rfit, dt = dt) # , pos = pos, di = di
    if (is.null(pos_max) || is.null(pos_min)) {
        warning("Can't find a complete growing season before trim!")
        # return(res)
    }
    
    if (rm.closed) {
        # 1. 移除过大的y_trough,  `y_trough <= rtrough_max*A + ylu[1]`
        pos_min <- pos_min[(val - ylu[1]) <= rtrough_max * A, ]
        cat("1. 移除过高的阈值", "\n")
        print(pos_min[(val - ylu[1]) > rtrough_max * A, ])

        pos <- rbind(pos_min, pos_max)[order(pos), ]
        # rm peak value if peak value smaller than the nearest trough values
        I <- !with(pos, (c(diff(val) > 0, FALSE) & c(diff(type) == -2, FALSE)) |
            (c(FALSE, diff(val) < 0) & c(FALSE, diff(type) == 2)))
        pos <- pos[I, ]

        pos <- adjust_pos(pos, rm.closed, ypred, A, y_min = r_min * A, minpeakdistance)
        pos$t <- t[pos$pos]
        if (nrow(pos) < 2) { # at least two points, begin and end
            warning("Can't find a complete growing season before!")
            return(res)
        }
        di <- check_GS_HeadTail(pos, ypred, minlen = nptperyear / 3, A)
    } else {
        di <- pos_max[, .(beg = left, peak = pos, end = right)]
    }

    # if (calendarYear) {
    #     # need to remove incomplete year
    #     dt <- season_calendar(years, south) # [2:(nyear-1)]
    # } else {
    # fix whole year data missing in FLUXNET data, di: beg, peak and end
    dt <- di2dt(di, t, ypred)
    if (!is.continuous) dt %<>% fixYearBroken(t, ypred)

    dt <- dt[len > 45 & len < 650, ] # mask too long and short gs
    if (.check_season) {
        dt <- check_season_dt(dt, rtrough_max = rtrough_max, r_min = r_min)
        if (!is.continuous) dt %<>% fixYearBroken(t, ypred)
    }
    if (is_empty(dt)) return(NULL)
    dt %<>% rename_season()
    listk(dt0, dt)
    # get the growing season year, not only the calendar year
    # if (south) dt[, year := year + as.integer(peak >= ymd(sprintf('%d0701', year))) - 1L]
    # dt[, `:=`(season = as.numeric(1:.N), flag = sprintf("%d_%d", year, 1:.N)), .(year)]
}

rename_season <- function(d) {
    names(d)[1:6] <- c("time_start", "time_peak", "time_end", "val_start", "val_peak", "val_end")
    d
}
