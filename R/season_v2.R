div_season <- function(ypred, t, r_max = 0.2, r_min = 0.05, minpeakdistance = 30, 
    ypeak_min = 0.1, rtrough_max = 0.6, 
    calendarYear = FALSE, is.continuous = TRUE, rm.closed = TRUE, 
    .check_season = TRUE, 
    ...) 
{
    ylu = range(ypred)
    A = diff(range(ypred))
    nptperyear = table(year(t)) %>% median()
    nyear = length(ypred)/nptperyear

    nups <- ifelse(nptperyear >= 100, 2, 1)

    info_peak = findpeaks_season(ypred, r_max*A, r_min*A, minpeakdistance,
                                 ypeak_min, nyear, nups = nups)

    npeak_PerYear   <- info_peak$npeak_PerYear
    ntrough_PerYear <- info_peak$ntrough_PerYear

    pos_max = info_peak$pos_max
    pos_min = info_peak$pos_min
    dt0 = rbind(pos_min, pos_max)#[order(pos), ]

    dt   <- di <- NULL
    ## rough curve fitting time-series
    # rfit <- as.data.table(c(list(t = t, y = y), yfits$ws, yfits$zs))
    # res  <- list(fit = rfit, dt = dt) # , pos = pos, di = di
    if (is.null(pos_max) || is.null(pos_min)){
        warning("Can't find a complete growing season before trim!")
        # return(res)
    }
    
    if (rm.closed) {
        # 1. 移除过大的y_trough,  `y_trough <= rtrough_max*A + ylu[1]`
        pos_min <- pos_min[(val - ylu[1]) <= rtrough_max*A, ]
        cat("1. 移除过高的阈值", "\n")
        print(pos_min[(val - ylu[1]) > rtrough_max*A, ])
        
        pos <- rbind(pos_min, pos_max)[order(pos), ]
        # rm peak value if peak value smaller than the nearest trough values
        I   <- !with(pos, (c(diff(val) > 0, FALSE) & c(diff(type) == -2, FALSE)) |
                         (c(FALSE, diff(val) < 0) & c(FALSE, diff(type) == 2)))
        pos <- pos[I, ]

        pos   <- adjust_pos(pos, rm.closed, ypred, A, y_min = r_min*A, minpeakdistance)
        pos$t <- t[pos$pos]
        if (nrow(pos) < 2){ # at least two points, begin and end
            warning("Can't find a complete growing season before!"); return(res)
        }

        ## 5. check head and tail break point, and reform breaks
        locals <- pos[, c("pos", "type")]
        ns     <- nrow(locals)
        
        # check the head and tail minimum values
        minlen <- nptperyear/3 #distance from peak point
        # tail: len( pos[end-1], nlen) > minlen, 并且ypred[end]足够小，则认为ypred[end]是trough
        if (last(pos$type) == 1 && (nlen - nth(pos$pos, -2)) > minlen &&
            abs(last(ypred) - nth(pos$val, -2)) < 0.15*A )
            locals %<>% rbind.data.frame(., data.frame(pos = nlen, type = -1))
        # head: len( pos[end-1], nlen) > minlen, 并且ypred[1]足够小，则认为ypred[1]是trough
        if (pos$type[1] == 1 && pos$pos[2] > minlen && abs(ypred[1] - pos$val[2]) < 0.15*A)
            locals %<>% rbind.data.frame(data.frame(pos = 1, type = -1), .)

        # a complete growing season, from minimum to minimum
        I      <- which(locals$type == -1)
        locals <- locals[I[1]:I[length(I)], ]

        s  <- locals$pos; ns <- length(s)
        if (ns < 3) {
            warning("Can't find a complete growing season!"); return(res)
        }

        locals %<>% mutate(val = ypred[pos], t = t[pos])
        pos_max <- subset(locals, type == 1)
        pos_min <- subset(locals, type == -1)

        di <- data.table(beg  = s[seq(1, ns-1, 2)], peak = s[seq(2, ns, 2)],
                         end  = s[seq(3, ns, 2)])
    } else {
        di = pos_max[, .(beg = left, peak = pos, end = right)]
    }

    # if (calendarYear) {
    #     # need to remove incomplete year
    #     dt <- season_calendar(years, south) # [2:(nyear-1)]
    # } else {
        # only used for fluxsites data, di: beg, peak and end
        dt = if (!is.continuous) {
            # fix whole year data missing
            fixYearBroken(di, t, ypred)
        } else di2dt(di, t, ypred)
        dt <- dt[dt$len > 45 & dt$len < 650, ] # mask too long and short gs
        if (.check_season) {
            dt = check_season_dt(dt, rtrough_max = rtrough_max, r_min = r_min)
            if (!is.continuous) dt %<>% fixYearBroken(t, ypred)
        }
    # }
    if (is.null(dt) || nrow(dt) == 0) return(NULL)
    dt %<>% set_names(c("time_start", "time_peak", "time_end", "val_start", "val_peak", "val_end", "len", "year"))
    listk(dt0, dt)
    # get the growing season year, not only the calendar year
    # if (south) dt[, year := year + as.integer(peak >= ymd(sprintf('%d0701', year))) - 1L]
    # dt[, `:=`(season = as.numeric(1:.N), flag = sprintf("%d_%d", year, 1:.N)), .(year)]
}
