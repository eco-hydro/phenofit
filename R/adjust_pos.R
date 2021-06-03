# ' @param pos data.table with the columns of `"val", "pos", "left", "right", "type"`.
# ' @param rm.closed Boolean
# ' - `TRUE` : remove both points (date or value of peak and trough too close)
# ' - `FALSE`: `pos_max` used
# ' @examples
# ' adjust_pos(pos, ypred, A, y_min, minpeakdistance, minPeakHeight)
# '
# ' @return A data.table with the columns of:
# ' `"val", "pos", "left", "right", "type"`
adjust_pos <- function(pos, rm.closed, ypred, A, y_min, minpeakdistance)
{
    # 1.2 remove both points (date or value of min and max too close)
    # I_del <- union(I_date, I_date+1)
    # I_del <- union(I_date + 1, I_val + 1)
    # A = diff(range(ypred))
    if (rm.closed) {
        for (i in 1:2) {
            if (i == 1) {
                # remove dates too close first, then rm continue max or min value
                I_del  <- which(diff(pos$pos) < minpeakdistance/2) + 1
            } else if(i == 2) {
                # remove values too close second, then rm continue max or min value again.
                # I_del = NULL
                I_del  <- which(abs(diff(pos$val)) < 0.05 * A) + 1
            }
            if (length(I_del) > 0) pos <- pos[-I_del, ]

            # 1.3 remove replicated
            pos$flag <- cumsum(c(1, diff(pos$type) != 0))
            has_dulicated <- duplicated(pos$flag) %>% any()
            if (has_dulicated) {
                pos %<>% group_by(flag) %>%
                    group_modify(~rm_duplicate(.x, y = ypred, threshold = y_min)) %>%
                    ungroup() %>%
                    select(-flag)
            }
            pos
        }
    } else {
        # 对于CRO, GRA, 变化比较剧烈，相邻极值，不进行调整和融合
        pos[type == 1, ]
    }
    pos
}

# rm duplicated max or min values
rm_duplicate <- function(d, y, threshold){
    # d <- d[, 1:5]
    if (nrow(d) > 1) {
        type <- d$type[1] #1(-1) represent max(min)
        # if range amplitude less than TRS, get median
        if (diff(range(d$val)) < threshold) {
            ## BUG found for CRO, 20191120
            #  修改限定，不移除min
            # if (d$type[1] == -1) return(d)
            I <- floor(median(d$pos))
            tibble(val = y[I], pos = I,
                       left = min(d$left),
                       right = max(d$right), type = type)
        } else {
            # else, get the local extreme value
            fun <- ifelse(type == 1, which.max, which.min)
            d[fun(d$val), ]
        }
    } else {
        d
    }
}


findpeaks_season <- function(ypred, y_max = 0, y_min = 0,
    minpeakdistance = 0, minpeakheight = 0,
    nyear = 1,
    nups = 1, ndowns = nups)
{
    # local minimum values
    # peak values is small for minimum values, so can't use r_min here
    peaks <- findpeaks(-ypred,
        y_max = y_max, y_min = y_min * 0,
        minpeakdistance = minpeakdistance, zero = "-", nups = 0)
    pos_min  <- peaks$X
    if (!is.null(pos_min)) {
        pos_min[, 1] %<>% multiply_by(-1)
        pos_min$type <- -1
    }
    ntrough_PerYear <- length(peaks$gregexpr) / nyear #max peaks

    # minpeakheight = 0.1*A + ylu[1]
    # local maximum values,
    peaks   <- findpeaks(ypred, zero = "+",
        y_max = y_max, y_min = y_min,
        minpeakdistance = minpeakdistance,
        minpeakheight = minpeakheight,
        nups = nups, ndowns = ndowns) #, ypeak_min
    pos_max <- peaks$X
    if (!is.null(pos_max)) pos_max$type <- 1

    npeak_PerYear <- length(peaks$gregexpr) / nyear # max peaks
    listk(pos_min, pos_max, ntrough_PerYear, npeak_PerYear)
}
