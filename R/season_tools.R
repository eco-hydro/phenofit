# ' @param t Date object
get_ylu.default <- function(y, t) {
    ylu_min <- aggregate(y, list(year = year(t)), min)$x %>% median()
    ylu_max <- aggregate(y, list(year = year(t)), max)$x %>% median()
    A <- ylu_max - ylu_min
    listk(ylu_min, ylu_max, A)
}

get_A <- function(x, na.rm = FALSE) {
    range(x, na.rm = na.rm) %>% diff()
}

default_nups <- function(nptperyear) {
    ifelse(nptperyear >= 100, 2, 1)
}

default_minpeakdistance <- function(nptperyear) {
    floor(nptperyear / 6)
}

guess_nyear <- function(INPUT) {
    nlen = length(INPUT$y)
    nptperyear = INPUT$nptperyear

    ## solution 1
    nyear = nlen / nptperyear
    ## solution 2
    date_year <- year(INPUT$t) + ((month(INPUT$t) >= 7) - 1) * INPUT$south
    info  <- table(date_year) # rm years with so limited obs
    years <- info[info > nptperyear*0.2] %>% {as.numeric(names(.))}
    # nyear <- length(years)
    # .[2:(length(.)-1)] # rm head and tail filled years

    ## solution 3
    # npt <- sum(INPUT$w > wmin)
    # if (npt == 0) npt = length(INPUT$y)
    # nyear <- ceiling(npt / nptperyear) # matter most for parameter adjustment
    # print(nyear)
    listk(nyear, year = years)
}


# ' @param di with the columns of `beg`, `peak` and `end`
di2dt <- function(di, t, ypred){
    dt = di %>% mutate(across(.fns = ~ ypred[.x], .names = "y_{.col}")) %>%
        mutate(across(beg:end, .fns = ~ t[.x]))

    if (is.Date(t)) {
        dt %>% mutate(len = as.integer(difftime(end, beg, units = "days") + 1), year = year(peak))
    } else {
        dt %>% mutate(len = end - beg + 1, year = NA_integer_)
    }
}

rename_season <- function(d) {
    names(d)[1:6] <- c("time_start", "time_peak", "time_end", "val_start", "val_peak", "val_end")
    d
}

#' rm too closed peaks or troughs
#'
#' @param pos data.table with the columns of `"val", "pos", "left", "right", "type"`.
#'
#' @examples
#' #removeClosedExtreme(pos, ypred, A, y_min)
#'
#' @keywords internal
#' @return A data.table with the columns of:
#' `"val", "pos", "left", "right", "type"`
#' @export
removeClosedExtreme <- function(pos, ypred, A = NULL, y_min)
{
    # 1. remove both points (date or value of min and max too close)
    # I_del <- union(I_date, I_date+1)
    # I_del <- union(I_date + 1, I_val + 1)
    if (is.null(A)) A = max(ypred) - min(ypred)

    I_del <- which(abs(diff(pos$val)) < 0.05 * A) + 1
    if (length(I_del) > 0) pos <- pos[-I_del, ]

    # 2. remove duplicated peaks or troughs
    pos$flag <- cumsum(c(1, diff(pos$type) != 0))
    has_duplicated <- duplicated(pos$flag) %>% any()
    if (has_duplicated) {
        pos %<>% group_by(flag) %>%
            group_modify(~ merge_duplicate(.x, y = ypred, threshold = y_min)) %>%
            ungroup() %>%
            select(-flag)
    }
    pos
}

# if (rm.closed) {
# for (i in 1:2) {
# if (i == 1) {
# remove dates too close first, then rm continue max or min value
# This step not works for sharply change. Removed in v0.3.4
# I_del  <- which(diff(pos$pos) < minpeakdistance/3) + 1
# } else if(i == 2) {
# remove values too close second, then rm continue max or min value again.
# I_del = NULL
# }
# } else {
#     # 对于CRO, GRA, 变化比较剧烈，相邻极值，不进行调整和融合
#     pos[type == 1, ]
# }
# pos

# 如果在指定阈值之内，则进行合并；
# 否则，仅保留最大的极值
# merge_duplicate duplicated max or min values
merge_duplicate <- function(d, y, threshold, max_gap = 180){
    # d <- d[, 1:5]
    if (nrow(d) > 1) {
        type <- d$type[1] #1(-1) represent max(min)
        # if range amplitude less than TRS, get median

        # 如果距离过远则不进行融合
        diff_val = diff(range(d$val))
        diff_gap = diff(range(d$pos))
        if ( diff_gap < max_gap) {
            if (diff_val < threshold) {
                ## BUG found for CRO, 20191120
                #  修改限定，不移除min
                # if (d$type[1] == -1) return(d)
                I <- floor(median(d$pos))
                data.table(
                    val = y[I], pos = I,
                    left = min(d$left),
                    right = max(d$right), type = type
                )
            } else {
                # else, get the local extreme value
                fun <- ifelse(type == 1, which.max, which.min)
                d[fun(d$val), ]
            }
        } else {
            # bug might exist
            d[1, ]
        }
    } else {
        d
    }
}

# fix across multi-year breaks points, when whole year data are missing
# 用于处理通量站数据中，年不连续的现象
# This function is only for fluxsites data, dt returned
#' @importFrom methods is
fixYearBroken <- function(di, t, ypred) {
    is_date <- is(di$beg[1], "Date")
    origin <- t[1]
    if (is_date) {
        # t is non-continuous sometimes
        I_beg <- match(di$beg, t)
        I_end <- match(di$end, t)
        I_peak <- match(di$peak, t)
    } else {
        I_beg <- di$beg
        I_end <- di$end
        I_peak <- di$peak
    }

    for (i in 1:nrow(di)) {
        I <- I_beg[i]:I_end[i]
        # Try to remove NA values at head and tail. Fialed, Na value may not
        # at head or tail.
        # I_nona <- which(!is.na(y[I])) %>% {I[first(.):last(.)]}
        I_nona <- I # checking year brk is enough
        ti <- t[I_nona]

        # check year brokens
        I_brkyear <- which(diff(ti) >= 365) # broken year
        nbrk <- length(I_brkyear)

        if (nbrk > 0 & nbrk <= 2) {
            if (nbrk == 1) {
                I_1 <- I_nona[1:I_brkyear]
                I_2 <- I_nona[(I_brkyear + 1):length(I_nona)]

                lst <- list(I_1, I_2)
            } else if (nbrk == 2) {
                I_1 <- I_nona[1:I_brkyear[1]]
                I_2 <- I_nona[(I_brkyear[1] + 1):I_brkyear[2]]
                I_3 <- I_nona[(I_brkyear[2] + 1):length(I_nona)]

                lst <- list(I_1, I_2, I_3)
            }
            # select the longest segment, 选择多个片段中最常一段的起始
            I_nona <- lst[[which.max(sapply(lst, length))]]

            # only update dates of beg and end
            I_beg[i] <- first(I_nona)
            I_end[i] <- last(I_nona)
            I_peak[i] <- I_nona[which.max(ypred[I_nona])]
        }
    }

    di <- data.table(beg = I_beg, peak = I_peak, end = I_end)
    di2dt(di, t, ypred)
}

# guess_nyear <- function(INPUT) {
#     npt <- sum(INPUT$w > wmin)
#     ceiling(npt / nptperyear) # matter most for parameter adjustment
# }

#' Check growing season head and tail minimum values
#'
#' @param pos data.frame, with the columns of `pos`, `type`, and `val`.
#' @param minlen `nptperyear/3`, distance from peak point.
#'
#' @keywords internal
#' @return di
check_GS_HeadTail <- function(pos, ypred, minlen, A = NULL, ...) {
    if (is.null(A)) A <- diff(range(ypred))

    nlen <- length(ypred)
    ## 5. check head and tail break point, and reform breaks
    locals <- pos[, c("pos", "type")]
    ns <- nrow(locals)

    # tail: len( pos[end-1], nlen) > minlen, 并且ypred[end]足够小，则认为ypred[end]是trough
    if (last(pos$type) == 1 && (nlen - nth(pos$pos, -2)) > minlen &&
        abs(last(ypred) - nth(pos$val, -2)) < 0.15 * A) {
        locals %<>% rbind.data.frame(., data.frame(pos = nlen, type = -1))
    }
    # head: len( pos[end-1], nlen) > minlen, 并且ypred[1]足够小，则认为ypred[1]是trough
    if (pos$type[1] == 1 && pos$pos[2] > minlen && abs(ypred[1] - pos$val[2]) < 0.15 * A) {
        locals %<>% rbind.data.frame(data.frame(pos = 1, type = -1), .)
    }

    # a complete growing season, from minimum to minimum
    I <- which(locals$type == -1)
    locals <- locals[I[1]:I[length(I)], ]

    s <- locals$pos
    ns <- length(s)
    if (ns < 3) {
        warning("Can't find a complete growing season!")
        return(NULL)
    }
    # locals %<>% mutate(val = ypred[pos], t = t[pos])
    # pos_max <- subset(locals, type == 1)
    # pos_min <- subset(locals, type == -1)
    data.table(
        beg = s[seq(1, ns - 1, 2)], peak = s[seq(2, ns, 2)],
        end = s[seq(3, ns, 2)]
    )
}
