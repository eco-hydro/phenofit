# fix across multi-year breaks points, when whole year data are missing
# 用于处理通量站数据中，年不连续的现象
# This function is only for fluxsites data, dt returned
#' @importFrom methods is
fixYearBroken <- function(di, t, ypred){
    is_date = is(di$beg[1], "Date")
    origin = t[1]
    if (is_date) {
        # t is non-continuous sometimes
        I_beg  = match(di$beg , t)
        I_end  = match(di$end , t)
        I_peak = match(di$peak, t)
    } else {
        I_beg  = di$beg
        I_end  = di$end
        I_peak = di$peak
    }
    
    for (i in 1:nrow(di)){
        I    <- I_beg[i]:I_end[i]
        # Try to remove NA values at head and tail. Fialed, Na value may not
        # at head or tail.
        # I_nona <- which(!is.na(y[I])) %>% {I[first(.):last(.)]}
        I_nona <- I #checking year brk is enough
        ti     <- t[I_nona]

        # check year brokens
        I_brkyear <- which(diff(ti) >= 365) # broken year
        nbrk      <- length(I_brkyear)

        if (nbrk > 0 & nbrk <= 2) {
            if (nbrk == 1) {
                I_1 <- I_nona[1:I_brkyear]
                I_2 <- I_nona[(I_brkyear + 1):length(I_nona)]

                lst <- list(I_1, I_2)
            } else if (nbrk == 2) {
                I_1 <- I_nona[1:I_brkyear[1]]
                I_2 <- I_nona[(I_brkyear[1]+1):I_brkyear[2]]
                I_3 <- I_nona[(I_brkyear[2]+1):length(I_nona)]

                lst <- list(I_1, I_2, I_3)
            }
            #select the longest segment, 选择多个片段中最常一段的起始
            I_nona <- lst[[which.max(sapply(lst, length))]]

            # only update dates of beg and end
            I_beg[i]  <- first(I_nona)
            I_end[i]  <- last(I_nona)
            I_peak[i] <- I_nona[which.max(ypred[I_nona])]
        }
    }

    di = data.table(beg = I_beg, peak = I_peak, end = I_end)
    di2dt(di, t, ypred)
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
