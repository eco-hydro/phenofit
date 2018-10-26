#'
#' Add one year data in the head and tail
#'
#' @param d A data.table, should have \code{t} (compositing date) or \code{date}
#' (image date) column which are (\code{Date} variable).
#' @inheritParams check_input
#' @inheritParams season
#' @param trs If tiny missing (nmissing < trs), than this year is include to 
#' extract phenology.
#' 
#' @return data.table
#' @importFrom lubridate ddays
#' @export
add_HeadTail <- function(d, south = FALSE, nptperyear, trs = 0.45){
    # date : image date
    # t    : compositing date
    bandname <- intersect(c("t", "date"), colnames(d))[1]

    if (missing(nptperyear)){
        nptperyear <- ceiling(365/as.numeric(difftime(d[[bandname]][2], d[[bandname]][1], units = "days")))
    }

    ntrs     <- nptperyear*trs
    
    ## can coop with years not continuous now
    ntime    <- nrow(d)

    if (ntime <= 1.2*nptperyear){
        # if only one year data, just triplicate it.
        d_tail <- d_head <- d
    } else {
        step   <- ceiling(365/nptperyear)

        deltaT <- ddays(181)*south
        tt <- d[[bandname]] - deltaT
        date_year <- year(tt) #+ ((month(t) >= 7)-1)*South

        n_head <- sum(date_year == first(date_year))
        n_tail <- sum(date_year == last(date_year))
        nmissing_head <- nptperyear - n_head
        nmissing_tail <- nptperyear - n_tail

        # if tiny missing, than this year is include to extract phenology
        # if too much, this year will be remove
        nyear_add_head <- ifelse (nmissing_head < ntrs, 1, 0)
        nyear_add_tail <- ifelse (nmissing_tail < ntrs, 1, 0)

        na_head <- nptperyear*nyear_add_head + nmissing_head
        na_tail <- nptperyear*nyear_add_tail + nmissing_tail

        I_head <- n_head %>% {seq(.+1, min(.+na_head, ntime))}
        I_tail <- (ntime - n_tail) %>% {seq(max(1, .-na_tail+1), .)}

        # head
        d_head <- d[I_head,]
        d_tail <- d[I_tail,]

    }

    deltaT_head <- first(d[[bandname]]) - last(d_head[[bandname]]) - 1
    deltaT_tail <- first(d_tail[[bandname]]) - last(d[[bandname]]) - 1

    d_head[[bandname]] %<>% `+`(deltaT_head)
    d_tail[[bandname]] %<>% `-`(deltaT_tail)

    # make sure date have no overlap
    # d_head <- d_head[t < date_beg]
    # d_tail <- d_tail[t > date_end]

    res <- rbind(d_head, d, d_tail)
    res # return
}

