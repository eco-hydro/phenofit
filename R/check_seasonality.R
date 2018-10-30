#' Check vegetation seasonality
#' @inheritParams season
#' @export
check_seasonality <- function(INPUT, IsPlot = F, plotdat = INPUT, nf = 3, ...){
    # y <- d$EVI
    # t <- as.numeric(d$date - ymd(20000101))
    # w <- d$w
    # nptperyear = 23
    # INPUT <- check_input(t, y, w)
    #
    # pdat <- as.list(d[, .(y = EVI, t = date, w = w)])
    # pdat$ylu <- INPUT$ylu
    fit <- wHANTS(INPUT$y, INPUT$t, INPUT$w, nf = nf, ylu = INPUT$ylu,
                  nptperyear = 23, iters = 3, wFUN = wTSM, wmin = 0.1)
    if (IsPlot){
        plotdata(plotdat, 23)
        colors <- c("red", "blue", "green")
        for (i in 1:(ncol(fit) - 1)){
            lines(plotdat$t, fit[[i+1]], col = colors[i], lwd = 2)
        }
    }
    stat  <- GOF(INPUT$y, dplyr::last(fit), INPUT$w, include.cv = T)
    stat # quickly return
}

#' REMOVE non-seasonality points according to HANTS smooth result
rmNonSeasonality <- function(dt, IsPlot = T, file = 'SI.pdf'){
    if (IsPlot){
        CairoPDF(file, width = 10, height = 2*6)
        par(mfrow = c(6, 1), mar = c(2, 3, 2, 1), mgp = c(1.5, 0.6, 0))
    }

    sites <- unique(dt$site)
    n     <- length(sites)
    res   <- list()
    stats <- list()

    for (i in 1:n){
        runningId(i)
        sitename  <- sites[i]
        d         <- dt[site == sitename]
        IGBP_code <- d$IGBPcode[1]#d$IGBP[1]
        IGBP_name <- IGBPnames[IGBP_code]
        lat       <- d$lat[1]

        # vc <- v_curve(INPUT$y, w = INPUT$w, llas = seq(-2, 2, by = 0.01), d = 2, show = F); lambda <- vc$lambda
        # lambda <- init_lambda(INPUT$y)
        # if ( sum(d$w >= 0.5, na.rm = T) > nrow(d) * 0.3){
            INPUT <- check_input(d$t, d$y, d$w, maxgap = nptperyear / 4, alpha = 0.02)

            INPUT_SI <- INPUT
            INPUT_SI$t %<>% {as.numeric(difftime(., ymd(20000101)))}

            stat       <- check_seasonality(INPUT_SI, IsPlot = FALSE)
            stat_str   <- stat[c("R2", "NSE", "cv")] %>% {paste(names(.), round(., 3), sep = "=", collapse = ", ")}
            stats[[i]] <- c(site = sitename, IGBPcode = IGBP_code, stat)
            str_title  <- sprintf("[%02d] %s, IGBP = %s, %s, lat = %.2f", i, sitename, IGBP_name, stat_str, lat)

            NSE <- stat[['NSE']]
            cv  <- stat[['cv']]

            if (NSE < 0 | (cv < 0.1 & NSE < 0.1)) {
                # plot data for check_SI
                pdat <- as.list(d[, .(y, t, w)]); pdat$ylu <- INPUT$ylu
                temp <- check_seasonality(INPUT_SI, IsPlot, pdat)
                if(IsPlot) title(str_title)
            }
        # } else{
        #     message(str_title)
        # }
    }
    if(IsPlot){
        dev.off()
        file.show(file)
    }

    info <- do.call(rbind, stats) %>% data.table()
    info[, IGBPname := IGBPnames[IGBPcode]]
    info <- info[order(IGBPcode, site)]
    info # quickly return
}

#'
#' Add one year data in the head and tail
#'
#' @param d A data.table, should have \code{t} (compositing date) column
#' (\code{Date} variable).
#' @inheritParams check_input
#' @inheritParams season
#' @param trs If tiny missing (nmissing < trs), than this year is include to 
#' extract phenology.
#' 
#' @return data.table
#' @importFrom lubridate ddays
#' @export
add_HeadTail <- function(d, south = FALSE, nptperyear, trs = nptperyear*0.45){
    if (missing(nptperyear)){
        nptperyear <- ceiling(365/as.numeric(difftime(d$t[2], d$t[1], units = "days")))
    }

    ## can coop with years not continuous now
    ntime    <- nrow(d)

    if (ntime <= 1.2*nptperyear){
        # if only one year data, just triplicate it.
        d_tail <- d_head <- d
    } else {
        step     <- ceiling(365/nptperyear)

        deltaT <- ddays(181)*south
        tt <- d$t - deltaT
        date_year <- year(tt) #+ ((month(t) >= 7)-1)*South

        n_head <- sum(date_year == first(date_year))
        n_tail <- sum(date_year == last(date_year))
        nmissing_head <- nptperyear - n_head
        nmissing_tail <- nptperyear - n_tail

        # if tiny missing, than this year is include to extract phenology
        # if too much, this year will be remove
        nyear_add_head <- ifelse (nmissing_head < trs, 1, 0)
        nyear_add_tail <- ifelse (nmissing_tail < trs, 1, 0)

        na_head <- nptperyear*nyear_add_head + nmissing_head
        na_tail <- nptperyear*nyear_add_tail + nmissing_tail

        I_head <- n_head %>% {seq(.+1, min(.+na_head, ntime))}
        I_tail <- (ntime - n_tail) %>% {seq(max(1, .-na_tail+1), .)}

        # head
        d_head <- d[I_head,]
        d_tail <- d[I_tail,]

    }

    deltaT_head <- first(d$t) - last(d_head$t) - 1
    deltaT_tail <- first(d_tail$t) - last(d$t) - 1

    d_head$t %<>% `+`(deltaT_head)
    d_tail$t %<>% `-`(deltaT_tail)

    # make sure date have no overlap
    # d_head <- d_head[t < date_beg]
    # d_tail <- d_tail[t > date_end]

    res <- rbind(d_head, d, d_tail)
    res # return
}

