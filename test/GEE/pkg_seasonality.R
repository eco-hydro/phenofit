#' REMOVE non-seasonality points according to HANTS smooth result
rm_nonSIPoints <- function(dt, IsPlot = T, file = 'SI.pdf'){
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

        if ( sum(d$w >= 0.5, na.rm = T) > nrow(d) * 0.3){
            INPUT <- check_input(d$t, d$y, d$w, trim = T, maxgap = nptperyear / 4, alpha = 0.02)

            INPUT_SI <- INPUT
            INPUT_SI$t %<>% {as.numeric(difftime(., ymd(20000101)))}

            stat       <- check_SI(INPUT_SI, FALSE)
            stat_str   <- stat[c("R2", "NSE", "cv")] %>% {paste(names(.), round(., 3), sep = "=", collapse = ", ")}
            stats[[i]] <- c(site = sitename, IGBPcode = IGBP_code, stat)
            # vc <- v_curve(INPUT$y, w = INPUT$w, llas = seq(-2, 2, by = 0.01), d = 2, show = F); lambda <- vc$lambda
            # lambda <- init_lambda(INPUT$y)
            str_title <- sprintf("[%02d] %s, IGBP = %s, %s, lat = %.2f", i, sitename, IGBP_name, stat_str, lat)

            NSE <- stat[['NSE']]
            cv  <- stat[['cv']]

            if (NSE < 0 | (cv < 0.1 & NSE < 0.1)) {
            # if (stat[['NSE']] < 0.2 & stat[['NSE']] > 0) {
            # if (stat[['NSE']] < 0) {
                # plot data for check_SI
                pdat <- as.list(d[, .(y, t, w)]); pdat$ylu <- INPUT$ylu
                stat <- check_SI(INPUT_SI, IsPlot, pdat)
                if(IsPlot) title(str_title)
            }
        } else{
            message(str_title)
        }
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


#' This function is designed for MODIS DATASET
whit_brks <- function(d, nptperyear = 23, IsPlot = T){
    # I_beg <- floor(yday(ymd(20000218))/16) # 0218 is the 4th 16-day
    ntime  <- nrow(d)
    d_head <- d[21:46, ]
    d_head$t %<>% `-`(., dyears(2) - ddays(1))

    d_tail <- d[(ntime - 2*nptperyear + 1):(ntime - nptperyear)]
    d_tail$t %<>% `+`(dyears(2))

    dnew  <- rbind(d_head, d, d_tail)
    INPUT <- check_input(dnew$t, dnew$y, dnew$w, trim = T, maxgap = nptperyear / 4, alpha = 0.02)

    years <- 2000:2016
    nyear <- length(years)

    brks <- list()
    wFUN <- wTSM
    ymax_min <- 0.05
    lat <- d$lat[1]

    for (i in 1:nyear){
        i = 15
        runningId(i)
        year <- years[i]
        I_beg = nptperyear*(i-1)+1
        I_end = nptperyear*(i+2)
        I <- I_beg:I_end
        input  <- lapply(INPUT[1:3], `[`, I) %>% c(INPUT[5])
        lambda <- init_lambda(input$y)
        params <- list(input, lambda = lambda, nptperyear, iters = 3, wFUN = wFUN,
                       rymin_less = 0.6, ymax_min = ymax_min,
                       threshold_min = 0.05, IsPlot = T,
                       max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5,
                       south = lat < 0, plotdat = d)
        brk    <- do.call(season, params)
        brk$dt %<>% subset(year == years[i])
        if (nrow(brk$dt) == 0){
            params$threshold_min = 0
            brk <- do.call(season, params)
        }
        brks[[i]] <- list(whit = brk$whit[(nptperyear+1):(2*nptperyear), ],
                          dt   = brk$dt)
    }
    brks %<>% purrr::transpose()
    brks$whit %<>% do.call(rbind, .)
    brks$dt   %<>% do.call(rbind, .)
    # brks # quickly return

    ## visualization
    if(IsPlot){
        pdat     <- as.list(d[, .(t, y, w)]) %>% c(INPUT[5])
        plotdata(pdat, nptperyear)

        whit <- brks$whit
        pos  <- brks$dt
        pos_max = pos$peak
        colors <- c("red", "blue", "green")
        for (i in 1:3){
            lines(whit$t, whit[[i+3]], col = colors[i], lwd = 2)
        }
        # subfunction to plot breaks points
        subplot <- function(t, ...) {
            I <- match(t, whit$t)
            points(t, whit$iter3[I], pch = 20, cex = 1.5, ...)
        }
        subplot(pos$peak           , col = "red")
        subplot(c(pos$beg, pos$end), col = "blue")
    }
}
