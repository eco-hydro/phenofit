#' Growing season dividing in the 3-year length moving window
#'
#' @inheritDotParams season INPUT nptperyear FUN
#' @param d must have columns: 'y', 't', 'w'
#' @param FUN Curve fitting functions, call be one of `sgfitw`, `whitsmw2` and `wHANTS`.
#' 
#' @export
season_3y <- function(INPUT, nptperyear = 23, FUN = whitsmw2,
    south = FALSE,
    IsPlot = T, plotdat = INPUT,
    partial = TRUE, ...){

    nlen  <- length(INPUT$t)
    years <- seq(year(first(INPUT$t)) + 1, year(last(INPUT$t)) - 1)

    nyear <- length(years)

    brks  <- list()
    wFUN  <- wTSM
    ymax_min  <- 0.05

    debug <- FALSE
    # debug <- TRUE
    i = 18; debug <- T
    params <- list(nptperyear = nptperyear,
            FUN = FUN, wFUN = wFUN, iters = 2,
            rymin_less = 0.6, ymax_min = ymax_min,
            threshold_min = 0.0, threshold_max = 0.3,
            IsPlot = debug,
            max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5,
            south = south, plotdat = plotdat)#, ...

    for (i in 1:nyear){
        runningId(i)
        year <- years[i]
        I_beg = nptperyear*(i-1)+1
        I_end = pmin ( nptperyear*(i+2), nlen)

        I      <- I_beg:I_end
        input  <- lapply(INPUT[1:3], `[`, I) %>% c(INPUT[5])
        lambda <- init_lambda(input$y)#*2

        params_i = c(list(INPUT = input, lambda = lambda), params)
        brk    <- do.call(season, params_i)

        brk$dt %<>% subset(year == years[i])
        if (is.null(brk$dt) || nrow(brk$dt) == 0){
            params_i$threshold_max = 0.2
            brk <- do.call(season, params_i)
            brk$dt %<>% subset(year == years[i])
        }
        brks[[i]] <- list(whit = brk$whit[(nptperyear+1):(2*nptperyear), ],
                          dt   = brk$dt)
    }
    brks %<>% purrr::transpose()
    brks$whit %<>% do.call(rbind, .)
    brks$dt   %<>% do.call(rbind, .)
    # brks # quickly return

    dt <- brks$dt
    I_fix <- dt$end[-1] - dt[2:.N]$beg

    if (is.null(brks$dt)){
        warning(sprintf("[site:%s] No data!", sitename))
        return(NULL)
    }
    ## statistics
    I    <- match(brks$whit$t, INPUT$t)
    stat      <- GOF(INPUT$y[I], dplyr::last(brks$whit), INPUT$w[I], include.cv = TRUE)
    stat_str  <- stat[c("R2", "NSE", "cv")] %>% {paste(names(.), round(., 3), sep = "=", collapse = ", ")}
    # str_title <- sprintf("[%s] IGBP = %s, %s, lat = %.2f", sitename, IGBP_name, stat_str, lat)
    str_title <- stat_str

    # lat       <- d$lat[1]
    # Id        <- d$Id[1]
    # IGBP_code <- d$IGBPcode[1]#d$IGBP[1]
    # IGBP_name <- ifelse (!is.null(IGBP_code), IGBPnames[IGBP_code], d$IGBPname)

    nseason           <- nrow(brks$dt)
    stat[['nseason']] <- nseason
    NSE <- stat[['NSE']]
    cv  <- stat[['cv']]
    # if (NSE < 0 | (cv < 0.1 & NSE < 0.1)) {}

    ## VISUALIZATION
    con <- ifelse(partial, IsPlot && (NSE < 0.3), IsPlot)
    if (con){
    # if(IsPlot && (NSE < 0 && cv < 0.2)){
        pdat     <- as.list(d[, .(t, y, w)]) %>% c(INPUT[5])
        plotdata(pdat, nptperyear)

        whit <- brks$whit
        pos  <- brks$dt
        pos_max <- pos$peak

        colors  <- c("red", "blue", "green")
        if (ncol(whit) - 3 < 3) colors <- c("red", "green")

        for (i in 1:(ncol(whit) - 3)){
            lines(whit$t, whit[[i+3]], col = colors[i], lwd = 2)
        }
        # subfunction to plot breaks points
        subplot <- function(t, ...) {
            I <- match(t, whit$t)
            points(t, last(whit)[I], pch = 20, cex = 1.5, ...)
        }
        subplot(pos$peak           , col = "red")
        subplot(c(pos$beg, pos$end), col = "blue")
        title(str_title)
    }
    # stat <- c(site = sitename, IGBPcode = IGBP_code, stat) # quickly return
    return(brks)
}