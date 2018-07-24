#' Growing season dividing in the 3-year length moving window
#'
#' @param INPUT A list object with the elements of 't', 'y', 'w', 'Tn' (option)
#' and 'ylu', returned by \code{check_input}.
#' @param nptperyear Integer, points per year.
#' @param south Boolean. In south hemisphere, growing year is 1 July to the
#' following year 31 June; In north hemisphere, growing year is 1 Jan to 31 Dec.
#' @param FUN Coarse curve fitting function, can be one of `sgfitw`, `whitsmw2`
#' and `wHANTS`.
#' @param wmin Double, minimum weigth value (i.e. weight for snow, ice and cloud).
#' @param IsPlot Boolean
#' @param plotdat A list or data.table, with 't', 'y' and 'w'. Only if 
#' IsPlot=true, plotdata will be used to plot. Known that y and w in \code{INPUT} 
#' have been changed, we suggest using the original data.table.
#' @param print Whether to print progress information
#' @param partial If true, only plot partial figures whose NSE < 0.3
#'
#' @return List object, list(whit, dt, stat)
#' @export
season_3y <- function(INPUT, nptperyear = 23, south = FALSE,
    FUN = whitsmw2, wFUN = wTSM, wmin = 0.2,
    IsPlot = T, plotdat = INPUT, print = TRUE,
    partial = TRUE, ...){

    nlen      <- length(INPUT$t)
    date_year <- year(INPUT$t)
    years     <- seq(year(first(INPUT$t)) + 1, year(last(INPUT$t)) - 1)
    nyear     <- length(years)

    ypeak_min <- 0.05
    width_ylu <- nptperyear*1

    debug <- FALSE
    # debug <- TRUE
    # i = 1;
    params <- list(nptperyear = nptperyear, south = south,
            FUN = FUN, wFUN = wFUN, iters = 2,
            rytrough_max = 0.8, ypeak_min = ypeak_min,
            threshold_min = 0.0, threshold_max = 0.3,
            MaxPeaksPerYear = 2.5, MaxTroughsPerYear = 3.5,
            IsPlot = debug, plotdat = plotdat)#, ...

    brks  <- list()
    for (i in 1:nyear){
        if (print) runningId(i, prefix = '\t[season_3y] ')

        year <- years[i]
        I_beg = nptperyear*(i-1)+1
        I_end = pmin ( nptperyear*(i+2), nlen)

        ylu <- get_ylu (INPUT$y, date_year, INPUT$w, width_ylu, I_beg, I_end, Imedian = TRUE, wmin)
        ylu <- merge_ylu(INPUT$ylu, ylu) # curvefits.R

        I      <- I_beg:I_end
        input  <- lapply(INPUT[1:3], `[`, I)
        input$ylu <- ylu

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
    dt <- do.call(rbind, brks$dt)

    if (is.null(dt)){
        warning( 'No growing season found!')
        return(NULL)
    }

    dt %<>% subset(len < 650 & len > 45) # mask too long and short gs
    phenofit:::fix_dt(dt) # c++ address operation, fix growing season overlap
    # after fix_dt, growing season length will become shorter
    dt %<>% subset(len < 650 & len > 45) # mask too long and short gs

    brks$dt <- dt

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

        whit <- brks$whit %>% {.[, contain(., "^t|ziter"), with = F]}
        pos  <- brks$dt
        pos_max <- pos$peak

        iters   <- ncol(whit) - 1 # first row is \code{t}
        colors  <- c("red", "blue", "green")
        if (iters < 3) colors <- c("red", "green")

        for (i in 1:iters){
            lines(whit$t, whit[[i+1]], col = colors[i], lwd = 2)
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
    brks$stat <- stat
    # stat <- c(site = sitename, IGBPcode = IGBP_code, stat) # quickly return
    return(brks)
}
