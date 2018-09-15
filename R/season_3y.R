#' Growing season dividing in the 3-year length moving window
#'
#' @inheritParams season
#' @param lambda If lambda is not null, \code{initial_lambda} will be not used.
#' @param titlestr string for title
#' @param IsOnlyPlotbad If true, only plot partial figures whose NSE < 0.3
#' @param ... Other parameters passed to `season`
#'
#' @return List object, list(whit, dt, stat)
#' @export
season_3y <- function(INPUT, south = FALSE,
    FUN = wWHIT, wFUN = wTSM, iters = 2, wmin = 0.1,
    lambda = NULL, nf  = 3, frame = floor(INPUT$nptperyear/5)*2 + 1,
    ...,
    IsPlot = T, plotdat = INPUT, print = TRUE, titlestr = "",
    IsOnlyPlotbad = TRUE)
{
    nptperyear <- INPUT$nptperyear
    t <- INPUT$t
    nlen      <- length(t)
    date_year <- year(t) + ((month(t) >= 7)-1)*south

    # yearly data count info
    info  <- table(date_year) #%>% as.data.table()
    years <- info[info > nptperyear*0.2] %>% {as.numeric(names(.))}
        #.[2:(length(.)-1)] # rm head and tail filled years
    nyear     <- length(years)

    ypeak_min <- 0.05
    width_ylu <- nptperyear*0 # already 3y group, moving window for ylu unnecessary

    debug <- FALSE
    # debug <- TRUE
    # i = 1;
    params <- list(south = south,
            FUN = FUN, wFUN = wFUN, iters = iters, wmin = wmin,
            nf  = nf, frame = frame,
            IsPlot = debug, plotdat = plotdat, ...)#

    has_lambda = !is.null(lambda)
    brks  <- list()

    # If data is not continuous, `season_3y` will be error!
    # Fixed at 20180915
    for (i in 2:(nyear-1)){
        if (print) runningId(i, prefix = '\t[season_3y] ')

        year_i <- years[i]
        I <- which(date_year %in% years[(i-1):(i+1)]) # 3y index

        ylu <- get_ylu (INPUT$y, date_year, INPUT$w, width_ylu, I, Imedian = TRUE, wmin)
        ylu <- merge_ylu(INPUT$ylu, ylu) # curvefits.R

        input  <- lapply(INPUT[1:3], `[`, I)
        input$ylu        <- ylu
        input$nptperyear <- nptperyear

        if (!has_lambda) lambda <- init_lambda(input$y)#*2

        params_i = c(list(INPUT = input, lambda = lambda), params)
        brk    <- do.call(season, params_i)

        brk$dt %<>% subset(year == year_i)
        if (is.null(brk$dt) || nrow(brk$dt) == 0){
            params_i$threshold_max = 0.2
            brk <- do.call(season, params_i)
            brk$dt %<>% subset(year == year_i) # bug found, need to fix for South Hemisphere
        }
        brks[[i-1]] <- list(whit = brk$whit[date_year[I] == year_i, ],
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
    str_title <- paste(titlestr, stat_str)

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
    con <- ifelse(IsOnlyPlotbad, IsPlot && (NSE < 0.3), IsPlot)
    if (con){
    # if(IsPlot && (NSE < 0 && cv < 0.2)){
        pdat     <- as.list(plotdat)[c("t", "y", "w")] %>% c(INPUT[5])
        plotdata(pdat)

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
