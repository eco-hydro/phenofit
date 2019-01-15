#' Growing season dividing in the 3-year length moving window
#'
#' Before using `season_3y`, INPUT should be added a year in the head and tail
#' first by \code{add_HeadTail}.
#' @inheritParams season
#' @param IsOptim_lambda Whether to optimize Whittaker's parameter lambda by
#' V-curve theory?
#' @param maxExtendMonth Previous and subsequent `maxExtendMonth` data were added
#' for every year curve fitting.
#' subsequent `maxExtendMonth` period.
#' @param titlestr string for title
#' @param IsPlot.vc Whether to plot V-curve optimized time-series.
#' @param IsPlot.OnlyBad If true, only plot partial figures whose NSE < 0.3.
#'
#' @return List object, list(whit, dt)
#'
#' @example inst/examples/ex-check_input.R
#' @example inst/examples/ex-season.R
#' 
#' @rdname season
#' @export
season_3y <- function(INPUT,
    rFUN = wWHIT, wFUN = wTSM, iters = 2, wmin = 0.1,
    IsOptim_lambda = FALSE,
    lambda = NULL, nf  = 3, frame = floor(INPUT$nptperyear/5)*2 + 1,
    maxExtendMonth = 12,
    ...,
    IsPlot = TRUE, IsPlot.vc = FALSE, IsPlot.OnlyBad = FALSE,
    plotdat = INPUT, print = TRUE, titlestr = "")
{
    nptperyear <- INPUT$nptperyear
    south      <- INPUT$south
    t          <- INPUT$t

    nlen      <- length(t)
    date_year <- year(t) + ((month(t) >= 7)-1)*south

    # yearly data count info
    info  <- table(date_year) # rm years with so limited obs
    years <- info[info > nptperyear*0.2] %>% {as.numeric(names(.))}
        #.[2:(length(.)-1)] # rm head and tail filled years
    nyear     <- length(years)

    # I_beg <- first(which(date_year == years[2]))
    # I_end <- first(which(date_year == last(years))) - 1
    # I0    <- I_beg:I_end

    ypeak_min <- 0.05
    width_ylu <- nptperyear*0 # already 3y group, moving window for ylu unnecessary

    debug <- FALSE #FALSE # IsPlot = debug
    # debug <- TRUE
    # i = 1;
    params <- list(
        rFUN = rFUN, wFUN = wFUN, iters = iters, wmin = wmin,
        nf  = nf, frame = frame,
        IsPlot = debug, plotdat = plotdat, ...)

    has_lambda = !is.null(lambda)
    brks  <- list()
    vcs   <- vector("list", nyear-2) %>% set_names(years[2:(nyear-1)])

    nextent   <- ceiling(maxExtendMonth/12*nptperyear)
    width_ylu <- nptperyear*2 # This is quite important, to make time-series continuous.
    # If data is not continuous, `season_3y` will be error!
    # Fixed at 20180915
    for (i in 2:(nyear-1)){
        if (print) runningId(i-1, prefix = '\t[season_3y] ')

        year_i <- years[i]
        # I <- which(date_year %in% years[(i-ny_extend):(i+ny_extend)]) # 3y index
        I   <- which(date_year %in% years[i]) # 3y index
        # `nextent` is not enough

        ylu <- get_ylu (INPUT$y, date_year, INPUT$w, width = width_ylu, I, Imedian = TRUE, wmin)
        ylu <- merge_ylu(INPUT$ylu, ylu) # curvefits.R

        # extend curve fitting period, for continuity.
        I <- seq( max(1, first(I) - nextent), min(last(I) + nextent, nlen) )

        input <- lapply(INPUT[c("t", "y", "w")], `[`, I) # y, t, w
        input <- c(input, list(ylu = ylu, nptperyear=nptperyear, south=south))

        if (!has_lambda) {
            if (IsOptim_lambda){
                y <- input$y %>% rm_empty() # should be NA values now

                # update 20181029, add v_curve lambda optimiazaiton in season_3y
                vc <- v_curve(input, lg_lambdas = seq(-1, 2, by = 0.005), d = 2,
                                  wFUN = wFUN, iters = iters,
                        IsPlot = IsPlot.vc)

                lambda <- vc$lambda
                vcs[[i-1]] <- vc
            } else {
                lambda <- init_lambda(input$y) #*2
            }
        }

        params_i = c(list(INPUT = input, lambda = lambda), params)
        brk    <- do.call(season, params_i)

        if (!is.null(brk$dt)){
            brk$dt %<>% subset(year == year_i)
            brk$dt$lambda <- lambda
        }

        if (is.null(brk$dt) || nrow(brk$dt) == 0){
            # if have no brks, try to decrease r_max
            params_i$r_max <- max(params_i$r_max-0.1, 0.05)
            brk <- do.call(season, params_i)
            # we need `rfit` time-series, so can't skip NULL brks.
            if (!is.null(brk$dt)){
                brk$dt %<>% subset(year == year_i)
                brk$dt$lambda <- lambda
            }
        }

        brks[[i-1]] <- list(whit = brk$whit[date_year[I] == year_i, ],
                          dt   = brk$dt)
    }
    brks <- rm_empty(brks) %>% purrr::transpose()
    brks$whit %<>% do.call(rbind, .)
    dt   <- do.call(rbind, brks$dt)

    if (is.null(dt)){
        warning( 'No growing season found!')
        return(NULL)
    }

    dt <- dt[dt$len > 45 & dt$len < 650, ] # mask too long and short gs
    # phenofit:::fix_dt(dt) # c++ address operation, fix growing season overlap
    fix_dt(dt) # c++ address operation, fix growing season overlap
    # after fix_dt, growing season length will become shorter
    dt <- dt[dt$len > 45 & dt$len < 650, ] # mask too long and short gs
    brks$dt <- dt

    ## VISUALIZATION
    if (IsPlot) plot_season(INPUT, brks, plotdat, ylu = INPUT$ylu, IsPlot.OnlyBad)
    if (IsOptim_lambda) brks$optim <- vcs
    return(brks)
}

# # triplicate HANTS test, 2018-09-19
# # Not perfect at all for regions with multiple growing season.
# # No one method can cope with all the situation.
# {
#     nextent <- length(I)
#     I_beg <- max(1, first(I) - nextent)
#     I_end <- min(last(I) + nextent, nlen)

#     yi <- INPUT$y[I]
#     yhead <- I_beg:(first(I) - 1) %>% { . - .[1] + 1} %>% yi[.]
#     ytail <- (last(I)+1):I_end %>% { . - .[1] + 1} %>% yi[.]
#     yi <- c(yhead, yi, ytail)
#     I <- I_beg:I_end
#     input$y <- yi
# }

#' statistics
#' @rdname season
stat_season <- function(INPUT, brks){
    d_org <- as.data.table(INPUT[c("t", "y", "w")])
    d_fit <- brks$whit %>% .[,.SD,.SDcols=c(1, ncol(.))] %>% set_colnames(c("t", "ypred"))
    d <- merge(d_org, d_fit, by = "t")

    stat <- with(d, GOF(y, ypred, w, include.cv = TRUE))# %>% as.list()
    stat['nseason'] <- nrow(brks$dt)

    # str_title <- sprintf("[%s] IGBP = %s, %s, lat = %.2f", sitename, IGBP_name, stat_str, lat)
    # str_title <- paste(titlestr, stat_txt)
    # NSE <- stat$NSE
    # cv  <- stat$cv
    stat
}

