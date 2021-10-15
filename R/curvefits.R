#' Fine Curve fitting
#'
#' Fine Curve fitting for INPUT time-series.
#'
#' @param INPUT A list object with the elements of 't', 'y', 'w', 'Tn' (optional)
#' and 'ylu', returned by `check_input`.
#' @param brks A list object with the elements of 'fit' and 'dt', returned by
#' `season` or `season_mov`, which contains the growing season
#' dividing information.
#'
#' @param options
#' - `methods` (default `c('AG', 'Beck', 'Elmore', 'Zhang')``):
#' Fine curve fitting methods, can be one or more of
#' `c('AG', 'Beck', 'Elmore', 'Zhang', 'Gu', 'Klos')`. Note that 'Gu' and 'Klos'
#' are very slow.
#'
#' - `wFUN` (default `wTSM`): Character or function, weights updating function
#' of fine fitting function.
#'
#' - `iters` (default 2): max iterations of fine fitting.
#'
#' - `wmin` (default 0.1): min weights in the weights updating procedure.
#'
#' - `use.rough` (default FALSE): Whether to use rough fitting smoothed time-series as input?
#' If `false`, smoothed VI by rough fitting will be used for Phenological metrics
#' extraction; If `true`, original input `y` will be used (rough fitting is used
#' to divide growing seasons and update weights.
#'
#' - `use.y0` (default TRUE): boolean. whether to use original `y0` as the input of `plot_input`,
#' note that not for curve fitting. `y0` is the original value before the process
#' of `check_input`.
#'
#' - `nextend` (default 2): Extend curve fitting window, until `nextend` good or
#' marginal element are found in previous and subsequent growing season.
#'
#' - `maxExtendMonth` (default 1): Search good or marginal good values in previous and
#' subsequent `maxExtendMonth` period.
#'
#' - `minExtendMonth` (default 2): Extending perid defined by `nextend` and `maxExtendMonth`,
#' should be no shorter than `minExtendMonth`.
#' When all points of the input time-series are good value, then the extending
#' period will be too short. In that situation, we can't make sure the connection
#' between different growing seasons is smoothing.
#'
#' - `minPercValid`: (default 0, not use). If the percentage of good
#' and marginal quality points is less than `minPercValid`, curve fiting result
#' is set to `NA`.
#'
#' - `minT`: (not used currently). If `Tn` not provided in `INPUT`, `minT` will not be used.
#' `minT` use night temperature Tn to define backgroud value (days with `Tn < minT`
#' treated as ungrowing season).
#'
#' @param ...  other parameters to [curvefit()]
#'
#' @return List of phenofit fitting object.
#' @seealso [FitDL()]
#'
#' @example inst/examples/ex-curvefits.R
#'
#' @export
curvefits <- function(
    INPUT, brks,
    # methods, wFUN,
    options = list(),
    #   iters = 2, wmin = 0.1,
    #   nextend = 2, maxExtendMonth = 2, minExtendMonth = 1,
    #   minT = 0,
    #   minPercValid = 0,
    #   use.rough = FALSE,
    #   use.y0 = TRUE,
    ...)
{
    if (all(is.na(INPUT$y))) return(NULL)

    .options$fitting %<>% modifyList(options)
    # if (missing(methods)) methods = .options$methods_fine
    # if (missing(wFUN)) wFUN = get(.options$wFUN_fine)
    .options$fitting$wFUN %<>% check_function()
    opt = .options$fitting

    QC_flag    <- INPUT$QC_flag
    nptperyear <- INPUT$nptperyear
    t          <- INPUT$t
    years <- year(t)
    n     <- length(t)

    # doys  <- as.numeric(t) # days from origin
    doys <- as.numeric(difftime(t, date.origin, units = "day")) # + 1

    # Tn for background module
    w  <- w0 <- INPUT$w
    y0 <- if (opt$use.y0) INPUT$y0 else INPUT$y
    # Tn <- INPUT$Tn # if has no Tn, NULL will be return
    # has_Tn <- ifelse(is_empty(Tn), FALSE, TRUE)

    # possible snow or cloud, replaced with Whittaker smoothing.
    I_all <- match(brks$fit$t, t) %>% rm_empty()

    if (opt$use.rough) {
        # if the range of t is smaller than `fit$t`
        INPUT$y[I_all] <- last(brks$fit)
    } else {
        I_fix <- which(w[I_all] == opt$wmin)
        I_y   <- I_all[I_fix]
        INPUT$y[I_y] <- last(brks$fit)[I_fix]
    }
    # use the weights of last iteration of rough fitting
    # w[I_all] <- brks$fit %>% {.[, contain(., "witer"), with = F]} %>% last()
    # w[I_fix] <- wmin + 0.1 # exert the function of whitaker smoother

    # growing season dividing
    di <- data.table( beg  = getDateId_before(brks$dt$beg, t),
                      peak = getDateId_before(brks$dt$peak, t),
                      end  = getDateId_after(brks$dt$end, t)) #%>% na.omit()

    width_ylu = nptperyear*2
    
    y    <- INPUT$y
    fits <- vector(nrow(di), mode = "list")
    for (i in 1:nrow(di)){
        if (opt$verbose) fprintf("  [curvefits] running %d ... \n", i)

        I     <- di$beg[i]:di$end[i]
        I_beg <- di$beg[i]
        I_end <- di$end[i]

        I_extend <- get_extentI(w0, I_beg, I_end, nptperyear)

        ## 2. input data
        ti   <- doys[I_extend]
        yi   <- INPUT$y[I_extend]
        wi   <- w[I_extend]
        # yi_good <- yi[w0[I_extend] > wmin]

        ## update ylu in a three year moving window
        ylu <- get_ylu(y, years, w0, width_ylu, I_beg:I_end, Imedian = TRUE, opt$wmin)
        ylu <- merge_ylu(INPUT$ylu, ylu)
        # yi[yi < ylu[1]] <- ylu[1] # update y value

        # if (has_Tn){
        #     # add background module here, 20180513
        #     Tni        <- Tn[I_extend]
        #     back_value <- backval(yi, ti, wi, Tni, opt$minT, nptperyear)
        #     if (!is.na(back_value)){
        #         I_back     <- yi < back_value
        #         yi[I_back] <- back_value
        #         wi[I_back] <- 0.5
        #     }
        # }
        beginI <- ifelse(i == 1, 1, 2) # make sure no overlap
        tout   <- doys[I] %>% {.[beginI]:last(.)} # make sure return the same length result.

        fFITs  <- curvefit(yi, ti, tout, nptperyear = nptperyear,
                         w = wi, ylu = ylu,
                         iters = opt$iters, methods = opt$methods, wFUN = opt$wFUN, 
                         ...)
        # add original input data here, global calculation can comment this line
        # `y` is original time-series without checked, This for plot
        data <- list(t = doys[I], y = y0[I], QC_flag = QC_flag[I]) %>% as.data.table()
        fFITs$data <- data

        # if too much missing values
        if (sum(wi > pmax(opt$wmin+0.1, 0.2))/length(wi) < opt$minPercValid){
            fFITs$model %<>% map(function(x){
                x$zs %<>% map(~.x*NA) # list()
                return(x)
            })
        }
        fits[[i]] <- fFITs
    }
    # L1:curve fitting method, L2:yearly flag
    fits %<>% set_names(brks$dt$flag)
    fits
    # return(list(tout = t[first(di$beg):last(di$end)],  # dates for OUTPUT curve fitting VI
    #             fits = fits))
}

getDateId <- function(dates, t){
    match(dates, t) #%>% rm_empty()
}
getDateId_before <- function(dates, t) {
    sapply(dates, function(date) {
        ans <- last(which(t <= date))
        if (is.na(ans))
            ans <- first(which(t >= date))
        ans
    })
}
getDateId_after  <- function(dates, t) {
    sapply(dates, function(date) {
        ans <- first(which(t >= date))
        if (is.na(ans))
            ans <- last(which(t <= date))
        ans
    })
}

############################## END OF CURVEFITS ################################
# HIDING FUNCTIONS
# extend curve fitting period
get_extentI <- function(w0, I_beg, I_end, nptperyear) {
    opt <- .options$fitting
    MaxExtendWidth = ceiling(nptperyear / 12 * opt$maxExtendMonth)
    MinExtendWidth = ceiling(nptperyear / 12 * opt$minExtendMonth)

    n <- length(w0)

    I_beg2 <- I_end2 <- NA
    # period <- floor(nptperyear/12*extend_month)
    if (!is.null(opt$nextend)) {
        I_beg_raw <- seq(max(1, I_beg - 1), max(1, I_beg - MaxExtendWidth))
        I_end_raw <- seq(min(n, I_end + 1), min(n, I_end + MaxExtendWidth))

        # at least `nextend` marginal point in previous and following season
        # I_beg2 and I_end2 have been constrained in the range of [1, n]
        I_beg2 <- I_beg_raw[which(w0[I_beg_raw] > opt$wmin)[opt$nextend]]
        I_end2 <- I_end_raw[which(w0[I_end_raw] > opt$wmin)[opt$nextend]]
    }

    # in case of previous and subsequent season good values too closed
    max_Beg <- max(1, I_beg - MinExtendWidth)
    min_End <- min(n, I_end + MinExtendWidth)

    I_beg2 <- ifelse(is.na(I_beg2), max_Beg, min(I_beg2, max_Beg))
    I_end2 <- ifelse(is.na(I_end2), min_End, max(I_end2, min_End))

    if (is_empty(I_beg2) || is_empty(I_end2)) {
        return(I_beg:I_end)
    } else {
        return(I_beg2:I_end2)
    }
}

# merge two limits
merge_ylu <- function(ylu_org, ylu_new){
    if (!is.na(ylu_new[1])) ylu_org[1] <- pmax(ylu_new[1], ylu_org[1])
    if (!is.na(ylu_new[2])) ylu_org[2] <- pmin(ylu_new[2], ylu_org[2])
    ylu_org # return
}

# get ylu in a three year moving window,
# Known issues:
# -------------
# If not complete year, new ylu will be questionable, especially for ylu_max.
# For the continuous of different years, ylu_min is enough.
get_ylu <- function(y, years, w0, width, I, Imedian = TRUE, wmin = 0.2){
    n <- length(w0)
    I_beg  <- I[1]
    I_end  <- last(I)

    # if less than 0.6*12 ≈ 8
    I_allin  <- pmax(1, I_beg - width) : pmin(n, I_end + width)

    w0_win <- w0[I_allin]
    I_allin  <- I_allin[which(w0_win > wmin)]

    if (is_empty(I_allin)) {
        ylu_max <- ylu_min <- NA
    } else {
        y_win  <- y[I_allin]

        ylu_min <- min(y_win)
        ylu_max <- max(y_win)

        if (Imedian){
            year_win  <- years[I_allin]
            # Assume peak in the middle of growing season, a few steps will be
            # ylu_min; while a long way to ylu_max; 20180918
            # length(I) reflects nptperyear
            # d = data.table(y = y_win, year = year_win)
            if (width > length(I)*2/12){
                # ylu_min = d[, min(y), .(year)]$V1 %>% median()
                ylu_min <- median(aggregate.data.frame(y_win, list(year = year_win), min)$x)
            }

            if (width > length(I)*7/12){
                ylu_max <- median(aggregate.data.frame(y_win, list(year = year_win), max)$x)
                # ylu_max <- d[, max(y), .(year)]$V1 %>% median()
            }
        }
    }
    c(ylu_min, ylu_max) # return
}
