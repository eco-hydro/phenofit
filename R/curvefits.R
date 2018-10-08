#' Fine Curve fitting
#'
#' Fine Curve fitting for INPUT time-series.
#'
#' @inheritParams season
#' @param INPUT A list object with the elements of 't', 'y', 'w', 'Tn' (option)
#' and 'ylu', returned by \code{check_input}.
#' @param brks A list object with the elements of 'fit' and 'dt', returned by
#' \code{season} or \code{season_3y}, which contains the growing season
#' dividing information.
#' @param nextent Extend curve fitting window, until \code{nextent} good or
#' marginal element are found in previous and subsequent growing season.
#' @param maxExtendMonth Search good or marginal good values in previous and
#' subsequent `maxExtendMonth` period.
#' @param minExtendMonth Extending perid defined by \code{nextent} and
#' \code{maxExtendMonth} should be no shorter than \code{minExtendMonth}.
#' When all points of the input time-series are good value, then the extending
#' period will be too short. In that situation, we can't make sure the connection
#' between different growing seasons is smoothing.
#' @param minT Double, use night temperature Tn to define backgroud value.
#' Tn < minT is treated as ungrowing season.
#' @param methods Character, find curve fitting names, can be one of
#' c("AG", "zhang", "beck", "elmore", "Gu").
#' @param qc Factor (optional), only suit for MOD13A1. SummaryQA code for PhenoExtracr.
#' Other dataset, just leave qc as the default.
#' @param minPercValid If valid percentage is less than \code{minPercValid}, the
#' fits are set to NA.
#' @param print Whether to print progress information?
#' @param ... Other parameters will be ignore.
#'
#' @return fits Multiple phenofit object.
#' @export
curvefits <- function(INPUT, brks,
                      wFUN = wTSM, iters = 2, wmin = 0.2,
                      nextent = 2, maxExtendMonth = 3, minExtendMonth = 1,
                      minT = 0,
                      methods = c('AG', 'zhang', 'beck', 'elmore', 'Gu'),
                      qc = INPUT$w*0+1, minPercValid = 0.2,
                      print = TRUE, ...)
{
    nptperyear <- INPUT$nptperyear
    t     <- INPUT$t
    years <- year(t)
    n    <- length(t)
    doys <- as.numeric(difftime(t, t[1], units = "day") + 1)#days from origin

    # Tn for background module
    Tn     <- INPUT$Tn #if has no Tn, NULL will be return
    has_Tn <- ifelse(is_empty(Tn), F, T)
    y0     <- INPUT$y0 #original y
    w <- w0 <- INPUT$w

    if (is.null(y0)) y0 <- INPUT$y
    if (all(is.na(INPUT$y))) return(NULL)

    # possible snow or cloud, replaced with whittaker smooth.
    I_w <- match(brks$whit$t, t) %>% rm_empty()

    I_fix <- which(w[I_w] == wmin)
    I_y   <- I_w[I_fix]
    INPUT$y[I_y] <- dplyr::last(brks$whit)[I_fix]
    w[I_w] <- brks$whit %>% {.[, contain(., "witer"), with = F]} %>% last()
    #w[I_fix]       <- wmin + 0.1 # exert the function of whitaker smoother

    di <- brks$di
    if (is.null(di)){
        getDateId <- function(dates) match(dates, t) #%>% rm_empty()
        getDateId_before <- function(dates) sapply(dates, function(date) last(which(t <= date)))
        getDateId_after  <- function(dates) sapply(dates, function(date) first(which(t >= date)))

        di <- data.table( beg  = getDateId_before(brks$dt$beg),
                          peak = getDateId_before(brks$dt$peak),
                          end  = getDateId_after(brks$dt$end)) %>% na.omit()
    }

    # plot(y, type = "b"); grid()
    # lines(brks$whit$iter3, col = "blue")
    # lines(INPUT$y        , col = "red")
    MaxExtendWidth = ceiling(nptperyear/12*maxExtendMonth)
    MinExtendWidth = ceiling(nptperyear/12*minExtendMonth)
    width_ylu    = nptperyear*2

    y    <- INPUT$y
    fits <- list()
    for (i in 1:nrow(di)){ #
        if (print) runningId(i, prefix = '\t[curvefits] ')

        I    <- di$beg[i]:di$end[i]
        I_beg <- di$beg[i]
        I_end <- di$end[i]

        I_extend <- get_extentI(w0, MaxExtendWidth, MinExtendWidth, I_beg, I_end, nextent, wmin)

        ### 2. input data
        ti   <- doys[I_extend]
        yi   <- INPUT$y[I_extend]
        wi   <- w[I_extend]
        # yi_good <- yi[w0[I_extend] > wmin]

        ### update ylu in a three year moving window
        ylu <- get_ylu (y, years, w0, width_ylu, I_beg:I_end, Imedian = TRUE, wmin)
        ylu <- merge_ylu(INPUT$ylu, ylu)

        # yi[yi < ylu[1]] <- ylu[1] # update y value

        # original weights, put in w0 incurvefitting is unwisdom, but for plot
        # w0   <- qc[I_extend] #INPUT$w
        # add background module here, 20180513
        if (has_Tn){
            Tni        <- Tn[I_extend]
            back_value <- backval(yi, ti, wi, Tni, minT, nptperyear)
            if (!is.na(back_value)){
                I_back     <- yi < back_value
                yi[I_back] <- back_value
                wi[I_back] <- 0.5
            }
        }
        beginI = ifelse(i == 1, 1, 2) # make sure no overlap
        tout <- doys[I] %>% {.[beginI]:last(.)} # make sure return the same length result.

        fit  <- curvefit(yi, ti, tout, nptperyear = nptperyear,
                         w = wi, ylu = ylu, iters = iters,
                         methods = methods, meth = 'BFGS', wFUN = wFUN, ...)
        # add original input data here, global calculation can comment this line
        data <- data.table(y = y0[I], t = doys[I], w = qc[I]) #INPUT$w[I]
        for (j in seq_along(fit)) fit[[j]]$data <- data
        # x <- fit$ELMORE
        # plot(y~t, x$data, type = "b"); grid()
        # lines(x$tout, x$fits$iter2)

        #if too much missing values
        if (sum(wi > pmax(wmin+0.1, 0.2))/length(wi) < minPercValid){
            fit %<>% map(function(x){
                x$fits %<>% map(~.x*NA) # list()
                return(x)
            })
        }
        fits[[i]] <- fit
    }
    # L1:curve fitting method, L2:yearly flag
    fits %<>% set_names(brks$dt$flag) %>% purrr::transpose()
    return(list(tout = t[first(di$beg):last(di$end)],  #dates for OUTPUT curve fitting VI
                fits = fits))
}

# extend curve fitting period
get_extentI <- function(w0, MaxExtendWidth, MinExtendWidth, I_beg, I_end, nextent = 1, wmin = 0.2){
    n <- length(w0)

    I_beg2  <- I_end2 <- NA
    # period <- floor(nptperyear/12*extend_month)
    if (!is.null(nextent)){
        I_beg_raw <- seq(max(1, I_beg-1), max(1, I_beg-MaxExtendWidth))
        I_end_raw <- seq(min(n, I_end+1), min(n, I_end+MaxExtendWidth))

        # at least `nextent` marginal point in previous and following season
        # I_beg2 and I_end2 have been constrained in the range of [1, n]
        I_beg2 <- I_beg_raw[ which(w0[I_beg_raw] > wmin)[nextent]  ]
        I_end2 <- I_end_raw[ which(w0[I_end_raw] > wmin)[nextent]  ]
    }

    # in case of previous and subsequent season good values too closed
    max_Beg <- max(1, I_beg - MinExtendWidth)
    min_End <- min(n, I_end + MinExtendWidth)

    I_beg2 <- ifelse( is.na(I_beg2), max_Beg, min(I_beg2, max_Beg))
    I_end2 <- ifelse( is.na(I_end2), min_End, max(I_end2, min_End))

    return( I_beg2:I_end2 )
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

    I_win  <- pmax(1, I_beg - width) : pmin(n, I_end + width)
    w0_win <- w0[I_win]
    I_win  <- I_win[which(w0_win > wmin)]

    if (is_empty(I_win)){
        ylu_max <- ylu_min <- NA
    } else{
        y_win  <- y[I_win]

        ylu_min <- min(y_win)
        ylu_max <- max(y_win)

        if (Imedian){
            year_win  <- years[I_win]
            # Assume peak in the middle of growing season, a few steps will be
            # ylu_min; while a long way to ylu_max; 20180918
            # length(I) reflects nptperyear
            if (width > length(I)*2/12){
                ylu_min <- aggregate(y_win, list(year = year_win), min)$x %>% median()
            }

            if (width > length(I)*7/12){
                ylu_max <- aggregate(y_win, list(year = year_win), max)$x %>% median()
            }
        }
    }
    c(ylu_min, ylu_max) # return
}
