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
#' @param methods Fine curve fitting methods, can be one or more of
#' \code{c('AG', 'Beck', 'Elmore', 'Gu', 'Klos', 'Zhang')}.
#' @param QC_flag Factor (optional) returned by \code{qcFUN}, levels should be
#' in the range of \code{c("snow", "cloud", "shadow", "aerosol", "marginal",
#' "good")}, others will be categoried into \code{others}. \code{QC_flag} is
#' used for visualization in \code{\link{PhenoExtract}} and
#' \code{\link{plot_phenofit}}.
#' @param minPercValid If the percentage of good and marginal quality points is
#' less than \code{minPercValid}, curve fiting result is set to \code{NA}.
#' @param print Whether to print progress information?
#' @param ... Other parameters will be ignore.
#'
#' @return fits Multiple phenofit object.
#'
#' @examples
#' library(phenofit)
#' data("MOD13A1")
#'
#' dt <- tidy_MOD13.gee(MOD13A1$dt)
#' st <- MOD13A1$st
#'
#' sitename <- dt$site[1]
#' d     <- dt[site == sitename, ] # get the first site data
#' sp    <- st[site == sitename, ] # station point
#' # global parameter
#' IsPlot = TRUE
#' print  = FALSE
#' nptperyear = 23
#' ypeak_min  = 0.05
#' wFUN = wTSM
#'
#' dnew     <- add_HeadTail(d, nptperyear = nptperyear) # add one year in head and tail
#' INPUT    <- check_input(dnew$t, dnew$y, dnew$w, nptperyear,
#'                         maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
#'
#' brks2 <- season_3y(INPUT,
#'     rFUN = wWHIT, wFUN = wFUN,
#'     plotdat = d, IsPlot = IsPlot, print = FALSE, IsOnlyPlotbad = FALSE)
#'
#' fit <- curvefits(
#'     INPUT, brks2,
#'     methods = c("AG", "Beck", "Elmore", "Zhang"), #,"klos", "Gu"
#'     verbose = FALSE,
#'     wFUN = wFUN,
#'     nextent = 2, maxExtendMonth = 2, minExtendMonth = 1,
#'     QC_flag = dnew$QC_flag, minPercValid = 0.2,
#'     print = TRUE)
#'
#' df_fit <- get_fitting(fit)
#' g <- plot_phenofit(df_fit, brks2)
#' grid::grid.newpage(); grid::grid.draw(g)
#' @export
curvefits <- function(INPUT, brks,
                      wFUN = wTSM, iters = 2, wmin = 0.2,
                      nextent = 2, maxExtendMonth = 3, minExtendMonth = 1,
                      minT = 0,
                      methods = c('AG', 'Beck', 'Elmore', 'Gu', 'Klos', 'Zhang'),
                      QC_flag = NULL,
                      minPercValid = 0.2,
                      print = TRUE, ...)
{
    if (all(is.na(INPUT$y))) return(NULL)

    nptperyear <- INPUT$nptperyear
    t     <- INPUT$t
    years <- year(t)
    n     <- length(t)

    # doys  <- as.numeric(t) # days from origin
    doys <- as.numeric(difftime(t, date.origin, units = "day") + 1)

    # Tn for background module
    w <- w0 <- INPUT$w
    y0     <- INPUT$y0 # original y
    Tn     <- INPUT$Tn # if has no Tn, NULL will be return
    has_Tn <- ifelse(is_empty(Tn), F, TRUE)

    # possible snow or cloud, replaced with Whittaker smoothing.
    I_w <- match(brks$whit$t, t) %>% rm_empty()

    I_fix <- which(w[I_w] == wmin)
    I_y   <- I_w[I_fix]
    INPUT$y[I_y] <- dplyr::last(brks$whit)[I_fix]
    w[I_w] <- brks$whit %>% {.[, contain(., "witer"), with = F]} %>% last()
    # w[I_fix] <- wmin + 0.1 # exert the function of whitaker smoother

    # growing season dividing
    getDateId <- function(dates) match(dates, t) #%>% rm_empty()
    getDateId_before <- function(dates) sapply(dates, function(date) last(which(t <= date)))
    getDateId_after  <- function(dates) sapply(dates, function(date) first(which(t >= date)))

    di <- data.table( beg  = getDateId_before(brks$dt$beg),
                      peak = getDateId_before(brks$dt$peak),
                      end  = getDateId_after(brks$dt$end)) %>% na.omit()

    MaxExtendWidth = ceiling(nptperyear/12*maxExtendMonth)
    MinExtendWidth = ceiling(nptperyear/12*minExtendMonth)
    width_ylu      = nptperyear*2

    y    <- INPUT$y
    fits <- list()
    for (i in 1:nrow(di)){
        if (print) runningId(i, prefix = '\t[curvefits] ')

        I     <- di$beg[i]:di$end[i]
        I_beg <- di$beg[i]
        I_end <- di$end[i]

        I_extend <- get_extentI(w0, MaxExtendWidth, MinExtendWidth, I_beg, I_end, nextent, wmin)

        ## 2. input data
        ti   <- doys[I_extend]
        yi   <- INPUT$y[I_extend]
        wi   <- w[I_extend]
        # yi_good <- yi[w0[I_extend] > wmin]

        ## update ylu in a three year moving window
        ylu <- get_ylu (y, years, w0, width_ylu, I_beg:I_end, Imedian = TRUE, wmin)
        ylu <- merge_ylu(INPUT$ylu, ylu)
        # yi[yi < ylu[1]] <- ylu[1] # update y value

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
        beginI <- ifelse(i == 1, 1, 2) # make sure no overlap
        tout   <- doys[I] %>% {.[beginI]:last(.)} # make sure return the same length result.

        fFITs  <- curvefit(yi, ti, tout, nptperyear = nptperyear,
                         w = wi, ylu = ylu, iters = iters,
                         methods = methods, wFUN = wFUN, ...)
        # add original input data here, global calculation can comment this line

        # \code{y} at here is original time-series without checked
        data <- list(t = doys[I], y = y0[I], QC_flag = QC_flag[I]) %>% as.data.table()
        fFITs$data <- data

        # if too much missing values
        if (sum(wi > pmax(wmin+0.1, 0.2))/length(wi) < minPercValid){
            fFITs %<>% map(function(x){
                x$zs %<>% map(~.x*NA) # list()
                return(x)
            })
        }
        fits[[i]] <- fFITs
    }
    # L1:curve fitting method, L2:yearly flag
    fits %<>% set_names(brks$dt$flag) # %>% purrr::transpose()
    fits
    # return(list(tout = t[first(di$beg):last(di$end)],  # dates for OUTPUT curve fitting VI
    #             fits = fits))
}


############################## END OF CURVEFITS ################################
# HIDING FUNCTIONS
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
        I_beg2 <- I_beg_raw[ which(w0[I_beg_raw] > wmin)[nextent] ]
        I_end2 <- I_end_raw[ which(w0[I_end_raw] > wmin)[nextent] ]
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
