#' @title Rough fitting
#' @name roughFit
#' 
#' @description
#' Divide growing seasons according to rough fitting (`rFUN`) result .
#' 
#' For `season`, rough fitting is applied for whole.
#' For `season_mov` rough fitting is applied in every year, during which
#' `maxExtendMonth` is extended.
#' 
#' @details
#' Before division growing season, `INPUT` should be added a year in head
#' and tail first by `add_HeadTail`.
#'
#' Finally, use [findpeaks()] to get local maximum and local minimum values.
#' Two local minimum define a growing season.
#' If two local minimum(maximum) are too closed, then only the smaller(biger)
#' is left.
#'
#' @param INPUT A list object with the elements of `t`, `y`, `w`,
#' `Tn` (optional) and `ylu`, returned by [check_input()].
#' @param options see details
#' @param  ... ignored.
#' 
#' @inheritSection season_mov options for season
#' @return
#' - `fit`: A data.table of Rough fitting result, with the columns of
#' (`t`, `y`, `witer1`, ..., `witerN`, `ziter1`, ..., `ziterN`).
#'
#' - `dt`: A data.table of growing season division information, with the columns
#' of (`beg`, `peak`, `end`, `y_beg`, `y_peak`, `y_end`, `len`, `year`,
#' `season`, `flag`).
#'
#' @seealso [findpeaks()].
#' @example R/examples/ex-season.R
#'
#' @keywords internal
#' @export
roughFit <- function(INPUT, options = list(),
    frame = floor(INPUT$nptperyear / 5) * 2 + 1, ...)
{
    if (all(is.na(INPUT$y))) return(NULL)

    south <- INPUT$south
    t <- INPUT$t
    y <- INPUT$y
    nlen <- length(t)

    set_options(season = options)
    .options$season %<>% modifyList(list(frame = frame))
    opt <- .options$season

    c(nyear, years) %<-% guess_nyear(INPUT)
    nptperyear <- INPUT$nptperyear
    list2env(opt[c("lambda", "nf", "frame")], environment())

    if (is.null(frame)) frame <- floor(nptperyear / 5) * 2 + 1
    if (is.null(lambda)) lambda <- max(nyear * frame, 15)

    ylu0 <- INPUT$ylu
    A0 <- diff(ylu0)
    ylu <- ylu0
    rFUN = check_function(opt$rFUN)

    nups <- default_nups(nptperyear)
    for (iloop in 1:3) {
        param <- c(
            INPUT[c("y", "t", "w", "ylu", "nptperyear")],
            opt[c("wFUN", "wmin", "iters")], listk(lambda, nf, frame)
        )
        yfits <- do.call(opt$rFUN, param)
        ypred <- last(yfits$zs) # as.numeric(runmed(ypred, frame))
        alpha <- 0.01
        # default is three year data input, median will be much better
        # This module was primarily designed for `season_mov`. It also works for
        # group large than 3-year. 2-year median will be underestimated.
        if (nyear >= 2.5) { # considering NA values, nyear of 3-year will be smaller.
            b <- get_ylu.default(ypred, t) # boundary
            # If multiple years are completely missing, ylu_min possiable equals ylu_max
            if (b$A > A0 * 0.2) {
                ylu <- c(pmax(b$ylu_min, INPUT$ylu[1]), # quantile(ypred, alpha/2)
                        pmin(b$ylu_max, INPUT$ylu[2]))
            }
        }
        INPUT$ylu <- ylu
        A <- diff(ylu)
        # minPeakHeight <- pmax(ypeak_min, A*0.1 + ylu[1])
        info_peak <- findpeaks_season(ypred,
            opt$r_max, opt$r_min,
            minpeakdistance = 0, minpeakheight = opt$ypeak_min,
            nups = nups, nyear = nyear)

        if (.options$season$verbose) {
            cat(sprintf(
                "iloop = %d: lambda = %.1f, ntrough_PerYear = %.2f, npeak_PerYear = %.2f\n",
                iloop, lambda, info_peak$ntrough_PerYear, info_peak$npeak_PerYear
            ))
        }
        if (!opt$adj.param) break
        pars <- adjustRoughParam(lambda, nf, frame, info_peak, nptperyear, opt$MaxPeaksPerYear, opt$MaxTroughsPerYear)
        if (pars$status == FALSE) break
        list2env(pars[c("lambda", "nf", "frame")], environment())
    }
    rfit <- as.data.table(c(list(t = t, y = y), yfits$ws, yfits$zs))
    # years is used get calendar break points
    info_peak %<>% c(listk(year = years, nptperyear, south))
    listk(fit = rfit, info_peak)
}

## This module will automatically update lambda, nf and wHANTS
adjustRoughParam <- function(lambda, nf, frame,
                             info_peak, nptperyear, MaxPeaksPerYear, MaxTroughsPerYear) {
    npeak_PerYear <- info_peak$npeak_PerYear
    ntrough_PerYear <- info_peak$ntrough_PerYear

    # auto_adjust <- function(nptperyear, npeak_PerYear, ntrough_PerYear, MaxPeaksPerYear, MaxTroughsPerYear, ) {
    delta_frame <- ceiling(nptperyear / 12)
    # adjust frame in the step of `month`
    status <- TRUE
    if (is.null(info_peak$npeak_PerYear) || is.null(ntrough_PerYear)) {
        status <- FALSE
    } else {
        if (npeak_PerYear > MaxPeaksPerYear | ntrough_PerYear > MaxTroughsPerYear) {
            lambda <- lambda * 2
            nf <- max(1, nf - 1)
            frame <- min(frame + delta_frame, nptperyear * 2)
        } else if (npeak_PerYear < 0.8 | ntrough_PerYear < 0.8) {
            lambda <- lambda / 2
            nf <- min(5, nf + 1)
            frame <- max(frame - delta_frame, delta_frame)
        } else {
            status <- FALSE
        }
    }
    listk(lambda, nf, frame, status = status)
}

#' get rough fitting
#'
#' @param brks returned by function [season_mov()]
#'
#' @keywords internal
#' @return
#' - `data`:
#'     + t
#'     + y
#'     + QC_flag
#' - `tout`:
#' - `zs`: list of iter1, ..., itern
#' - `ws`: list of iter1, ..., itern
#' @export
brks2rfit <- function(brks) {
    dt = brks$dt
    fit = brks$fit
    data = fit[, .(t, y)]

    # doys = difftime(data$t, data$t[1])
    t = data$t #%>% as.integer()
    tout = seq(data$t[1], last(data$t), by = "day") #%>% as.integer()
    zs = ldply(fit %>% select(starts_with("ziter")),
        ~approx2(t, .x, tout, na.rm = FALSE)$y)
    structure(listk(data, dt, tout, zs), class = "rfit")
}

approx2 <- function(...) suppressWarnings(approx(...))
