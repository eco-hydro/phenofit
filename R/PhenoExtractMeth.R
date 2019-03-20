#' @name PhenoExtractMeth
#' @title Phenology Extraction methods
#'
#' @inheritParams D
#'
#' @param fFIT \code{fFIT} object returned by \code{\link{optim_pheno}}.
#' @param approach to be used to calculate phenology metrics.
#' 'White' (White et al. 1997) or 'Trs' for simple threshold.
#' @param trs threshold to be used for approach "Trs", in (0, 1).
#' @param IsPlot whether to plot?
#' @param show.lgd whether show figure lelend?
#' @param ... other parameters to PhenoPlot
#'
#' @examples
#' library(phenofit)
#' # simulate vegetation time-series
#' fFUN = doubleLog.Beck
#' par  = c(
#'     mn  = 0.1,
#'     mx  = 0.7,
#'     sos = 50,
#'     rsp = 0.1,
#'     eos = 250,
#'     rau = 0.1)
#' t    <- seq(1, 365, 8)
#' tout <- seq(1, 365, 1)
#' y <- fFUN(par, t)
#'
#' methods <- c("AG", "Beck", "Elmore", "Gu", "Zhang") # "Klos" too slow
#' fFITs <- curvefit(y, t, tout, methods)
#' fFIT  <- fFITs$fFIT$AG
#'
#' par(mfrow = c(2, 2))
#' PhenoTrs(fFIT)
#' PhenoDeriv(fFIT)
#' PhenoGu(fFIT)
#' PhenoKl(fFIT)
NULL

#' @description
#' \itemize{
#' \item \code{PhenoTrs} Threshold method
#' \item \code{PhenoDeriv} Derivative method
#' \item \code{PhenoGu} Gu method
#' \item \code{PhenoKl} Inflection method
#' }
#'
#' @param asymmetric If true, background value in spring season and autumn season
#' is regarded as different.
#'
#' @rdname PhenoExtractMeth
#' @export
PhenoTrs <- function(fFIT, approach = c("White", "Trs"), trs = 0.5, #, min.mean = 0.1
    asymmetric = TRUE,
    IsPlot = TRUE, ...)
{
    metrics <- c(sos = NA, eos = NA)

    t      <- fFIT$tout
    values <- last(fFIT$zs)
    n      <- length(t)

    # get peak of season position
    half.season <- median(which.max(values)) %>% round() # + 20, half season + 20 was unreasonable
    pop <- t[half.season]

    if (all(is.na(values))) return(metrics)
    if (half.season < 5 || half.season > (n - 5)) return(metrics)

    # get statistical values
    # n    <- t[length(t)]
    # avg  <- mean(x, na.rm = TRUE)

    # avg2 <- mean(x2[x2 > min.mean], na.rm = TRUE)
    peak <- max(values, na.rm = TRUE)

    if (asymmetric) {
        mn_a <- min(values[1:half.season], na.rm = T)
        mn_b <- min(values[-(1:half.season)], na.rm = T)

        mn <- c(rep(mn_a, half.season),
            rep(mn_b, n - half.season))
    } else {
        mn   <- min(values, na.rm = TRUE)
    }

    ampl <- peak - mn

    # select (or scale) values and thresholds for different methods
    approach <- approach[1]
    if (approach == "White") {
        # scale annual time series to 0-1
        ratio   <- (values - mn)/ampl
        # trs   <- 0.5
        trs.low <- trs - 0.05
        trs.up  <- trs + 0.05
    }
    if (approach == "Trs") {
        ratio   <- values
        a       <- diff(range(ratio, na.rm = TRUE)) * 0.05
        trs.low <- trs - a
        trs.up  <- trs + a
    }

    greenup <- .Greenup(ratio)
    ## distinguish the first half year and second year
    # select time where SOS and EOS are located (around trs value)
    bool <- ratio >= trs.low & ratio <= trs.up

    # get SOS, EOS, LOS
    # fixed 2017-01-04, according to TP phenology property
    # sos <- round(median(sose os[greenup & bool], na.rm = TRUE))
    # eos <- round(median(soseos[!greenup & bool], na.rm = TRUE))

    sos <- round(median(t[ greenup & bool & t < pop], na.rm = TRUE))
    eos <- round(median(t[!greenup & bool & t > pop], na.rm = TRUE))
    # los <- eos - sos#los < 0 indicate that error
    # los[los < 0] <- n + (eos[los < 0] - sos[los < 0])

    # get MGS, MSP, MAU
    # mgs <- mean(x[ratio > trs], na.rm = TRUE)
    # msp <- mau <- NA
    # if (!is.na(sos)) {
    #     id <- (sos - 10):(sos + 10)
    #     id <- id[(id > 0) & (id < n)]
    #     msp <- mean(x[which(index(x) %in% id == TRUE)], na.rm = TRUE)
    # }
    # if (!is.na(eos)) {
    #     id <- (eos - 10):(eos + 10)
    #     id <- id[(id > 0) & (id < n)]
    #     mau <- mean(x[which(index(x) %in% id == TRUE)], na.rm = TRUE)
    # }
    # metrics <- c(sos = sos, eos = eos, los = los, pop = pop, mgs = mgs,
    #   rsp = NA, rau = NA, peak = peak, msp = msp, mau = mau)
    metrics <- c(sos = sos, eos = eos)#, los = los

    if (IsPlot){
        main   <- ifelse(all(par("mar") == 0), "", sprintf("TRS%d", trs*10))
        PhenoPlot(t, values, main = main, ...)

        lines(t, trs*ampl + mn, lwd = linewidth)
        lines(t, trs.low*ampl + mn, lty = 2, lwd = linewidth)
        lines(t, trs.up*ampl + mn, lty = 2, lwd = linewidth)

        abline(v = metrics, col = colors[c(1, 4)], lwd = linewidth)
        text(metrics[1] - 5, min(trs + 0.15, 1)*ampl[1] + mn[1], "SOS", col = colors[1], adj = c(1, 0))
        text(metrics[2] + 5, min(trs + 0.15, 1)*last(ampl) + last(mn), "EOS", col = colors[4], adj = c(0, 0))
    }
    return(metrics)
    ### The function returns a vector with SOS, EOS, LOS, POP, MGS, rsp, rau, PEAK, MSP and MAU. }
}

# identify greenup or dormancy(brown) period
#' @export
.Greenup <- function(x, ...) {
    ratio.deriv <- c(NA, diff(x))
    greenup     <- rep(NA, length(x))
    greenup[ratio.deriv > 0] <- TRUE
    greenup[ratio.deriv < 0] <- FALSE
    return(greenup)
}

#' @inheritParams D
#' @inheritParams PhenoTrs
#'
#' @rdname PhenoExtractMeth
#' @export
PhenoDeriv <- function(fFIT,
    analytical = TRUE, smoothed.spline = FALSE,
    IsPlot = TRUE, show.lgd = TRUE, ...)
{
    PhenoNames <- c("SOS", "POP", "EOS")
    metrics <- setNames(rep(NA, 3), c("sos", "pop", "eos")) # template

    t      <- fFIT$tout
    values <- last(fFIT$zs)
    n      <- length(t)

    # get peak of season position
    half.season <- median(which.max(values)) # deal with multiple pop values
    pop <- t[half.season]

    if (all(is.na(values))) return(metrics)
    if (half.season < 5 || half.season > (n - 5)) return(metrics)

    der1   <- D1.fFIT(fFIT, analytical, smoothed.spline)
    # get SOS and EOS according to first order derivative
    # fixed 20180510, ±5 to make sure sos and eos are not in the POP.
    # I_sos <- median(which.max(der1[1:(half.season - 5)]))
    # I_eos <- median(which.min(der1[(half.season+5):length(der1)])) + half.season

    # der.sos (eos) is impossible to occur near the pop.
    nmin  <- 5
    I_sos <- findpeaks( der1[1:(half.season - 5)], nups=nmin,
        ndowns=nmin, npeaks=1, sortstr=TRUE)$X$pos
    I_eos <- findpeaks(-der1[(half.season+5):length(der1)],
        nups=nmin, ndowns=nmin, npeaks=1, sortstr=TRUE)$X$pos + half.season + 4

    # if half.season > length(der1), error will be occur
    sos <- t[I_sos]
    eos <- t[I_eos]
    if (is_empty(sos)) sos <- NA
    if (is_empty(eos)) eos <- NA

    metrics <- c(sos = sos, pop = pop, eos = eos)#, los = los

    if (IsPlot){
        main  <- ifelse(all(par("mar") == 0), "", "DER")
        PhenoPlot(t, values, main = main, ...)
        if (show.lgd) legend('topright', c("f(t)'"), lty = 2, col = "black", bty='n')

        abline(v = c(sos, eos), col = colors[c(1, 4)], lwd = linewidth)
        abline(v = pop, col ="darkgreen", lty = 1, lwd = linewidth)

        A <- diff(range(values))
        I_metrics <- match(metrics, t)
        if (all(is.na(I_metrics))) {
            ylons <- I_metrics
        }else{
            ylons <- values[I_metrics] + c(1, -1, 1)*0.1*A
        }
        xlons <- metrics + c(-1, 1, 1)*5
        xlons[xlons < min(t)] <- min(t)
        xlons[xlons > max(t)] <- max(t)

        I <- c(1); text(xlons[I], ylons[I], PhenoNames[I], col = colors[I], adj = c(1, 0))
        I <- 2:3 ; text(xlons[I], ylons[I], PhenoNames[I], col = c("darkgreen", colors[3]), adj = c(0, 0))

        #der1 last plot
        op <- par(new = TRUE)
        plot(t, der1, type= "l", lty = 2, lwd = linewidth,
             col = "black", axes = FALSE)
    }
    return(metrics)
}


#' @inheritParams PhenoTrs
#' @importFrom dplyr last
#' @rdname PhenoExtractMeth
#' @export
PhenoGu <- function(fFIT,
    analytical = TRUE, smoothed.spline = FALSE,
    IsPlot = TRUE, ...)
{
    PhenoNames <- c("UD", "SD", "DD", "RD")
    metrics <- setNames(rep(NA, 4), c("UD", "SD", "DD", "RD"))

    t      <- fFIT$tout
    values <- last(fFIT$zs)
    n      <- length(t)

    # get peak of season position
    half.season <- median(which.max(values)) # deal with multiple pop values
    pop <- t[half.season]

    if (all(is.na(values))) return(metrics)
    if (half.season < 5 || half.season > (n - 5)) return(metrics)

    der1   <- D1.fFIT(fFIT, analytical, smoothed.spline)
    # get SOS and EOS according to first order derivative
    sos.index <- median(which.max(der1[1:(half.season-5)]))
    eos.index <- median(which.min(der1[(half.season+5):length(der1)])) + half.season
    sos <- t[sos.index]
    eos <- t[eos.index]

    rsp <- der1[sos.index] # rate of spring, also known as peak recovery rate
    rau <- der1[eos.index] # rate of autumn, peak senescence rate
    VI.rsp <- values[sos.index] # VI index of corresponding date
    VI.rau <- values[eos.index]

    # Gu Phenology Extraction method also quite rely on first order derivative
    rl.b <- VI.rsp - rsp*sos # interception of recovery line
    sl.b <- VI.rau - rau*eos # interception of recovery line

    baseline <- min(values, na.rm = TRUE)
    maxline  <- max(values, na.rm = TRUE)

    ## y = kx + b; x = (y - b)/k
    UD <- (baseline - rl.b)/rsp # upturn day, intersection of rl and x axis
    SD <- (maxline  - rl.b)/rsp # stabilization day, intersection of maxline and rl
    DD <- (maxline  - sl.b)/rau # downturn day, intersection of maxline and sl
    RD <- (baseline - sl.b)/rau # recession day, intersection of sl and x axis

    ## subset data between SD and DD
    sub.time <- t[t >= SD & t <= DD]
    sub.gcc  <- values[t >= SD & t <= DD]

    ## compute a linear fit
    if (length(sub.time) > 3) {
        X <- cbind(1, sub.time)
        plateau.lm        <- .lm.fit(X, sub.gcc)
        plateau.slope     <- plateau.lm$coefficients[2]
        plateau.intercept <- plateau.lm$coefficients[1]
        # y1 = rau*t + sl.b
        # y2 = k2t + b2
        # so, t = (sl.b - b2)/(k2 -rau)
        DD <- (sl.b - plateau.intercept) / (plateau.slope - rau)
    }else{
        plateau.slope     <- NA
        plateau.intercept <- NA
    }
    ## calculate area under the curve
    # cut.x <- days[which(days>=UD & days<=RD)]
    # cut.y <- offset.y[which(days>=UD & days<=RD)]
    # the.fun <- function(t) {eval(retrieved.formula, envir=as.list(params))}
    metrics   <- round(c(UD, SD, DD, RD)) %>%
        sapply(function(x){ ifelse(x < min(t) || x > max(t), NA, x) }) %>%
        set_names(PhenoNames)
    # c("UD", "SD", "DD", "RD", "maxline", "baseline", "rsp", "rau", "plateau.slope")
    # c(pheno, maxline, baseline, rsp, rau, plateau.slope)

    if (IsPlot){
        main   <- ifelse(all(par("mar") == 0), "", "Gu")
        PhenoPlot(t, values, main = main, ...)

        abline(rl.b, rsp, col = "blue", lty=  2, lwd = linewidth,)
        abline(sl.b, rau, col = "red" , lty=  2, lwd = linewidth,)
        # any na values input to abline will lead to error
        if (all(!is.na(c(plateau.slope, plateau.intercept))))
            abline(a = plateau.intercept, b = plateau.slope,
                lty = 2, lwd = linewidth, col = "darkgreen")

        abline(h = c(maxline, baseline), lty = 2, lwd = linewidth,)
        abline(v = metrics[1:4],  col = colors, lwd = linewidth,)

        A <- diff(range(values))
        I_metrics <- match(metrics, t)
        if (all(is.na(I_metrics))) {
            ylons <- I_metrics
        }else{
            ylons <- values[I_metrics] + c(1, -3, -3, 1)*0.1*A
        }

        xlons <- metrics[1:4] + c(-1, -1, 1, 1)*5
        xlons[xlons < min(t)] <- min(t)
        xlons[xlons > max(t)] <- max(t)
        I <- c(1, 2); text(xlons[I], ylons[I], PhenoNames[I], col = colors[I], adj = c(1, 0))
        I <- c(3, 4); text(xlons[I], ylons[I], PhenoNames[I], col = colors[I], adj = c(0, 0))
    }
    return(metrics)
}

#' @inheritParams PhenoTrs
#' @rdname PhenoExtractMeth
#' @export
PhenoKl <- function(fFIT,
    analytical = TRUE, smoothed.spline = FALSE, 
    IsPlot = TRUE, show.lgd = TRUE, ...)
{
    PhenoNames <- c("Greenup", "Maturity", "Senescence", "Dormancy")
    metrics <- setNames(rep(NA, 4), PhenoNames)

    t      <- fFIT$tout
    values <- last(fFIT$zs)
    n      <- length(t)
    xlim   <- range(t)

    # get peak of season position
    half.season <- median(which.max(values)) # + 20, half season + 20 was unreasonable
    pop <- t[half.season]

    if (all(is.na(values))) return(metrics)
    if (half.season < 5 || half.season > (n - 5)) return(metrics)

    derivs <- curvature.fFIT(fFIT, analytical, smoothed.spline)
    k      <- derivs$k
    # define cutoff date for spline functions

    # x <- ifelse(uncert==TRUE, x$uncertainty, x$fit)
    # if (length(which(is.na(k) == TRUE)) != 0 | length(which(is.infinite(k) == TRUE)) != 0) {

    # If have no NA and infinite values
    if (!(any(is.na(k)) || any(is.infinite(k)))) {
        spline.k <- smooth.spline(k, df = 0.1 * length(k))
        der.k    <- predict(spline.k, d = 1)$y
        der.k2   <- predict(spline.k, d = 2)$y

        # der.k  <- c(NA, diff(k))
        # der.k2 <- c(NA, NA, diff(k, differences = 2))
        ## find maxima of derivative of k ## split season
        dist_fromPeak <- 1 # days
        asc.k   <- try(der.k[1:(half.season - dist_fromPeak)]) # k   of first half year
        asc.k.d <- try(t[1:(half.season - dist_fromPeak)])
        # doy of first half year
        des.k   <- try(der.k[(half.season + dist_fromPeak):length(k)])
        des.k.d <- try(t[(half.season + dist_fromPeak):length(k)])

        # first half maximum local values of k'
        pos <- findpeaks(asc.k, minpeakdistance = 15, npeaks = 2, sortstr = TRUE)$X$pos
        pos <- sort(pos)
        pos <- c(rep(NA, 2 - length(pos)), pos) #at least two values
        I_asc <- asc.k.d[pos]
        if (all(is.na(pos))){
            I_asc <- pos
        }else{
            I_asc <- asc.k.d[pos]
        }

        # second half minimum local values of k'
        pos <- findpeaks(-des.k, minpeakdistance = 15, npeaks = 2, sortstr = TRUE)$X$pos
        pos <- sort(pos);
        pos <- c(pos, rep(NA, 2 - length(pos)))
        if (all(is.na(pos))){
            I_dec <- pos
        }else{
            I_dec <- des.k.d[pos]
        }
        metrics <- c(I_asc, I_dec)
    }

    if (IsPlot){
        main   <- ifelse(all(par("mar") == 0), "", "Zhang (Curvature Rate)")

        A          <- diff(range(der.k))
        I_metrics <- match(metrics, t)
        if (all(is.na(I_metrics))) {
            ylons <- I_metrics
        }else{
            ylons <- values[I_metrics] + c(1, -1, -1, 1)*0.2*A
        }
        xlons  <- metrics + c(1, -1, 1, -1) * 5
        xlons[xlons < min(t)] <- min(t)
        xlons[xlons > max(t)] <- max(t)
        # plotrix::twoord.plot(t, values, t, der.k,
        #                      main = main, type= c("p", "b"), xlim = xlim,
        #                      rcol = "black", lcol = "grey60", lpch = 20, rpch = 1,
        #                      rytickpos = NULL)
        PhenoPlot(t, values, main = main, ...)
        if (show.lgd){
            legend('topright', c("K'"), lty = c(3), col = c("black"), bty='n') ##pch =c(20, 1),
        }

        pop     <- t[half.season]
        abline(v = metrics, col = colors, lwd = linewidth)
        # abline(v = pop, col ="darkgreen", lty = 1, lwd = linewidth)
        # abline(v = pop + 20, col ="darkgreen", lty = 2, lwd = linewidth)
        I <- c(1, 3); text(xlons[I], ylons[I], PhenoNames[I], col = colors[I], adj = c(0, 0))
        I <- c(2, 4); text(xlons[I], ylons[I], PhenoNames[I], col = colors[I], adj = c(1, 0))
        # der.k last plot
        op <- par(new = TRUE)
        plot(t, der.k, xlim = xlim, type= "l",
             lty = 3, lwd = linewidth, col = "black", axes = TRUE) # cex = 1, pch = 20,
    }
    return(setNames(metrics, PhenoNames))
}
