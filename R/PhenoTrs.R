# colors    <- c("blue", "green3", "orange", "red")
colors <- c("green3", "darkgreen", "darkorange3", "red")
linewidth <- 1.2


# identify greenup or dormancy(brown) period
.Greenup <- function(x, ...) {
    ratio.deriv <- c(NA, diff(x))
    greenup <- rep(NA, length(x))
    greenup[ratio.deriv > 0] <- TRUE
    greenup[ratio.deriv < 0] <- FALSE
    return(greenup)
}

#' Phenology extraction in Threshold method (TRS)
#' 
#' @param x numeric vector, or `fFIT` object returned by [curvefit()].
#' @param t `doy` vector, corresponding doy of vegetation index.
#' @param approach to be used to calculate phenology metrics.
#' 'White' (White et al. 1997) or 'Trs' for simple threshold.
#' @param trs threshold to be used for approach "Trs", in (0, 1).
#' @param IsPlot whether to plot?
#' @param ... other parameters to PhenoPlot
#' @param asymmetric If true, background value in spring season and autumn season
#' is regarded as different.
#' 
#' @example R/examples/ex-PhenoTrs.R
#' @seealso [PhenoDeriv()], [PhenoGu()], [PhenoKl()]
#' @export
PhenoTrs <- function(x, t = NULL, 
    approach = c("White", "Trs"), trs = 0.5, #, min.mean = 0.1
    asymmetric = TRUE, IsPlot = TRUE, ...) 
{
    UseMethod("PhenoTrs", x)
}

#' @rdname PhenoTrs
#' @export
PhenoTrs.fFIT <- function(x, t = NULL, ...) {
    if (!is.null(x$tout)) t <- x$tout
    values <- last2(x$zs)
    PhenoTrs.default(values, t, ...)
}

#' @export
#' @rdname PhenoTrs
PhenoTrs.default <- function(x, t = NULL,
                             approach = c("White", "Trs"), trs = 0.5, # , min.mean = 0.1
                             asymmetric = TRUE,
                             IsPlot = TRUE, ...) {
    metrics <- c(sos = NA, eos = NA)
    n <- length(x)

    # get peak of season position
    half.season <- median(which.max(x)) %>% round() # + 20, half season + 20 was unreasonable
    pos <- t[half.season]

    if (all(is.na(x))) {
        return(metrics)
    }
    if (half.season < 5 || half.season > (n - 5)) {
        return(metrics)
    }

    # get statistical x
    # avg  <- mean(x, na.rm = TRUE)

    # avg2 <- mean(x2[x2 > min.mean], na.rm = TRUE)
    peak <- max(x, na.rm = TRUE)

    if (asymmetric) {
        mn_a <- min(x[1:half.season], na.rm = T)
        mn_b <- min(x[-(1:half.season)], na.rm = T)

        mn <- c(
            rep(mn_a, half.season),
            rep(mn_b, n - half.season)
        )
    } else {
        mn <- min(x, na.rm = TRUE)
    }

    ampl <- peak - mn

    # select (or scale) x and thresholds for different methods
    approach <- approach[1]
    if (approach == "White") {
        # scale annual time series to 0-1
        ratio <- (x - mn) / ampl
        # trs   <- 0.5
        trs.low <- trs - 0.05
        trs.up <- trs + 0.05
    }
    if (approach == "Trs") {
        ratio <- x
        a <- diff(range(ratio, na.rm = TRUE)) * 0.05
        trs.low <- trs - a
        trs.up <- trs + a
    }

    greenup <- .Greenup(ratio)
    ## distinguish the first half year and second year
    # select time where SOS and EOS are located (around trs value)
    bool <- ratio >= trs.low & ratio <= trs.up

    # get SOS, EOS, LOS
    # fixed 2017-01-04, according to TP phenology property
    # sos <- round(median(sose os[greenup & bool], na.rm = TRUE))
    # eos <- round(median(soseos[!greenup & bool], na.rm = TRUE))

    sos <- round(median(t[greenup & bool & t < pos], na.rm = TRUE))
    eos <- round(median(t[!greenup & bool & t > pos], na.rm = TRUE))

    # greenup <- .Greenup(ratio)
    if (is.na(sos + eos)) {
        bool2 <- ratio >= (trs - 0.1) & ratio <= (trs + 0.1)
        if (is.na(sos)) sos <- round(median(t[greenup & bool2 & t < pos], na.rm = TRUE))
        if (is.na(eos)) eos <- round(median(t[!greenup & bool2 & t > pos], na.rm = TRUE))

        # if still is NA
        if (is.na(sos)) sos <- round((t[1] + pos) / 2)
        if (is.na(eos)) eos <- round((t[n] + pos) / 2)
    }
    # monotonous = length(unique(greenup[2:n]) %>% rm_empty()) == 1
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
    # metrics <- c(sos = sos, eos = eos, los = los, pos = pos, mgs = mgs,
    #   rsp = NA, rau = NA, peak = peak, msp = msp, mau = mau)
    metrics <- c(sos = sos, eos = eos) # , los = los

    if (IsPlot) {
        main <- ifelse(all(par("mar") == 0), "", sprintf("TRS%d", trs * 10))
        PhenoPlot(t, x, main = main, ...)

        lines(t, trs * ampl + mn, lwd = linewidth)
        # lines(t, trs.low*ampl + mn, lty = 2, lwd = linewidth)
        # lines(t, trs.up*ampl + mn, lty = 2, lwd = linewidth)
        abline(v = metrics, col = colors[c(1, 4)], lwd = linewidth)
        text(metrics[1] - 5, min(trs + 0.15, 1) * ampl[1] + mn[1], "SOS", col = colors[1], adj = c(1, 0))
        text(metrics[2] + 5, min(trs + 0.15, 1) * last(ampl) + last(mn), "EOS", col = colors[4], adj = c(0, 0))
    }
    return(metrics)
    ### The function returns a vector with SOS, EOS, LOS, POS, MGS, rsp, rau, PEAK, MSP and MAU. }
}
