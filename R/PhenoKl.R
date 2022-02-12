
#' Phenology extraction in Inflection method (Zhang)
#'
#' @inheritParams PhenoDeriv
#' @param fFIT object return by [curvefit()]
#'
#' @inherit PhenoTrs examples
#' @references
#' 1. Zhang, X., Friedl, M. A., Schaaf, C. B., Strahler, A. H., Hodges, J. C. F.
#'    F., Gao, F., … Huete, A. (2003). Monitoring vegetation phenology using
#'    MODIS. Remote Sensing of Environment, 84(3), 471–475.
#'    \doi{10.1016/S0034-4257(02)00135-9}
#' @return A numeric vector, with the elements of:
#' `Greenup`, `Maturity`, `Senescence`, `Dormancy`.
#'
#' @export
PhenoKl <- function(fFIT, t = NULL,
                    analytical = FALSE, smoothed.spline = FALSE,
                    IsPlot = TRUE, show.legend = TRUE, ...) {
    PhenoNames <- c("Greenup", "Maturity", "Senescence", "Dormancy")
    metrics <- setNames(rep(NA, 4), PhenoNames)

    if (!is.null(fFIT$tout)) t <- fFIT$tout
    x <- last2(fFIT$zs)
    n <- length(t)
    xlim <- range(t)

    # get peak of season position
    half.season <- median(which.max(x)) # + 20, half season + 20 was unreasonable
    pos <- t[half.season]

    if (all(is.na(x))) {
        return(metrics)
    }
    if (half.season < 5 || half.season > (n - 5)) {
        return(metrics)
    }

    derivs <- curvature.fFIT(fFIT, t, analytical, smoothed.spline)
    k <- derivs$k
    # define cutoff date for spline functions

    # x <- ifelse(uncert==TRUE, x$uncertainty, x$fit)
    # if (length(which(is.na(k) == TRUE)) != 0 | length(which(is.infinite(k) == TRUE)) != 0) {

    # If have no NA and infinite x
    if (!(any(is.na(k)) || any(is.infinite(k)))) {
        spline.k <- smooth.spline(k, df = 0.1 * length(k))
        der.k <- predict(spline.k, d = 1)$y
        der.k2 <- predict(spline.k, d = 2)$y

        # der.k  <- c(NA, diff(k))
        # der.k2 <- c(NA, NA, diff(k, differences = 2))
        ## find maxima of derivative of k ## split season
        dist_fromPeak <- 0 # days
        asc.k <- try(der.k[1:(half.season - dist_fromPeak)]) # k   of first half year
        asc.k.d <- try(t[1:(half.season - dist_fromPeak)])
        # doy of first half year
        des.k <- try(der.k[(half.season + dist_fromPeak):length(k)])
        des.k.d <- try(t[(half.season + dist_fromPeak):length(k)])

        A <- range(der.k, na.rm = TRUE) %>% diff()
        # first half maximum local x of k'
        # minpeakdistance in the unit of day
        pos <- findpeaks(asc.k,
            npeaks = 2, ndowns = 0, sortstr = TRUE,
            minpeakdistance = 15, minpeakheight = A * 0.025 + min(asc.k)
        )$X$pos
        pos <- sort(pos)
        # pos <- c(rep(NA, 2 - length(pos)), pos) # at least two x
        pos <- fill_NA_kl(pos, length(asc.k))

        I_asc <- asc.k.d[pos]
        if (all(is.na(pos))) {
            I_asc <- pos
        } else {
            I_asc <- asc.k.d[pos]
        }

        # second half minimum local x of k'
        pos <- findpeaks(-des.k,
            npeaks = 2, sortstr = TRUE,
            minpeakdistance = 15, minpeakheight = A * 0.025 + min(-des.k)
        )$X$pos
        pos <- sort(pos)
        # pos <- c(pos, rep(NA, 2 - length(pos)))
        pos <- fill_NA_kl(pos, length(des.k))
        
        if (all(is.na(pos))) {
            I_dec <- pos
        } else {
            I_dec <- des.k.d[pos]
        }
        metrics <- c(I_asc, I_dec)
    }

    if (IsPlot) {
        main <- ifelse(all(par("mar") == 0), "", "Zhang (Curvature Rate)")
        A <- diff(range(x))

        I_metrics <- match(metrics, t)
        if (all(is.na(I_metrics))) {
            ylons <- I_metrics
        } else {
            ylons <- x[I_metrics] + c(-1, -1, -1, 1) * 0.1 * A
        }

        xlons <- metrics + c(1, -1, 1, -1) * 5
        xlons[xlons < min(t)] <- min(t)
        xlons[xlons > max(t)] <- max(t)
        # plotrix::twoord.plot(t, x, t, der.k,
        #                      main = main, type= c("p", "b"), xlim = xlim,
        #                      rcol = "black", lcol = "grey60", lpch = 20, rpch = 1,
        #                      rytickpos = NULL)
        PhenoPlot(t, x, main = main, ...)
        if (show.legend) {
            legend("topright", c("K'"), lty = c(3), col = c("black"), bty = "n") ## pch =c(20, 1),
        }

        pos <- t[half.season]
        abline(v = metrics, col = colors, lwd = linewidth)
        # abline(v = pos, col ="darkgreen", lty = 1, lwd = linewidth)
        # abline(v = pos + 20, col ="darkgreen", lty = 2, lwd = linewidth)
        PhenoNames2 <- c("Greenup", "Maturity", "Senescence", "Dormancy")
        # PhenoNames2 <- c("G", "M", "S", "D")

        I <- c(1, 3)
        text(xlons[I], ylons[I], PhenoNames2[I], col = colors[I], adj = c(0, 0))
        I <- c(2, 4)
        text(xlons[I], ylons[I], PhenoNames2[I], col = colors[I], adj = c(1, 0))
        # der.k last plot
        op <- par(new = TRUE)
        plot(t, der.k,
            xlim = xlim, type = "l",
            lty = 3, lwd = linewidth, col = "black", axes = TRUE
        ) # cex = 1, pch = 20,
    }
    return(setNames(metrics, PhenoNames))
}


# not sure in the current result
fill_NA_kl <- function(pos, half.season) {
    # guess the position of NA value
    if (length(pos) == 1) {
        fill_right = pos < (half.season - pos)

        if (fill_right) {
            # if pos near the pos
            pos = c(pos, rep(NA, 2 - length(pos)))
        } else {
            # if pos near the end
            pos = c(rep(NA, 2 - length(pos)), pos)
        }
    } else {
        pos <- c(pos, rep(NA, 2 - length(pos)))
    }
    pos
}
