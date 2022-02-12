#' Phenology extraction in Derivative method (DER)
#' 
#' @inheritParams D
#' @inheritParams PhenoTrs
#' @param der1 the first order difference
#' 
#' @references 
#' 1. Filippa, G., Cremonese, E., Migliavacca, M., Galvagno, M., Forkel, M.,
#'    Wingate, L., … Richardson, A. D. (2016). Phenopix: A R package for
#'    image-based vegetation phenology. Agricultural and Forest Meteorology,
#'    220, 141–150. \doi{10.1016/j.agrformet.2016.01.006}
#' 
#' @seealso [PhenoTrs()], [PhenoGu()], [PhenoKl()]
#' @export
PhenoDeriv <- function(x, t, ...) UseMethod("PhenoDeriv", x)

#' @rdname PhenoDeriv
#' @export
PhenoDeriv.fFIT <- function(x, t = NULL,
    analytical = FALSE, smoothed.spline = FALSE, ...) {
    if (!is.null(x$tout)) t <- x$tout
    values <- last2(x$zs)
    der1 <- D1.fFIT(x, t, analytical, smoothed.spline)
    PhenoDeriv.default(values, t, der1, ...)
}

#' @rdname PhenoDeriv
#' @param show.legend whether show figure lelend?
#' @export
PhenoDeriv.default <- function(x, t, der1,
    IsPlot = TRUE, show.legend = TRUE, ...)
{
    PhenoNames <- c("SOS", "POS", "EOS")
    metrics <- setNames(rep(NA, 3), c("sos", "pos", "eos")) # template
    n      <- length(t)

    # get peak of season position
    half.season <- median(which.max(x)) # deal with multiple pos x
    pos <- t[half.season]

    if (all(is.na(x))) return(metrics)
    if (half.season < 5 || half.season > (n - 5)) return(metrics)

    # get SOS and EOS according to first order derivative
    # fixed 20180510, ±5 to make sure sos and eos are not in the POP.
    # I_sos <- median(which.max(der1[1:(half.season - 5)]))
    # I_eos <- median(which.min(der1[(half.season+5):length(der1)])) + half.season

    # der.sos (eos) is impossible to occur near the pop.
    nmin  <- 5
    I_sos <- findpeaks( der1[1:(half.season - 5)],
        nups=nmin, ndowns=nmin, npeaks=1, sortstr=TRUE)$X$pos
    I_eos <- findpeaks(-der1[(half.season+5):length(der1)],
        nups=nmin, ndowns=nmin, npeaks=1, sortstr=TRUE)$X$pos + half.season + 4

    # if half.season > length(der1), error will be occur
    sos <- t[I_sos]
    eos <- t[I_eos]
    if (is_empty(sos)) {
        I_sos <- findpeaks( der1[1:(half.season - 5)],
            nups=nmin, ndowns=0, npeaks=1, sortstr=TRUE)$X$pos
        sos <- t[I_sos]
        if (is_empty(sos)) sos <- round((t[1] + pos)/2)
    }
    if (is_empty(eos)) {
        I_eos <- findpeaks(-der1[(half.season+5):length(der1)],
            nups=nmin, ndowns=0, npeaks=1, sortstr=TRUE)$X$pos + half.season + 4
        eos <- t[I_eos]
        if (is_empty(eos)) eos <- round((t[n] + pos)/2)
    }
    # greenup <- .Greenup(x)
    # if (is.na(sos)) {
    #     if (all(!rm_empty(greenup))) sos = pos
    # }
    # if (is.na(eos)) {
    #     if (all(rm_empty(greenup)))  eos = pos
    # }
    metrics <- c(sos = sos, pos = pos, eos = eos)#, los = los

    if (IsPlot){
        main  <- ifelse(all(par("mar") == 0), "", "DER")
        PhenoPlot(t, x, main = main, ...)
        if (show.legend) legend('topright', c("f(t)'"), lty = 2, col = "black", bty='n')

        abline(v = c(sos, eos), col = colors[c(1, 4)], lwd = linewidth)
        abline(v = pos, col ="blue", lty = 1, lwd = linewidth)

        A <- diff(range(x))
        I_metrics <- match(metrics, t)
        if (all(is.na(I_metrics))) {
            ylons <- I_metrics
        }else{
            ylons <- x[I_metrics] + c(1, -2, 1)*0.1*A
        }
        xlons <- metrics + c(-1, 1, 1)*5
        xlons[xlons < min(t)] <- min(t)
        xlons[xlons > max(t)] <- max(t)

        I <- c(1); text(xlons[I], ylons[I], PhenoNames[I], col = colors[I], adj = c(1, 0))
        I <- 2:3 ; text(xlons[I], ylons[I], PhenoNames[I], col = c("blue", colors[4]), adj = c(0, 0))

        #der1 last plot
        op <- par(new = TRUE)
        plot(t, der1, type= "l", lty = 2, lwd = linewidth,
             col = "black", axes = FALSE)
    }
    return(metrics)
}
