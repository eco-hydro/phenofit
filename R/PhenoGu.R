#' Phenology extraction in GU method (GU)
#'
#' @inheritParams PhenoDeriv 
#' @inherit PhenoTrs examples
#' @param ... other parameters to [PhenoGu.default()] or [PhenoGu.fFIT()]
#' 
#' @references
#' 1. Gu, L., Post, W. M., Baldocchi, D. D., Black, T. A., Suyker, A. E., Verma,
#'    S. B., … Wofsy, S. C. (2009). Characterizing the Seasonal Dynamics of
#'    Plant Community Photosynthesis Across a Range of Vegetation Types. In A.
#'    Noormets (Ed.), Phenology of Ecosystem Processes: Applications in Global
#'    Change Research (pp. 35–58). New York, NY: Springer New York.
#'    \doi{10.1007/978-1-4419-0026-5_2}
#' 2. Filippa, G., Cremonese, E., Migliavacca, M., Galvagno, M., Forkel, M.,
#'    Wingate, L., … Richardson, A. D. (2016). Phenopix: A R package for
#'    image-based vegetation phenology. Agricultural and Forest Meteorology,
#'    220, 141–150. \doi{10.1016/j.agrformet.2016.01.006}
#' 
#' @return A numeric vector, with the elements of:
#' - `UD`: upturn date
#' - `SD`: stabilisation date
#' - `DD`: downturn date
#' - `RD`: recession date 
#' @export
PhenoGu <- function(x, t, ...) UseMethod("PhenoGu", x)

#' @rdname PhenoGu
#' @export
PhenoGu.fFIT <- function(x, t = NULL,
    analytical = TRUE, smoothed.spline = FALSE, ...) 
{
    if (!is.null(x$tout)) t <- x$tout
    values <- last2(x$zs)
    der1 <- D1.fFIT(x, t, analytical, smoothed.spline)
    PhenoGu.default(values, t, der1, ...)
}

#' @inheritParams PhenoTrs
#' @rdname PhenoGu
#' @export
PhenoGu.default <- function(x, t, der1,
    IsPlot = TRUE, ...)
{
    PhenoNames <- c("UD", "SD", "DD", "RD")
    metrics <- setNames(rep(NA, 4), c("UD", "SD", "DD", "RD"))

    n      <- length(t)

    # get peak of season position
    half.season <- median(which.max(x)) # deal with multiple pos x
    pos <- t[half.season]

    if (all(is.na(x))) return(metrics)
    if (half.season < 5 || half.season > (n - 5)) return(metrics)

    # get SOS and EOS according to first order derivative
    sos.index <- median(which.max(der1[1:(half.season-5)]))
    eos.index <- median(which.min(der1[(half.season+5):length(der1)])) + half.season
    sos <- t[sos.index]
    eos <- t[eos.index]

    rsp <- der1[sos.index] # rate of spring, also known as peak recovery rate
    rau <- der1[eos.index] # rate of autumn, peak senescence rate
    VI.rsp <- x[sos.index] # VI index of corresponding date
    VI.rau <- x[eos.index]

    # Gu Phenology Extraction method also quite rely on first order derivative
    rl.b <- VI.rsp - rsp*sos # interception of recovery line
    sl.b <- VI.rau - rau*eos # interception of recovery line

    baseline <- min(x, na.rm = TRUE)
    maxline  <- max(x, na.rm = TRUE)

    ## y = kx + b; x = (y - b)/k
    UD <- (baseline - rl.b)/rsp # upturn day, intersection of rl and x axis
    SD <- (maxline  - rl.b)/rsp # stabilization day, intersection of maxline and rl
    DD <- (maxline  - sl.b)/rau # downturn day, intersection of maxline and sl
    RD <- (baseline - sl.b)/rau # recession day, intersection of sl and x axis

    ## subset data between SD and DD
    sub.time <- t[t >= SD & t <= DD]
    sub.gcc  <- x[t >= SD & t <= DD]

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

    ## Avoid NA x by using `clamp` instead
    # sapply(function(x){ ifelse(x < min(t) || x > max(t), NA, x) }) %>%
    metrics   <- round(c(UD, SD, DD, RD)) %>%
        sapply(clamp, lims = c(t[1], t[n])) %>%
        set_names(PhenoNames)
    # c("UD", "SD", "DD", "RD", "maxline", "baseline", "rsp", "rau", "plateau.slope")
    # c(pheno, maxline, baseline, rsp, rau, plateau.slope)

    if (IsPlot){
        main   <- ifelse(all(par("mar") == 0), "", "Gu")
        PhenoPlot(t, x, main = main, ...)

        abline(rl.b, rsp, col = "blue", lty=  2, lwd = linewidth,)
        abline(sl.b, rau, col = "red" , lty=  2, lwd = linewidth,)
        # any na x input to abline will lead to error
        if (all(!is.na(c(plateau.slope, plateau.intercept))))
            abline(a = plateau.intercept, b = plateau.slope,
                lty = 2, lwd = linewidth, col = "darkgreen")

        abline(h = c(maxline, baseline), lty = 2, lwd = linewidth,)
        abline(v = metrics[1:4],  col = colors, lwd = linewidth,)

        A <- diff(range(x))
        I_metrics <- match(metrics, t)
        if (all(is.na(I_metrics))) {
            ylons <- I_metrics
        }else{
            ylons <- x[I_metrics] + c(1, -3, -3, 1)*0.1*A
        }

        xlons <- metrics[1:4] + c(-1, -1, 1, 1)*5
        xlons[xlons < min(t)] <- min(t)
        xlons[xlons > max(t)] <- max(t)

        I <- c(1, 2); text(xlons[I], ylons[I], PhenoNames[I], col = colors[I], adj = c(1, 0))
        I <- c(3, 4); text(xlons[I], ylons[I], PhenoNames[I], col = colors[I], adj = c(0, 0))
    }
    return(metrics)
}
