colors <- c("blue", "green3", "orange", "red")

#' PhenoPlot
#' @examples
#' PhenoPlot(t, x, main = main)
#' @export
PhenoPlot <- function(t, x, main = "", ...){
    plot(t, x, main = main, ...,
             type= "p", cex = 1.4, col = "grey60", pch = 20) #pch = 20,
    grid()
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

#' @export
PhenoTrs <- function(fit, approach = c("White", "Trs"), trs = 0.5, min.mean = 0.1,
    IsPlot = TRUE, ...) {
    t      <- fit$tout
    values <- last(fit$fits)

    if (all(is.na(values))) return(c(sos = NA, eos = NA))

    # get statistical values
    # n    <- t[length(t)]
    # avg  <- mean(x, na.rm = TRUE)
    x2   <- na.omit(values)
    # avg2 <- mean(x2[x2 > min.mean], na.rm = TRUE)
    peak <- max(x2)
    mn   <- min(x2)
    ampl <- peak - mn

    # get peak of season position
    pop <- median(t[which.max(values)])
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
    # if (plot) {
    #     if (approach == 'White') PlotPhenCycle(x, metrics=metrics, trs=trs, ...)
    #     if (approach == 'Trs') PlotPhenCycle(ratio, metrics=metrics, trs=trs, ...)
    # }
    if (IsPlot){
        main   <- ifelse(all(par("mar") == 0), "", sprintf("TRS%d", trs*10))
        PhenoPlot(t, ratio, main = main)

        abline(h = trs)
        abline(h = c(trs.low, trs.up), lty = 2)
        abline(v = metrics, col = colors[c(1, 4)])
        text(metrics + c(-1, 1)*24, min(trs + 0.15, 1), c("SOS", "EOS"), col = colors)
    }
    return(metrics)
    ### The function returns a vector with SOS, EOS, LOS, POP, MGS, rsp, rau, PEAK, MSP and MAU. }
}

#'
#' PhenoDeriv
#' @description rechecked in 2017-11-17, with no problems
#' @export
PhenoDeriv <- function(fit, IsPlot = TRUE, smspline = TRUE, ...){
    t      <- fit$tout
    values <- last(fit$fits)

    if (all(is.na(values))) return( setNames(rep(NA, 3), c("sos", "pop", "eos")) )

    der1   <- D1.phenofit(fit, smspline = smspline)

    # get peak of season position
    half.season <- median(which.max(values)) #deal with multiple pop values
    pop <- t[half.season]

    # get SOS and EOS according to first order derivative
    # fixed 20180510, ±5 to make sure sos and eos are not in the POP.
    # I_sos <- median(which.max(der1[1:(half.season - 5)]))
    # I_eos <- median(which.min(der1[(half.season+5):length(der1)])) + half.season

    nmin  <- 5
    I_sos <- findpeaks( der1[1:(half.season - 5)]         , nups=nmin, ndowns=nmin, npeaks=1, sortstr=TRUE)$X$pos
    I_eos <- findpeaks(-der1[(half.season+5):length(der1)], nups=nmin, ndowns=nmin, npeaks=1, sortstr=TRUE)$X$pos + half.season

    # if half.season > length(der1), error will be occur
    sos <- t[I_sos]
    eos <- t[I_eos]
    if (is_empty(sos)) sos <- NA
    if (is_empty(eos)) eos <- NA

    if (IsPlot){
        main  <- ifelse(all(par("mar") == 0), "", "DER")
        PhenoPlot(t, values, main = main)
        op <- par(new = T)
        plot(t, der1, type= "b", cex = 1, col = "black", axes = FALSE, pch = 16) #pch = 20,
        legend('bottomleft', c("Fitting", "f'"), lty = c(0, 1), pch =c(20, 1),
            col = c("grey60", "black"), bty='n')
        abline(v = c(sos, eos), col = colors[c(1, 4)])
    }
    metrics <- c(sos = sos, pop = pop, eos = eos)#, los = los
    return(metrics)
}

#' PhenoGu
#' @export
#' @importFrom dplyr last
PhenoGu <- function(fit, IsPlot = TRUE, smspline = TRUE, ...) {
    t      <- fit$tout
    values <- last(fit$fits)

    metrics <- setNames(rep(NA, 4), c("UD", "SD", "DD", "RD"))
    if (all(is.na(values))) return(metrics)

    der1   <- D1.phenofit(fit, smspline = smspline)
    # get peak of season position
    half.season <- median(which.max(values)) #deal with multiple pop values
    pop <- t[half.season]
    if (half.season >= length(der1)) return(metrics)

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
    
    baseline <- min(values, na.rm = T)
    maxline  <- max(values, na.rm = T)

    ## y = kx + b; x = (y - b)/k
    #  upturn day, intersection of rl and x axis
    UD <- (baseline - rl.b)/rsp
    #  stabilization day, intersection of maxline and rl
    SD <- (maxline  - rl.b)/rsp
    #  downturn day, intersection of maxline and sl
    DD <- (maxline  - sl.b)/rau
    #  recession day, intersection of sl and x axis
    RD <- (baseline - sl.b)/rau

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
        setNames(., c("UD", "SD", "DD", "RD"))
    # c("UD", "SD", "DD", "RD", "maxline", "baseline", "rsp", "rau", "plateau.slope")
    # c(pheno, maxline, baseline, rsp, rau, plateau.slope)

    if (IsPlot){
        main   <- ifelse(all(par("mar") == 0), "", "Gu")
        PhenoPlot(t, values, main = main)

        abline(rl.b, rsp, col = "blue", lty=  2)
        abline(sl.b, rau, col = "red" , lty=  2)
        # any na values input to abline will lead to error
        if (all(!is.na(c(plateau.slope, plateau.intercept))))
            abline(a = plateau.intercept, b = plateau.slope, lty = 2, col = "darkgreen")

        abline(h = c(maxline, baseline), lty = 2)
        abline(v = metrics[1:4],  col = colors)

        A <- diff(range(values))
        ylons  <- c(values[metrics[1]] + 0.1*A,
                    values[metrics[2]] - 0.3*A,
                    values[metrics[2]] - 0.3*A, #pred[p3] - 0.3*A,
                    values[metrics[4]] + 0.1*A)
        text(metrics[1:4] + c(-1, 1, -1, 1)*15, ylons, c("UD", "SD", "DD", "RD"),
            col = colors, offset = 0)
    }
    return(metrics)
}

#' @export
PhenoKl <- function(fit, IsPlot = TRUE, ...) {
    PhenoNames <- c("Greenup", "Maturity", "Senescence", "Dormancy")
    metrics <- setNames(rep(NA, 4), PhenoNames)

    t      <- fit$tout
    xlim   <- range(t)
    values <- last(fit$fits)
    half.season <- which.max(values)  # + 20, half season + 20 was unreasonable
    
    if (all(is.na(values))) return(setNames(rep(NA, 4), PhenoNames))

    derivs <- curvature.phenofit(fit, smspline = TRUE)
    k      <- derivs$k
    # define cutoff date for spline functions

    # x <- ifelse(uncert==TRUE, x$uncertainty, x$fit)
    # if (length(which(is.na(k) == TRUE)) != 0 | length(which(is.infinite(k) == TRUE)) != 0) {

    # If have no NA and infinite values
    if (!(any(is.na(k)) || any(is.infinite(k)))) {
        spline.k <- smooth.spline(k, df = 0.1 * length(k))
        der.k    <- predict(spline.k, d = 1)$y
        der.k2   <- predict(spline.k, d = 2)$y

        ## find maxima of derivative of k ## split season
        asc.k   <- try(der.k[1:(half.season-5)]) # k   of first half year
        asc.k.d <- try(t[1:(half.season-5)])
        # doy of first half year
        des.k   <- try(der.k[(half.season+5):length(k)])
        des.k.d <- try(t[(half.season+5):length(k)])

        # first half minimum local values of k'
        pos <- findpeaks(-asc.k, minpeakdistance = 15, npeaks = 2, sortstr = TRUE)$X$pos
        pos <- sort(pos)
        pos <- c(rep(NA, 2 - length(pos)), pos) #at least two values
        I_asc <- asc.k.d[pos]
        if (all(is.na(pos))){
            I_asc <- pos
        }else{
            I_asc <- asc.k.d[pos]
        }

        # second half maximum local values of k'
        pos <- findpeaks(des.k, minpeakdistance = 15, npeaks = 2, sortstr = TRUE)$X$pos
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

        A      <- diff(range(der.k))
        ylons  <- c(der.k[metrics[1]] + 0.1*A,
                    der.k[metrics[2]] + 0.1*A,
                    der.k[metrics[3]] - 0.1*A,
                    der.k[metrics[4]] - 0.1*A)

        xlons  <- metrics + c(-1, -1, 1, 1) * 5
        xlons[xlons < min(t)] <- min(t)
        xlons[xlons > max(t)] <- max(t)
        # plotrix::twoord.plot(t, values, t, der.k,
        #                      main = main, type= c("p", "b"), xlim = xlim,
        #                      rcol = "black", lcol = "grey60", lpch = 20, rpch = 1,
        #                      rytickpos = NULL)
        PhenoPlot(t, values, main = main, xlim = xlim)
        op <- par(new = T)
        plot(t, der.k, xlim = xlim, type= "b", cex = 1, col = "black", pch = 16, axes = T) #pch = 20,
        legend('bottomleft', c("Fitting", "K'"), lty = c(0, 1), pch =c(20, 1),
               col = c("grey60", "black"), bty='n')

        pop     <- t[half.season]

        abline(v = metrics, col = colors)
        abline(v = pop, col ="darkgreen", lty = 1)
        abline(v = pop + 20, col ="darkgreen", lty = 2)
        text(xlons[1:2], ylons[1:2], PhenoNames[1:2], col = colors[1:2], adj = c(1, 0))
        text(xlons[3:4], ylons[3:4], PhenoNames[3:4], col = colors[3:4], adj = c(0, 0))
    }
    return(setNames(metrics, PhenoNames))
}
