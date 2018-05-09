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
    sos.index <- median(which.max(der1[1:half.season]))
    eos.index <- median(which.min(der1[(half.season+1):length(der1)])) + half.season
    # if half.season > length(der1), error will be occur
    sos <- t[sos.index]
    eos <- t[eos.index]

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
    sos.index <- median(which.max(der1[1:half.season]))
    eos.index <- median(which.min(der1[(half.season+1):length(der1)])) + half.season
    sos <- t[sos.index]
    eos <- t[eos.index]

    rsp <- der1[sos.index] # rate of spring, also known as peak recovery rate
    rau <- der1[eos.index] # rate of autumn, peak senescence rate
    VI.rsp <- values[sos.index] # VI index of corresponding date
    VI.rau <- values[eos.index]

    # Gu Phenology Extraction method also quite rely on first order derivative
    rl.y  <- rsp * (t - sos) + VI.rsp
    rl.eq <- lm(rl.y ~ t)                        # recovery line

    sl.y  <- rau * (t - eos) + VI.rau
    sl.eq <- lm(sl.y ~ t)                     # senenscence line

    baseline <- min(values, na.rm = T)
    maxline  <- max(values, na.rm = T)

    ## y = kx + b; x = (y - b)/k
    #  upturn day is the intersection between rl and x axis
    UD <- (baseline - rl.eq$coefficients[[1]])/rl.eq$coefficients[[2]]
    #  stabilization day, intersection between maxline and rl
    SD <- (maxline  - rl.eq$coefficients[[1]])/rl.eq$coefficients[[2]]
    #  downturn day, intersection between maxline and sl
    DD <- (maxline  - sl.eq$coefficients[[1]])/sl.eq$coefficients[[2]]
    #  recession day is the intersection between sl and x axis
    RD <- (baseline - sl.eq$coefficients[[1]])/sl.eq$coefficients[[2]]

    ## subset data between SD and DD
    sub.time <- t[t >= SD & t <= DD]
    sub.gcc  <- values[t >= SD & t <= DD]

    ## compute a linear fit
    if (length(sub.time) > 3) {
        plateau.lm   <- lm(sub.gcc ~ sub.time)
        M            <- matrix(c(coef(plateau.lm)[2], coef(sl.eq)[2], -1, -1), nrow = 2, ncol = 2)
        intercepts   <- as.matrix(c(coef(plateau.lm)[1], coef(sl.eq)[1]))
        interception <- -solve(M) %*% intercepts
        DD           <- interception[1, 1]
        plateau.slope     <- plateau.lm$coefficients[2]
        plateau.intercept <- plateau.lm$coefficients[1]
    }else{
        plateau.slope     <- NA
        plateau.intercept <- NA
    }
    ## calculate area under the curve
    # cut.x <- days[which(days>=UD & days<=RD)]
    # cut.y <- offset.y[which(days>=UD & days<=RD)]
    # the.fun <- function(t) {eval(retrieved.formula, envir=as.list(params))}
    #
    metrics   <- round(c(UD, SD, DD, RD)) %>%
        sapply(function(x){ ifelse(x < min(t) || x > max(t), NA, x) }) %>%
        setNames(., c("UD", "SD", "DD", "RD"))
    # c("UD", "SD", "DD", "RD", "maxline", "baseline", "rsp", "rau", "plateau.slope")
    # c(pheno, maxline, baseline, rsp, rau, plateau.slope)

    if (IsPlot){
        main   <- ifelse(all(par("mar") == 0), "", "Gu")
        PhenoPlot(t, values, main = main)

        abline(rl.eq, col = "blue", lty=  2)
        abline(sl.eq, col = "red" , lty=  2)
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
    values <- last(fit$fits)
    xlim   <- range(t)

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
        half.season    <- which.max(values)  # + 20, half season + 20 was unreasonable
        asc.k   <- try(der.k[1:half.season]) # k   of first half year
        asc.k.d <- try(t[1:half.season])
        # doy of first half year
        des.k   <- try(der.k[half.season:length(k)])
        des.k.d <- try(t[half.season:length(k)])

        I_asc <- asc.k.d[localMinMax(asc.k, 'max')]
        I_dec <- des.k.d[localMinMax(des.k, 'min')]
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
        xlons[xlons < min(t)]   <- min(t)
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

        pop     <- t[which.max(values)]

        abline(v = metrics, col = colors)
        abline(v = pop, col ="darkgreen", lty = 1)
        abline(v = pop + 20, col ="darkgreen", lty = 2)
        text(xlons[1:2], ylons[1:2], PhenoNames[1:2], col = colors[1:2], adj = c(1, 0))
        text(xlons[3:4], ylons[3:4], PhenoNames[3:4], col = colors[3:4], adj = c(0, 0))
    }
    return(setNames(metrics, PhenoNames))
}

#'
#' Local maxima and minima based on Tommy's solution. And only kept the best if
#' two points was so close.
#'
#' maxval: f' = 0, f'_left > 0 & f'_right < 0
#' minval: f' = 0, f'_left < 0 & f'_right > 0
#'
#' https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
#' If having NA values means phenology result was nquestionable.
#'
#' @param x numeric vector
#' @param meth A string parameter, only can be 'min' or 'max'.
#' @param MinSep The Minimum gap of two local minimum (maximum).
#'
#' @export
localMinMax <- function(x, meth = c("min", "max"), MinSep = 15) {
    if (length(x) < 3) return(c(NA, NA))
    meth <- meth[1]
    # 1. find the local values
    # Use -Inf instead if x is numeric (non-integer)
    val_begin <- ifelse(meth == "max", -1, 1) * .Machine$integer.max
    y         <- diff(c(val_begin, x)) > 0L #& abs(x) > 1e-6
    # rle(y)$lengths
    run <- rle(y)
    len <- cumsum(run$lengths)
    if (meth == "max"){
        I   <- len[run$values]
        fun <- which.max
        decreasing <- T
    }else if (meth == "min"){
        I   <- len[!run$values]
        fun <- which.min
        decreasing <- F
    }
    if (x[[1]] == x[[2]]) I <- I[-1]

    # plot(x)
    # abline(v = I, lty = 2)
    # I <- c(1, 4, 6, 7, 10, 20,24, 28, 35)

    #2. If two local extreme values too close, kept the best
    if (length(I) >= 2){
        run <- rle(diff(I) <= MinSep)
        len <- run$length
        Id <- c(0, cumsum(len)) # %>% {.[seq(1, length(.), 2)]}

        I_opt <- list()
        # I_dup <- which(run$values)
        for (i in seq_along(run$values)){
            Ii     <- (Id[i]+1):(Id[i+1]+1)
            # If diff <= MinSep, then in those values select the best
            if (run$values[i]) Ii <- Ii[ fun(x[I[Ii]]) ]
            # If diff > MinSep, then accept those values
            I_opt[[i]] <- Ii
        }
        I <- I[unlist(I_opt)]
    }

    I <- setdiff(I, c(1, length(x))) #remove lcoal minima and maxima value at tail and head
    if (length(I) > 2) {
        #if local values number large than 2, one the extremest first two kept
        I <- sort(I[order(x[I], decreasing = decreasing)[1:2]])
    }
    if (length(I) < 2) {
        # if local max value number less than 2, add NA at tail
        # if local min value number less than 2, add NA at head
        if (meth == "min"){
            I <- c(rep(NA, 2 - length(I)), I) #at least two values
        }else if (meth == "max"){
            I <- c(I, rep(NA, 2 - length(I))) #at least two values
        }
    }
    return(I)
}
