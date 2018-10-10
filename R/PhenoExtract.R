colors    <- c("blue", "green3", "orange", "red")
linewidth <- 1.2

#' PhenoPlot
#'
#' @inheritParams check_input
#' @param main figure title
#' @param ... ignored parameters
#'
#' @export
PhenoPlot <- function(t, y, main = "", ...){
    plot(t, y, main = main, ...,
             type= "l", cex = 2, col = "black", lwd = linewidth) #pch = 20,
    # grid(nx = NA)
    grid(ny = 4, nx = NA)
}

#' ExtractPheno
#'
#' Get yearly vegetation phenological metrics of a curve fitting method
#'
#' @param fits Multiple \code{phenofit} object.
#' @param TRS Threshold for \code{PhenoTrs}.
#' @param IsPlot Boolean. Whether to plot figure?
#'
#' @export
ExtractPheno <- function(fits, TRS = c(0.1, 0.2, 0.5, 0.6), IsPlot = FALSE){
    names <- names(fits)
    pheno_list <- list()
    methods    <- c(paste0("TRS", TRS*10),"DER","GU", "ZHANG")
    TRS_last   <- last(TRS) # only display last threshold figure

    if (IsPlot)
        op <- par(mfrow = c(length(fits), 5),
            oma = c(1, 2, 3, 1), mar = rep(0, 4), yaxt = "n", xaxt = "n")
    ylim <- NULL
    for (i in seq_along(fits)){
        fit    <- fits[[i]]
        ypred  <- last(fit$fits)
        all_na <- all(is.na(ypred))
        # 1. show curve fitting RMSE
        # need to fix here, about status variable. 31 Jan, 2018
        show.lgd = FALSE
        if (IsPlot && !all_na){
            ti <- fit$data$t
            yi <- fit$data$y # may have NA values.
            # constrain plot ylims
            ylim0    <- c( pmin(min(yi, na.rm = T), min(ypred)),
                           pmax(max(yi, na.rm = T), max(ypred)))
            A = diff(ylim0);
            ylim     <- ylim0 + c(-1, 0.2) * 0.05 *A
            ylim_trs <- (ylim - ylim0) / A # TRS:0-1

            PhenoPlot(fit$tout, ypred, ylim = ylim)
            lines(ti, yi, lwd = 1, col = "grey60")
            # pch = 19, col = "grey60"
            wi <- as.numeric(fit$data$w)

            ## Just designed for MOD13A1
            # Levels:  good  margin  snow&ice  cloud
            labels <- c(" good", " margin", " snow&ice", " cloud")
            colors <- c("grey60", "#00BFC4", "#F8766D", "#C77CFF")
            pch <- c(19, 15, 4, 17)
            for (j in 1:4){
                ind = which(wi == j)
                if (!is_empty(ind)) points(ti[ind], yi[ind], pch = pch[j], col = colors[j])
            }

            if (i == 1){
                show.lgd = TRUE
                legend('topright', c('y', "f(t)"), lty = c(1, 1), pch =c(1, NA), bty='n')
            }
            stat     <- statistic.phenofit(fit)
            stat_txt <- sprintf("  R=%.2f, p=%.3f\n RMSE=%.3f\nNSE=%.2f\n",
                                stat[['R']], stat[['pvalue']], stat[['RMSE']], stat[['NSE']])
            legend('topleft', stat_txt, adj = c(0.2, 0.2), bty='n', text.col = "red")
            mtext(names[i], side = 2)
        }
        if (i == 1 && IsPlot) mtext("fitting")

        p_TRS <- map(TRS, function(trs) {
            PhenoTrs(fit, approach = "White", trs = trs, IsPlot = FALSE)
        })

        if (IsPlot && !all_na) {
            trs6 <- PhenoTrs(fit, approach="White", trs=TRS_last, IsPlot = IsPlot, ylim = ylim)
            if (i == 1) mtext(sprintf('TRS%d', TRS_last*10))
        }

        param_common  <- list(fit, IsPlot, ylim = ylim)
        param_common2 <- list(fit, IsPlot, ylim = ylim, show.lgd = show.lgd)

        der   <- do.call(PhenoDeriv, param_common2);   if (i == 1 && IsPlot) mtext("DER")
        gu    <- do.call(PhenoGu, param_common)[1:4];  if (i == 1 && IsPlot) mtext("GU")
        zhang <- do.call(PhenoKl, param_common2);      if (i == 1 && IsPlot) mtext("ZHANG")
        pheno_list[[i]] <- c(p_TRS, list(der, gu, zhang)) %>% set_names(methods)
    }
    pheno_list %<>% set_names(names)
    return(pheno_list)
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

#' Phenology extraction functions
#'
#'
#' @inheritParams D1.phenofit
#' 
#' @param approach to be used to calculate phenology metrics. 
#' 'White' (White et al. 1997) or 'Trs' for simple threshold.
#' @param trs threshold to be used for approach "Trs", in (0, 1).
#' @param IsPlot whether to plot?
#' @param show.lgd whether show figure lelend?
#' @param ... other parameters to PhenoPlot
#' @param IsSmoothed Boolean. If false, positive derivative in spring and negative 
#' derivative will be not applied. 
#' 
#' @rdname PhenoExtractMeth
#' @export
PhenoTrs <- function(fit, approach = c("White", "Trs"), trs = 0.5, #, min.mean = 0.1
    IsPlot = TRUE, IsSmoothed = T, ...) {
    metrics <- c(sos = NA, eos = NA)

    t      <- fit$tout
    values <- last(fit$fits)
    n      <- length(t)

    # get peak of season position
    half.season <- median(which.max(values)) # + 20, half season + 20 was unreasonable
    pop <- t[half.season]

    if (all(is.na(values))) return(metrics)
    if (half.season < 5 || half.season > (n - 5)) return(metrics)

    # get statistical values
    # n    <- t[length(t)]
    # avg  <- mean(x, na.rm = TRUE)
    x2   <- na.omit(values)
    # avg2 <- mean(x2[x2 > min.mean], na.rm = TRUE)
    peak <- max(x2)
    mn   <- min(x2)
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
    # if (plot) {
    #     if (approach == 'White') PlotPhenCycle(x, metrics=metrics, trs=trs, ...)
    #     if (approach == 'Trs') PlotPhenCycle(ratio, metrics=metrics, trs=trs, ...)
    # }
    if (IsPlot){
        main   <- ifelse(all(par("mar") == 0), "", sprintf("TRS%d", trs*10))
        PhenoPlot(t, values, main = main, ...)

        abline(h = trs*ampl + mn, lwd = linewidth)
        abline(h = c(trs.low, trs.up)*ampl + mn, lty = 2, lwd = linewidth)
        abline(v = metrics, col = colors[c(1, 4)], lwd = linewidth)
        text(metrics[1] - 5, min(trs + 0.15, 1)*ampl + mn, "SOS", col = colors[1], adj = c(1, 0))
        text(metrics[2] + 5, min(trs + 0.15, 1)*ampl + mn, "EOS", col = colors[4], adj = c(0, 0))
    }
    return(metrics)
    ### The function returns a vector with SOS, EOS, LOS, POP, MGS, rsp, rau, PEAK, MSP and MAU. }
}


#' PhenoDeriv
#' 
#' @inheritParams PhenoTrs
#' 
#' @rdname PhenoExtractMeth
#' @export
PhenoDeriv <- function(fit, IsPlot = TRUE, smspline = TRUE, show.lgd = T, ...){
    PhenoNames <- c("SOS", "POP", "EOS")
    metrics <- setNames(rep(NA, 3), c("sos", "pop", "eos")) # template

    t      <- fit$tout
    values <- last(fit$fits)
    n      <- length(t)

    # get peak of season position
    half.season <- median(which.max(values)) # deal with multiple pop values
    pop <- t[half.season]

    if (all(is.na(values))) return(metrics)
    if (half.season < 5 || half.season > (n - 5)) return(metrics)

    der1   <- D1.phenofit(fit, smspline = smspline)
    # get SOS and EOS according to first order derivative
    # fixed 20180510, ±5 to make sure sos and eos are not in the POP.
    # I_sos <- median(which.max(der1[1:(half.season - 5)]))
    # I_eos <- median(which.min(der1[(half.season+5):length(der1)])) + half.season

    # der.sos (eos) is impossible to occur near the pop.
    nmin  <- 5
    I_sos <- findpeaks( der1[1:(half.season - 5)]         , nups=nmin, ndowns=nmin, npeaks=1, sortstr=TRUE)$X$pos
    I_eos <- findpeaks(-der1[(half.season+5):length(der1)], nups=nmin, ndowns=nmin, npeaks=1, sortstr=TRUE)$X$pos + half.season + 4

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
        op <- par(new = T)
        plot(t, der1, type= "l", lty = 2, lwd = linewidth,
             col = "black", axes = FALSE)
    }
    return(metrics)
}


#' PhenoGu
#' 
#' @inheritParams PhenoTrs
#' @importFrom dplyr last
#' @rdname PhenoExtractMeth
#' @export
PhenoGu <- function(fit, IsPlot = TRUE, smspline = TRUE, ...) {
    PhenoNames <- c("UD", "SD", "DD", "RD")
    metrics <- setNames(rep(NA, 4), c("UD", "SD", "DD", "RD"))

    t      <- fit$tout
    values <- last(fit$fits)
    n      <- length(t)

    # get peak of season position
    half.season <- median(which.max(values)) # deal with multiple pop values
    pop <- t[half.season]

    if (all(is.na(values))) return(metrics)
    if (half.season < 5 || half.season > (n - 5)) return(metrics)

    der1   <- D1.phenofit(fit, smspline = smspline)
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

#' PhenoKl
#' 
#' @inheritParams PhenoTrs
#' @rdname PhenoExtractMeth
#' @export
PhenoKl <- function(fit, IsPlot = TRUE, show.lgd = T, ...) {
    PhenoNames <- c("Greenup", "Maturity", "Senescence", "Dormancy")
    metrics <- setNames(rep(NA, 4), PhenoNames)

    t      <- fit$tout
    values <- last(fit$fits)
    n      <- length(t)
    xlim   <- range(t)

    # get peak of season position
    half.season <- median(which.max(values)) # + 20, half season + 20 was unreasonable
    pop <- t[half.season]

    if (all(is.na(values))) return(metrics)
    if (half.season < 5 || half.season > (n - 5)) return(metrics)

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
        op <- par(new = T)
        plot(t, der.k, xlim = xlim, type= "l",
             lty = 3, lwd = linewidth, col = "black", axes = T) # cex = 1, pch = 20,
    }
    return(setNames(metrics, PhenoNames))
}
