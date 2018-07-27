# Functions for working with the V-curve

v_point = function(y, w = 0 * y + 1, lambda = 100, d = 2) {
    # Compute the value of the normalized V-curve for one value of lambda
    # Prepare for smoothing
    n = length(y)
    E = diag.spam(n)
    D = diff(E, diff = d)
    P = t(D) %*% D

    # Smooth for  log-lambdas to the left and to the right
    z     = whit2(y, lambda, w)
    pz    = P %*% z #D' * D * z
    zgrad = lambda * log(10) * whit2(- pz/ w, lambda, w) #whit2(- pz * lambda, lambda, 1)
    # zgrad1 = whit2(-lambda * pz, lambda, w)

    fit   = sum(w * (y - z) ^ 2)
    dlfit = 2 * sum(-zgrad * w * (y - z)) / fit
    pen   = sum(z * pz)
    dlpen = 2 * sum(pz * zgrad) / pen

    # Take distance
    v = sqrt(dlfit ^ 2 + dlpen ^ 2)
    return(v)
}

# sometimes not converge
v_opt = function(y, w = 0 * y + 1, d = 2, lambdas = c(0, 4), tol = 0.01) {
    # Locate the optimal value of log10(lambda) with optimizer
    # Specify bounds of search range for log10(lambda) in paramter 'lambdas'

    v_fun = function(lla, y, w, d) v_point(y, w, 10 ^ lla, d)
    op = optimize(v_fun, lambdas, y, w, d, tol = tol)
    return(op$minimum)
}

#' v_curve
#'
#' V-curve is used to optimize Whittaker parameter lambda.
#' Update 20180605 add weights updating to whittaker lambda selecting
#'
#' @inheritParams whitsmw2
#' @param d difference order
#' @param IsPlot Boolean. Whether to plot figure?
#'
#' @export
v_curve = function(INPUT, nptperyear, lambdas,  d = 2, IsPlot = F,
    wFUN = wTSM, iters=2) {
    # Compute the V-cure
    y <- INPUT$y
    w <- INPUT$w

    if (length(unique(y)) == 0) return(NULL)

    param <- c(INPUT, nptperyear = nptperyear, wFUN = wFUN, iters=iters,
        second = FALSE, lambdas=NA)

    fits = pens = NULL
    for (lla in lambdas) {
        # param$lambdas <- 10^lla
        # z    <- do.call(whitsmw2, param)$zs %>% dplyr::last()

        z    = whit2(y, 10 ^ lla, w)
        fit  = log(sum(w * (y - z) ^ 2))
        pen  = log(sum(diff(z, diff = d) ^2))
        fits = c(fits, fit)
        pens = c(pens, pen)
    }

    # Construct V-curve
    dfits   = diff(fits)
    dpens   = diff(pens)
    llastep = lambdas[2] - lambdas[1]
    v       = sqrt(dfits ^ 2 + dpens ^ 2) / (log(10) * llastep)

    nla     = length(lambdas)
    lamids  = (lambdas[-1] + lambdas[-nla]) / 2
    k       = which.min(v)
    lambda  = 10 ^ lamids[k]

    # param$lambdas <- lambda
    # fit <- do.call(whitsmw2, param)
    # d_sm <- fit %$% c(ws, zs) %>% as.data.table() %>% cbind(t = INPUT$t, .)

    z    <- whit2(y, lambda, w)
    d_sm <- data.table(z, t = INPUT$t)

    if (IsPlot) {
        ylim = c(0, max(v))
        plot(lamids, v, type = 'l', col = 'blue', ylim = ylim,
           xlab = 'log10(lambda)')
        points(lamids, v, pch = 16, cex = 0.5, col = 'blue' )
        abline(h = 0, lty = 2, col = 'gray')
        abline(v = lamids[k], lty = 2, col = 'gray', lwd = 2)
        title(sprintf("v-curve, lambda = %5.2f", lambda))
        grid()
    }

    return( list(fit = d_sm, lambdas = lamids, v = v, lambda = lambda, vmin = v[k]))
}

# According to v-curve theory, get the optimal lambda value.
# 
# Whittaker balanced the fidelity and smooth. The agreement index maybe poor 
# than others. But it'is much smoothing.
# 
optim_lambda <- function(sitename, df, deltaT, extent = T,
    IsPlot = F, IsSave = F, file = "test_whit_lambda.pdf",
    wFUN = wBisquare, iters = 2){
    # sitename <- sites[i]#; grp = 1
    nperiod <- ceiling(length(2000:2017)/deltaT)

    d     <- df[site == sitename]
    dnew  <- add_HeadTail(d) #
    INPUT <- check_input(dnew$t, dnew$y, dnew$w, maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
    if (length(unique(INPUT$y)) <= 5) return(NULL)

    years <- year(ymd(dnew$t))
    cat(sprintf('site: %s ...\n', sitename))

    # res <- numeric(nperiod)*NA_real_
    res <- list()

    for (i in 1:nperiod){
        year      <- (i - 1)*deltaT + 2000
        year_beg  <- year
        year_end  <- min(year + deltaT - 1, 2017)

        year_beg_ext <- ifelse(extent, year_beg-1, year_beg)
        year_end_ext <- ifelse(extent, year_end+1, year_end)

        I     <- which(years >= year_beg & years <= year_end)
        I_ext <- which(years >= year_beg_ext & years <= year_end_ext)

        INPUT_i <- lapply(INPUT[1:3], `[`, I_ext) %>% c(INPUT[5])

        res[[i]] <- tryCatch({
            if (IsPlot) par(mfrow = c(2, 1), mar = c(2.5, 2.5, 1, 0.2),
                mgp = c(1.3, 0.6, 0), oma = c(0, 0, 0.5, 0))

            vc <- v_curve(INPUT_i, nptperyear, lambdas = seq(-1, 3, by = 0.01), d = 2,
                          wFUN = wFUN, iters = iters,
                IsPlot = IsPlot)

            ind <- match(I, I_ext)
            vc$fit <- vc$fit[ind, ]

            # coefficient to construct Whittaker lambda formula
            vc$coef <- dnew[I_ext, .(mean = mean(y),
                            sd = sd(y),
                            kurtosis = kurtosis(y, type = 2),
                            skewness = skewness(y, type = 2))] %>% as.list()

            if (IsPlot){
                plotdata(INPUT_i, nptperyear, wmin = 0.1)
                colors <- c("blue", "red")
                lines(vc$fit$t, last(vc$fit), col = "blue", lwd = 1.2)
                # lines(vc$fit$t, vc$fit$ziter2, col = "red" , lwd = 1.2)
            }
            vc
            # listk(lambda = vc$lambda) #, vc
        }, error = function(e){
            message(sprintf("[e] %s, %d: %s", as.character(sitename), i, e$message))
            #return(NA)
        })
    }

    tryCatch({
        ## visualization
        res %<>% rm_empty()
        lambda <- map_dbl(res, "lambda")
        df_sm  <- res %>% map_df("fit") %>% data.table()
        coefs  <- res %>% map_df("coef")

        # browser(), iter2
        info <- merge(df_sm, d, by = "t") %$% GOF(y, z) %>% as.list()

        if (IsSave){
            titlestr <- info[c(1:3, 5, 7)] %>%
                {sprintf("%s = %.2f", names(.), .) %>% paste(collapse = ", ") }
            cairo_pdf(file, 10, 4)
            par(mar = c(2.5, 2.5, 1, 0.2),
                mgp = c(1.3, 0.6, 0), oma = c(0, 0, 0.5, 0))
            plotdata(d, nptperyear)
            lines(ziter1~t, df_sm, col = "blue", lwd = 1.2)
            lines(ziter2~t, df_sm, col = "red", lwd = 1.2)
            title(titlestr)
            dev.off()
            file.show(file)
        }
        listk(data = df_sm, coef = coefs, gof = info, lambda) # return
    }, error = function(e){
        message(sprintf("[e] %s, %d: %s", as.character(sitename), i, e$message))
        #return(NULL)
    })
}

#' Initial lambda value of whittaker taker
#' @param y Numeric vector
#' @export
init_lambda  <- function(y){
    y        <- y[!is.na(y)] #rm NA values
    mean     <- mean(y)
    sd       <- sd(y)
    cv       <- sd/mean
    skewness <- skewness(y, type = 2)
    kurtosis <- kurtosis(y, type = 2)
    # lambda was transformed by log10
    # lambda   <- 0.555484 + 1.671514*mean - 3.434064*sd - 0.052609*skewness + 0.009057*kurtosis
    # lambda   <- 0.555465 + 1.501239*mean - 3.204295*sd - 0.031902*skewness # Just three year result

    lambda <- 0.831120 + 1.599970*mean - 4.094027*sd - 0.035160*cv - 0.063533*skewness # all year
    # lambda <- 0.817783 + 1.803588*mean - 4.263469*sd - 0.038240*cv - 0.066914*skewness - 0.011289*kurtosis  #3y
    return(10^lambda)
}

## All year togather
# Call:
# lm_fun(formula = lambda ~ mean + sd + cv + skewness, data = d,
#     na.action = na.exclude)

# Residuals:
#     Min      1Q  Median      3Q     Max
# -2.4662 -0.4267  0.1394  0.4144  2.5824

# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
# (Intercept)  0.831120   0.021897  37.956  < 2e-16 ***
# mean         1.599970   0.089914  17.794  < 2e-16 ***
# sd          -4.094027   0.168844 -24.247  < 2e-16 ***
# cv          -0.035160   0.008459  -4.157 3.25e-05 ***
# skewness    -0.063533   0.007966  -7.976 1.62e-15 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.5851 on 15135 degrees of freedom
# Multiple R-squared:  0.1572,    Adjusted R-squared:  0.1569
# F-statistic: 705.5 on 4 and 15135 DF,  p-value: < 2.2e-16

#          term    estimate   std.error  statistic       p.value
# 1 (Intercept)  0.83112003 0.021897138  37.955647 3.306152e-301
# 2        mean  1.59997023 0.089914112  17.794428  4.039120e-70
# 3          sd -4.09402663 0.168843860 -24.247412 1.876008e-127
# 4          cv -0.03516045 0.008458787  -4.156677  3.246863e-05
# 5    skewness -0.06353256 0.007965580  -7.975886  1.620571e-15

## Three year fitting result
# Call:
# lm_fun(formula = lambda ~ (mean + sd + cv + skewness + kurtosis),
#     data = d, na.action = na.exclude)

# Residuals:
#     Min      1Q  Median      3Q     Max
# -2.6601 -0.4564  0.0490  0.4418  2.7551

# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
# (Intercept)  0.817783   0.010569  77.379  < 2e-16 ***
# mean         1.803588   0.042016  42.926  < 2e-16 ***
# sd          -4.263469   0.081937 -52.033  < 2e-16 ***
# cv          -0.038240   0.004041  -9.462  < 2e-16 ***
# skewness    -0.066914   0.003762 -17.785  < 2e-16 ***
# kurtosis     0.011289   0.001506   7.496 6.62e-14 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.6639 on 90639 degrees of freedom
# Multiple R-squared:  0.1541,    Adjusted R-squared:  0.1541
# F-statistic:  3303 on 5 and 90639 DF,  p-value: < 2.2e-16

#          term    estimate   std.error  statistic      p.value
# 1 (Intercept)  0.81778299 0.010568590  77.378628 0.000000e+00
# 2        mean  1.80358830 0.042016401  42.925816 0.000000e+00
# 3          sd -4.26346883 0.081937136 -52.033413 0.000000e+00
# 4          cv -0.03823967 0.004041253  -9.462331 3.080333e-21
# 5    skewness -0.06691403 0.003762368 -17.785083 1.216689e-70
# 6    kurtosis  0.01128888 0.001505916   7.496350 6.621350e-14

#' skewness and kurtosis
#'
#' Inherit from package `e1071`
#' @param x a numeric vector containing the values whose skewness is to be
#' computed.
#' @param na.rm a logical value indicating whether NA values should be stripped
#' before the computation proceeds.
#' @param type an integer between 1 and 3 selecting one of the algorithms for
#' computing \code{\link[e1071]{skewness}}.
#'
#' @export
kurtosis <- function (x, na.rm = FALSE, type = 3) {
    if (any(ina <- is.na(x))) {
        if (na.rm) {
            x <- x[!ina]
        } else {
            return(NA_real_)
        }
    }
    if (!(type %in% (1:3))) stop("Invalid 'type' argument.")
    n <- length(x)
    x <- x - mean(x)
    r <- n * sum(x^4)/(sum(x^2)^2)

    y <- if (type == 1) {
        r - 3
    } else if (type == 2) {
        if (n < 4) {
            warning("Need at least 4 complete observations.")
            return(NA_real_)
        }
        ((n + 1) * (r - 3) + 6) * (n - 1)/((n - 2) * (n - 3))
    } else{
        r * (1 - 1/n)^2 - 3
    }
    y
}

#' @rdname kurtosis
#' @export
skewness <- function (x, na.rm = FALSE, type = 3) {
    if (any(ina <- is.na(x))) {
        if (na.rm) {
            x <- x[!ina]
        } else {
            return(NA_real_)
        }
    }
    if (!(type %in% (1:3))) stop("Invalid 'type' argument.")
    n <- length(x)
    x <- x - mean(x)
    y <- sqrt(n) * sum(x^3)/(sum(x^2)^(3/2))

    if (type == 2) {
        if (n < 3){
            warning("Need at least 3 complete observations.")
            return(NA_real_)
        }
        y <- y * sqrt(n * (n - 1))/(n - 2)
    } else if (type == 3){
        y <- y * ((1 - 1/n))^(3/2)
    }
    y
}
