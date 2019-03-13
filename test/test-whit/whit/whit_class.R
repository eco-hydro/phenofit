whitV <- function(d, nptperyear, IsPlot = TRUE){
    lambdas <- 10^(seq(0.1, 6, 0.1))
    INPUT <- check_input(d$date, d$y, d$w, trim = T, maxgap = nptperyear / 4)
    fits  <- whitsmw2(INPUT$y, INPUT$w, INPUT$ylu, nptperyear, iters=1, lambda=lambdas, validation = T)

    pen <- map_dbl(fits, "pen")
    fit <- map_dbl(fits, "fit")

    dfit = diff(fit);
    dpen = diff(pen);
    v    = sqrt(dfit ^ 2 + dpen ^ 2);

    nla = length(lambdas)
    sel = 2:nla
    # midla = sqrt(lambdas[sel] * lambdas[sel - 1])
    midla = 10^( (log10(lambdas[sel]) + log10(lambdas[sel-1]))/2 )
    k <- which.min(v)
    lambda = midla[k]

    if (IsPlot){
        # plot
        par(mfrow=c(2,2))
        plot(fit, pen, xlab = 'log10(Fits)', ylab = 'log10(Penalties)', main = 'L-curve')
        points(fit[k], pen[k], col = "red"); grid()

        plot(log10(midla), v, xlab = 'log10(lambda)', ylab = 'Distance', main = 'V-curve')
        points(log10(lambda), v[k], col = "red"); grid()

        # map_dbl(fits, "cv")  %>% plot(log10(lambdas), ., main = "cv"); grid()
        # map_dbl(fits, "df")  %>% plot(log10(lambdas), ., main = "freedom of parameters"); grid()
        map_dbl(fits, "gcv") %>% plot(log10(lambdas), ., main = "gcv"); grid()
    }
    # yfit <- whitsmw2(INPUT$y, INPUT$w, iters=3, lambdas=lambda, validation = T)
    try(brks  <- season(INPUT, lambda, nptperyear, iters = 3, wFUN = wTSM, IsPlot = TRUE))
    title(x$site[1])
    # fitopt = yfit$fit
    # penopt = yfit$pen
    return(c(lambda = lambda, ny = length(INPUT$y), nw = sum(INPUT$w)))
}

whitsmw <- function (y, w, nptperyear, ylu, wFUN = wTSM, iters = 1, lambda = 100,
    ..., d = 2, validation = FALSE) {
    if (all(is.na(y))) return(y)
    # whittaker smooth
    n  <- length(y)
    E  <- Diagonal(x = rep(1, n))
    D  <- Matrix::diff(E, lag = 1, differences = d)
    D2 <- Matrix::t(D) %*% D

    fits <- list()
    for (i in 1:iters){
        W <- Diagonal(x = w)
        C <- chol(W + (lambda * D2))
        z <- solve(C, solve(t(C), Matrix(w*y)))[, 1]

        if (!validation){
            w <- wFUN(y, z, w, iters, nptperyear, ...)
            z <- check_ylu(z, ylu)
        }
        fits[[i]] <- z
    }

    fits %<>% set_names(paste0('iter', 1:iters))
    if (validation){
        H <- solve(C, solve(t(C), W))
        h <- diag(H)

        n   <- sum(w)
        df  <- sum(h)
        r   <- (y - z)/(1 - h)
        cv  <- sqrt( sum( r^2*w ) /n )
        gcv <- sqrt( sum( (r/(n-df))^2*w ))
        LV  <- whit_V(y, z, w, D) #L curve
        return(c(list(data = as_tibble(c(list(w = w), fits)),
            df = df, cv = cv, gcv = gcv), LV))
    }else{
        return(c(list(data = as_tibble(c(list(w = w), fits))), LV))
    }
}

whit_parts <- function(y, z, w){
    pen <- log( sum( (y - z)^2 * w ))
    # fit <- log10( sum( (D * z)^2 ))
    fit <- log( sum( diff(z, diff = 2)^2 ))
    return(list(pen = pen, fit = fit))
}

whit_optim <- function(y, W, nyear){
    ## select the best lambda
    # nyear <- length(y)/nptperyear
    # lambdas <- 10^seq(1, 6, 0.1)
    error_whitdf <- function(lambda){
        abs(whit_df(lambda) - df)
    }
    lambda <- bisection(.error, 10, 10^6, ln = T, toly = 1e-1)
}

#return: df, cv, gcv
# .error <- function(y, w, nyear, df.year = 8, lambda){
#     f_whitdf <- function(lambda){
#         fit1 <- whitsmw2(y, w, iters=1, lambda=lambda, validation = T)
#         fit2 <- whitsmw (y, w, iters=1, lambda=lambda, validation = T)

#         print(fit1$df)
#         print(fit2$df)
#         df.year*nyear - fit2$df#return
#     }
#     return(f_whitdf)
# }
