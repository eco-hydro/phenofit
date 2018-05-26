CairoPDF(sprintf("optim_lambda, nptperyear=%d_season_bad_3.pdf", nptperyear),
         width = 10, height = 2*6)

par(mfrow = c(6, 1), mar = c(2, 3, 2, 1), mgp = c(1.5, 0.6, 0))

res <- list()
stat <- list()
for (i in seq_along(sites)){
    runningId(i)
    sitename <- sites[i]#; grp = 1
    d    <- df[site == sitename]
    d <- d[1:(23*3), ]
    t    <- d$t
    y    <- d$y
    w    <- d$w
    lat  <- d$lat[1]
    IGBP <- d$IGBP[1]
    str_title <- sprintf("[%02d] %s, IGBP = %s, log10(lambda) = %.3f", i, sitename, IGBP, log10(vc$lambda))

    if ( sum(d$w >= 0.5) > nrow(d) * 0.3){
        INPUT <- check_input(d$t, d$y, d$w, trim = T, maxgap = nptperyear / 4, alpha = 0.02)

        vc <- v_curve(INPUT$y, w = INPUT$w, llas = seq(-2, 2, by = 0.01), d = 2, show = F)

        brks <- season(INPUT, lambda = vc$lambda, nptperyear, iters = 3, wFUN = wFUN, IsPlot = F,
                       south = d$lat[1] < 0,
                       Aymin_less = 0.6, ymax_min = ymax_min,
                       max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5,
                       plotdat = d) #, ...
        if (is.null(brks)) next()

        stat <- GOF(d$y, brks$whit$iter3, INPUT$w)

        if (stat[['NASH']] < 0.4){
            brks <- season(INPUT, lambda = vc$lambda, nptperyear, iters = 3, wFUN = wFUN, IsPlot = T,
                           south = d$lat[1] < 0,
                           Aymin_less = 0.6, ymax_min = ymax_min,
                           max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5) #, ...
            title(str_title)
        }
        # Plot data and smooth
        # plot_qc(t, d$y, d$w, xlab = '', ylab = 'EVI')
        # plot(t, INPUT$y, type = 'l', col = 'darkgrey', xlab = '', ylab = 'LAI');
        # abline(v = t[seq(1, length(INPUT$y), nptperyear)], col = "grey60", lty = 2)
        # lines(t, vc$z, col = 'blue', lwd = 1.5)
        res[[i]]  <- listk(site = sitename, lat, IGBP, lambda = vc$lambda, vc)

    } else{
        message(str_title)
    }

}

dev.off()
res %<>% set_names(sites)
