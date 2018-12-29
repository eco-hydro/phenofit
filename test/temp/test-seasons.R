{
    # test break points
    x <- df[site == "AU-Rig"]

    plot(GPPobs~date, x, type = "b"); grid()
    lm_fit <- lm(GPPobs~date, x)

    dtest <- davies.test(lm_fit, ~date, k = 10)
    sg <- segmented(lm_fit, ~date, psi = list(date=NA),
                    control=seg.control(stop.if.error=FALSE, n.boot = 0, it.max = 1000))
    abline(v = sg$psi[, 2], col = "red", lty = 2)


    piece_fit <- segmented(lm_fit, seg.Z = ~x, psi = psi_init,
                           control = seg.control(stop.if.error = FALSE, n.boot = 0, it.max = 100))

    bp <- breakpoints(x$GPPobs~x$date, breaks = 9)
    plot(bp)
    abline(v = x$date[bp$breakpoints], col = "red", lty = 2)

    nptperyear <- 46
    y <- na.approx(x$GPPobs)

    plot(y, type = "b")
    # points(x$GPPobs, col = "blue")

    bcp.1a <- bcp(x$GPPobs)
    plot(bcp.1a, main="Univariate Change Point Example")
    legacyplot(bcp.1a)
}

sitename <- sites[i]
file <- sprintf('[%03d]%s.pdf', i, sitename)
fprintf("%s\n", file)
x <- df[site == sitename]

x <- x[, .(GPP = median(GPP, na.rm = T), date = first(date)), .(site,  year, d8)]

nptperyear = 46

INPUT <- check_input(x$date, x$GPP, trim = T, maxgap = nptperyear/4)
brks  <- season(INPUT, lambda=15, nptperyear, iters = 3, wFUN = wTSM, IsPlot = TRUE)


yfits <- whitsmw(INPUT$y, INPUT$w, nptperyear, INPUT$ylu, wFUN, iters = 3, lambda = 2000)

plot_input(INPUT, nptperyear)
colors <- c("red", "blue", "green")
for (i in 1:(ncol(yfits) - 1)){
    lines(INPUT$t, yfits[, i+1, drop = T], col = colors[i])
}

plot(GPP~date, a, type = "b")
