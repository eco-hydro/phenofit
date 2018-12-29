################################################################################


Cairo::CairoPDF("Fig3_GPP_whittaker_v2_wTSM_1e4_1mon_3.pdf", width = 11, height = 6)
op <- par(mfrow = c(2, 3),
          oma = c(1, 2, 2, 1), mar = c(3, 2, 1, 1)) #, yaxt = "n", xaxt = "n"

yfits <- dlply(df, .(site), function(x) {
    tryCatch({
        whitV(x)
        # sitename <- sites[i]
        # file <- sprintf('[%03d]%s.pdf', i, sitename)
        # fprintf("%s\n", file)
        # x <- df[site == sitename]
        # dev.off()
    }, error = function(e) {
        # print(x)
        message(sprintf("[%s]:%s", x$site[1], e$message))
    })
}, .progress = "text")
dev.off()

plot(INPUT$y, type = "l")
lines(yfit$data$iter1, col = "red")
lambda <- 10^(c(1, 2, 3, 4, 5, 6))
I <- match(lambda, lambdas)


# yfit <- whitsmw2(INPUT$y, INPUT$w, INPUT$ylu, nptperyear, iters=3, lambdas=5e4, validation = F)
# plot_input(INPUT, nptperyear)
# colors <- c("red", "blue", "green")
# for (i in 1:(ncol(yfit$data) - 1)){
#     lines(INPUT$t, yfit$data[, i+1, drop = T], col = colors[i], lwd = 2)
# }
try(brks  <- season(INPUT, lambda=1e4, nptperyear, iters = 3, wFUN = wTSM, IsPlot = TRUE))

# f <- .error(INPUT$y, INPUT$w, 11, df.year = 8, lambda)
# bisection(f, 10, 1e6, toly =1e-1)

# f <- .error(INPUT$y, rep(1, length(INPUT$y)), 11, df.year = 8, lambda)
# bisection(f, 10, 1e6, toly =1e-1)

df[site %in% c("AR-Vir", "AU-DaP", "AU-Rig", "CA-NS1", "CH-Cha"), .(site, date, GPP_NT, lat, IGBP)] %>%
    split(.$site) %>%
    writelist_ToXlsx("fluxsites_GPP.xlsx")
