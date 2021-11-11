data("CA_NS6")
d <- CA_NS6

nptperyear <- 23
INPUT <- check_input(d$t, d$y, d$w,
    QC_flag = d$QC_flag,
    nptperyear = nptperyear, south = FALSE,
    maxgap = nptperyear / 4, alpha = 0.02, wmin = 0.2
)
# plot_input(INPUT)

wFUN <- "wTSM"
# all year as a whole
options = list(rFUN = "smooth_wWHIT", wFUN = wFUN, lambda = 10)
brks <- season(INPUT, lambda = 10)
plot_season(INPUT, brks, d)

brks2 = season_input(INPUT, options)
all.equal(brks2, brks)

c(d_fit, info_peak) %<-% roughFit(INPUT)
d_season = find_season.peaks(d_fit, info_peak)

c(t, ypred) %<-% d_fit[, .(t, ziter2)]
d_season = find_season.default(ypred, t)
all.equal(brks$dt, d_season)

# opt <- .options$season
# brks$fit - d_fit # function passed test

# curve fitting by year
brks_mov <- season_mov(INPUT,
    options = list(
        rFUN = "smooth_wWHIT", wFUN = wFUN,
        lambda = 10,
        r_min = 0.05, ypeak_min = 0.05,
        verbose = TRUE
    )
)
plot_season(INPUT, brks_mov)
