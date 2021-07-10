data("CA_NS6")
d = CA_NS6

nptperyear <- 23
INPUT <- check_input(d$t, d$y, d$w, QC_flag = d$QC_flag,
     nptperyear = nptperyear, south = FALSE,
     maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
# plot_input(INPUT)

wFUN = "wTSM"
# all year as a whole
brks  <- season(INPUT,
    rFUN = smooth_wWHIT, wFUN = wFUN, lambda = 10)
plot_season(INPUT, brks, d)

# curve fitting by year
brks2 <- season_mov(INPUT,
    options = list(
        rFUN = "smooth_wWHIT", wFUN = wFUN,
        lambda = 10,
        r_min = 0.05, ypeak_min = 0.05,
        verbose = TRUE
    ))
plot_season(INPUT, brks2)
