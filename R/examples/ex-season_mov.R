data("CA_NS6")
d <- CA_NS6

nptperyear <- 23
INPUT <- check_input(d$t, d$y, d$w,
    QC_flag = d$QC_flag,
    nptperyear = nptperyear, south = FALSE,
    maxgap = nptperyear / 4, alpha = 0.02, wmin = 0.2
)

# curve fitting by year
brks_mov <- season_mov(INPUT,
    options = list(
        rFUN = "smooth_wWHIT", wFUN = "wTSM",
        lambda = 10,
        r_min = 0.05, ypeak_min = 0.05,
        verbose = TRUE
    )
)
plot_season(INPUT, brks_mov)

rfit <- brks2rfit(brks_mov)
# Phenological Metrics from rough fitting
r <- get_pheno(rfit)
