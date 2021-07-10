data("CA_NS6")
d = CA_NS6

nptperyear <- 23
INPUT <- check_input(d$t, d$y, d$w, QC_flag = d$QC_flag,
     nptperyear = nptperyear, south = FALSE,
     maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
# plot_input(INPUT)

# Rough fitting and growing season dividing
wFUN <- "wTSM"
brks2 <- season_mov(INPUT,
    options = list(
        rFUN = smooth_wWHIT, wFUN = wFUN,
        r_min = 0.05, ypeak_min = 0.05,
        lambda = 10,
        verbose = FALSE
    ))
# plot_season(INPUT, brks2, d)
# Fine fitting
fit <- curvefits(
    INPUT, brks2,
    options = list(
        methods = c("AG", "Beck", "Elmore", "Zhang"), #,"klos", "Gu"
        wFUN = wFUN,
        nextend = 2, maxExtendMonth = 2, minExtendMonth = 1, minPercValid = 0.2
    )
)

r_param = get_param(fit)
r_pheno = get_pheno(fit)
r_gof = get_GOF(fit)
d_fit = get_fitting(fit)

g <- plot_curvefits(d_fit, brks2)
grid::grid.newpage(); grid::grid.draw(g)
