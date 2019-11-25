# Rough fitting and growing season dividing
brks2 <- season_mov(INPUT,
    rFUN = smooth_wWHIT, wFUN = wFUN,
    plotdat = d, IsPlot = IsPlot, print = FALSE, IsPlot.OnlyBad = FALSE)
# Fine fitting
fit <- curvefits(
    INPUT, brks2,
    methods = c("AG", "Beck", "Elmore", "Zhang"), #,"klos", "Gu"
    wFUN = wFUN,
    nextend = 2, maxExtendMonth = 2, minExtendMonth = 1, minPercValid = 0.2,
    print = TRUE, verbose = FALSE)
## visualization
df_fit <- get_fitting(fit)
g <- plot_phenofit(df_fit, brks2)
grid::grid.newpage(); grid::grid.draw(g)
