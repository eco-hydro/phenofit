# all year as a whole
brks  <- season(INPUT,
    rFUN = wWHIT, wFUN = wFUN,
    lambda = 10,
    plotdat = d, IsPlot = IsPlot, print = FALSE, IsPlot.OnlyBad = FALSE)
# curve fitting by year
brks2 <- season_3y(INPUT,
    rFUN = wWHIT, wFUN = wFUN,
    lambda = 10,
    plotdat = d, IsPlot = IsPlot, print = FALSE, IsPlot.OnlyBad = FALSE)
