context("season3y")
source('helper_MOD13A1.R')

# lambda   <- init_lambda(INPUT$y) # lambda for whittaker
# param = listk(
#     INPUT, nptperyear,
#     FUN = whitsmw2, wFUN = wBisquare, iters = 2,
#     lambda,
#     IsPlot = IsPlot, plotdat = d,
#     south = sp$lat[1] < 0,
#     rymin_less = 0.6, ypeak_min = ypeak_min,
#     max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5
# )

param = listk(
    INPUT, south = sp$lat[1] < 0, 
    rFUN = wWHIT, wFUN = wTSM, iters = 3,
    r_min = 0,
    IsPlot = IsPlot, plotdat = d, print = FALSE, IsOnlyPlotbad = F
)

param$rFUN <- wSG
# brks <- do.call(season_3y, param)

test_that("`season_3y` with wWHIT", {
    expect_silent(brks <- do.call(season_3y, param))
})

test_that("`season_3y` with wHANTS", {
    param$rFUN <- wHANTS
    expect_silent(brks <- do.call(season_3y, param))
})

test_that("`season_3y` with wSG", {
    param$rFUN <- wSG
    expect_silent(brks <- do.call(season_3y, param))
})
