context("season3y")
# source('helper_MOD13A1.R')

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
    INPUT,
    rFUN = smooth_wWHIT, wFUN = wTSM, iters = 3,
    r_min = 0,
    IsPlot = IsPlot, plotdat = d, IsOnlyPlotbad = F
)

set_options(verbose_season_mov = FALSE)
param$rFUN <- smooth_wSG
# brks <- do.call(season_mov, param)

test_that("`season` with smooth_wWHIT and calendarYear", {
    param$calendarYear <- TRUE
    expect_silent(brks <- do.call(season_mov, param))
})

test_that("`season_mov` with smooth_wWHIT", {
    param$calendarYear <- FALSE
    expect_silent(brks <- do.call(season_mov, param))
})

test_that("`season_mov` with wHANTS", {
    param$rFUN <- smooth_wHANTS
    expect_silent(brks <- do.call(season_mov, param))
})

test_that("`season_mov` with smooth_wSG", {
    param$rFUN <- smooth_wSG
    expect_silent(brks <- do.call(season_mov, param))
})
