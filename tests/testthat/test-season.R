context("season")

source('helper_MOD13A1.R')
# source('tests/testthat/helper_MOD13A1.R')
lambda   <- init_lambda(INPUT$y) # lambda for whittaker

param = listk(
    INPUT,
    FUN = wWHIT, wFUN = wBisquare, iters = 2,
    lambda,
    IsPlot = IsPlot, plotdat = d,
    south = sp$lat[1] < 0,
    rymin_less = 0.6, ypeak_min = ypeak_min,
    max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5
)

test_that("`season` with wWHIT", {
    expect_silent(brks <- do.call(season, param))
})

test_that("`season` with wHANTS", {
    param$FUN <- wHANTS
    expect_silent(brks <- do.call(season, param))
})

test_that("`season` with wSG", {
    param$FUN <- wSG
    expect_silent(brks <- do.call(season, param))
})
