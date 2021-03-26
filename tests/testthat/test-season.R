context("season")

source('helper_MOD13A1.R')
# source('tests/testthat/helper_MOD13A1.R')
lambda   <- init_lambda(INPUT$y) # lambda for whittaker

# `season` produce a pdf file here
param = listk(
    INPUT,
    rFUN = smooth_wWHIT, wFUN = wBisquare, iters = 2,
    lambda,
    IsPlot = IsPlot, plotdat = d,
    rymin_less = 0.6, ypeak_min = ypeak_min,
    max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5
)


test_that("`season` with smooth_wWHIT and calendarYear", {
    param$calendarYear <- TRUE
    expect_silent(brks <- do.call(season, param))
})

test_that("`season` with smooth_wWHIT", {
    param$calendarYear <- FALSE
    expect_silent(brks <- do.call(season, param))
})

test_that("`season` with wHANTS", {
    param$rFUN <- smooth_wHANTS
    expect_silent(brks <- do.call(season, param))
})

test_that("`season` with smooth_wSG", {
    param$rFUN <- smooth_wSG
    expect_silent(brks <- do.call(season, param))
})
