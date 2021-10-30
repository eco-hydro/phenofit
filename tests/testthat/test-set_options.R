test_that("set_options works", {
    wmin = 0.12
    set_options(season = list(wmin = 0.12))
    expect_equal(get_options("season")$wmin, wmin)
})
