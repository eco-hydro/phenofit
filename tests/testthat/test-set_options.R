test_that("set_options works", {
    set_options(verbose_season_mov = T, season = list(lambda = 2))
    expect_equal(get_options("season")$lambda, 2)
})
