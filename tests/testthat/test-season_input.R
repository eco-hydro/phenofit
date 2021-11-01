test_that("movmean works", {
    # set global options 
    set_options(
        wFUN = wTSM, wmin = 0.2,
        verbose = FALSE,
        season = list(
            # rFUN = "smooth_wWHIT", lambda = 0.5,
            rFUN = "smooth_wHANTS", nf = 6,
            maxExtendMonth = 12, r_min = 0.05
        )
    )
    brk = season_input(input_single)
    # plot_season(input_single, brk)
    expect_equal(nrow(brk$dt), 2)

    brk = season_input(input_single, options = list(r_min = 0.02))
    expect_equal(nrow(brk$dt), 2)

    brk = season_input(input_single, options = list(r_min = 0.0))
    expect_equal(nrow(brk$dt), 2)
})
