test_that("fFITs functions work", {
    fFUN = doubleLog.Beck
    par  = c( mn  = 0.1 , mx  = 0.7 , sos = 50 , rsp = 0.1 , eos = 250, rau = 0.1)

    t    <- seq(1, 365, 8)
    tout <- seq(1, 365, 1)
    y <- fFUN(par, t)

    methods <- c("AG", "Beck", "Elmore", "Gu", "Zhang") # "Klos" too slow
    fFITs <- curvefit(y, t, tout, methods)

    print(fFITs$model[1])
    expect_silent(plot(fFITs))
    # expect_equal(2 * 2, 4)
})
