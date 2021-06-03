context("test-v_curve")

data("MOD13A1")

dt <- tidy_MOD13(MOD13A1$dt)
st <- MOD13A1$st

sitename <- dt$site[1]
d     <- dt[site == sitename, ] # get the first site data
sp    <- st[site == sitename, ] # station point
# global parameter
IsPlot = TRUE
nptperyear = 23

dnew     <- add_HeadTail(d, nptperyear = nptperyear) # add one year in head and tail
INPUT    <- check_input(dnew$t, dnew$y, dnew$w, nptperyear,
                        maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
# INPUT$y0 <- dnew$y   # raw time-series, for visualization

test_that("v_curve works", {
    expect_silent({
        lg_lambdas <- seq(0, 3, 0.1)
        r <- v_curve(INPUT, lg_lambdas, d = 2, IsPlot = TRUE)
    })
})
