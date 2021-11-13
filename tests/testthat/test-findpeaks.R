context("findpeaks")

test_that("findpeaks works", {

    x <- seq(0, 1, len = 1024)
    pos <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81)
    hgt <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
    wdt <- c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005)
    pSignal <- numeric(length(x))
    for (i in seq(along=pos)) {
        pSignal <- pSignal + hgt[i]/(1 + abs((x - pos[i])/wdt[i]))^4
    }

    expect_equal( nrow(findpeaks(pSignal, y_max=4, sortstr=TRUE)$X), 5,
                  info = 'r_min works') # 5 points amplitude > 4
    expect_equal( nrow(findpeaks(pSignal, y_max=3, sortstr=TRUE)$X), 8,
                  info = 'r_min works') # 8 points amplitude > 3
    expect_equal( nrow(findpeaks(pSignal, y_min=4.1, sortstr=TRUE)$X), 3,
                  info = 'r_min works')

    ## add test for Henan multiple growing season
    ypred = c(0.329, 0.341, 0.356, 0.476, 0.539, 0.569, 0.645, 0.602, 0.515, 0.329, 
        0.128, 0.596, 0.535, 0.483, 0.355, 0.355, 0.128, 0.148, 0.176, 0.244)
    expect_equal(nrow(findpeaks(-ypred, y_min = 0.0)$X), 2)
})
