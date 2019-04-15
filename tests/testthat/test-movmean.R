context("test-movmean")

test_that("movmean works", {
    x <- 1:100
    x[50] <- NA
    x[80] <- Inf
    halfwin <- 10
    s1 <- movmean(x, 2, SG_style = FALSE)

    expect_true(all(is.finite(s1)))

    x <- c(rep(Inf, 4), 1, 2, 3, NA)
    expect_silent(s2 <- movmean(x, 2, SG_style = TRUE))
    expect_equal(movmean(x, 2, SG_style = FALSE)[1:3], c(NA, NA, 1))
})
