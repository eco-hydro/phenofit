context("test-movmean")

test_that("movmean works", {
    x <- 1:100
    x[50] <- NA
    x[80] <- Inf
    halfwin <- 10
    s1 <- movmean(x, 2)

    expect_true(all(is.finite(s1)))
})
