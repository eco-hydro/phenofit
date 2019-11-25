context("test-rcpp_wSG")

test_that("rcpp_wSG works", {
    ## 1. Seq value
    x <- 1:100
    x[50] <- NA
    x[80] <- Inf
    halfwin <- 10
    w <- is.finite(x)

    s1 <- rcpp_wSG(x, halfwin, 1)
    s2 <- rcpp_wSG(x, halfwin, 1, w)

    s2_1 <- rcpp_wSG(x, halfwin, 2)
    s2_2 <- rcpp_wSG(x, halfwin, 2, w)

    expect_equal(s1, 1:100)
    expect_equal(s2, 1:100)
    expect_equal(s2_1, s2_2)

    ## 2. random value
    x2 <- rnorm(100)
    x2[50] <- NA
    x2[80] <- Inf
    w <- is.finite(x)
    s3 <- rcpp_wSG(x2, halfwin, 2)
    s4 <- rcpp_wSG(x2, halfwin, 2, w)
    expect_equal(s3, s4)
})

test_that("rcpp_wSG works", {
    x <- 1:100
    s1 <- rcpp_SG(x, halfwin = 3, d = 2)
    expect_equal(x, s1)
})
