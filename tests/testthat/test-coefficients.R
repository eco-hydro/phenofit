context("test-coefficients")

test_that("R2_sign works", {
    expect_silent(R2_critical <- R2_sign(30, NumberOfPredictor = 2, alpha = 0.05))
    expect_equal(names(R2_critical), c("F", "R2"))

    expect_silent(R2_critical <- R2_sign(30, NumberOfPredictor = 2, alpha = c(.01, .05, .1)))
})

test_that("CV, skewness and kurtosis works", {
    set.seed(100)
    x = rnorm(100)

    expect_silent({
        coefs <- cv_coef(x)

        coef_kurtosis <- kurtosis(x)
        coef_skewness <- skewness(x)

        coef_kurtosis <- kurtosis(x, type = 1)
        coef_skewness <- skewness(x, type = 1)

        coef_kurtosis <- kurtosis(x, type = 2)
        coef_skewness <- skewness(x, type = 2)
    })

    # Test for NA values
    x <- c(x, NA)
    coef_kurtosis <- kurtosis(x)
    coef_skewness <- skewness(x)
    expect_equal(rep(NA_real_, 2), c(coef_kurtosis, coef_skewness))
})



