context("ifelse2")

test_that("ifelse2 works", {
    expect_equal(1:4, ifelse2(TRUE, 1:4, 1:10))
    expect_equal(1:10, ifelse2(FALSE, 1:4, 1:10))
})
