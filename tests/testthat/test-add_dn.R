context("test-add_dn")

test_that("add_dn works", {

    date = seq.Date(as.Date("2010-01-01"), as.Date("2010-12-31"), by = "day")
    d <- data.frame(date)
    days <- c(8, 16)
    dnew <- add_dn(d, days)

    # check colnames of \code{days}
    expect_true(all(is.element(paste0("d", days), colnames(dnew))))
})
