context("test-qcFUN")

set.seed(100)
QA <- as.integer(runif(100, 0, 2^7))

r1 <- qc_summary(QA, wmin = 0.2, wmid = 0.5, wmax = 1)
r2 <- qc_StateQA(QA, wmin = 0.2, wmid = 0.5, wmax = 1)
r3 <- qc_5l(QA, wmin = 0.2, wmid = 0.5, wmax = 1)
r4 <- qc_NDVIv4(QA, wmin = 0.2, wmid = 0.5, wmax = 1)


test_that("qcFUN works", {
    expect_true(is.list(r1) && is.list(r2))
    expect_true(is.numeric(r3) && is.numeric(r4))
})
