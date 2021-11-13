context("test-qcFUN")


test_that("qcFUN works", {
    set.seed(100)
    QA <- as.integer(runif(100, 0, 2^7))

    r1 <- qc_summary(QA, wmin = 0.2, wmid = 0.5, wmax = 1)
    r2 <- qc_StateQA(QA, wmin = 0.2, wmid = 0.5, wmax = 1)

    r3 <- qc_5l(QA, wmin = 0.2, wmid = 0.5, wmax = 1)
    r4 <- qc_NDVIv4(QA, wmin = 0.2, wmid = 0.5, wmax = 1)

    r_lai <- qc_FparLai(QA, wmin = 0.2, wmid = 0.5, wmax = 1)
    r_spot <- qc_SPOT(QA, wmin = 0.2, wmid = 0.5, wmax = 1)
    r_s2   <- qc_sentinel2(QA, wmin = 0.2, wmid = 0.5, wmax = 1)

    expect_true(is.list(r1) && is.list(r2))
    expect_equal(names(r_lai), names(r_spot))
    expect_true(is.numeric(r3) && is.numeric(r4))
})
