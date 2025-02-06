# context("test-curvefit")

# simulate vegetation time-series
t <- seq(1, 365, 8)
par <- c(mn = 0.1, mx = 0.7, sos = 50, rsp = 0.1, eos = 250, rau = 0.1)
y <- doubleLog.Beck(par, t)
tout <- seq(1, 365, 1)

methods <- c("AG", "Beck", "Elmore", "Gu", "Zhang", "Klos")[1:5]

# weird: cpp is much slower
system.time(suppressWarnings(fit <- curvefit(y, t, tout = tout, methods, use.cpp = FALSE)))
system.time(suppressWarnings(fit_cpp <- curvefit(y, t, tout = tout, methods, use.cpp = TRUE)))

test_that("curvefit works", {
    p1 = get_param(fit_cpp)[1:5]
    p2 = get_param(fit)[methods][1:5]
    expect_true(all.equal(p1, p2, tolerance = 1e-5))
    expect_silent(r <- get_param(list(fit)))

    # For Klos, the result of C++ is slightly different from that of R version.
    diff = get_fitting(fit_cpp)$ziter2 - get_fitting(fit)$ziter2
    expect_true(max(abs(diff)) <= 1e-5)
    expect_silent(dfit <- get_param(list(fit)))
})

# Klos is not convergent.
test_that("plot.fFITs works", {
    expect_silent(plot.fFITs(fit))
})

test_that("get_GOF works", {
    expect_silent({info <- get_GOF.fFITs(fit)})
    print(info)
})

# rbenchmark::benchmark(
#     fit1 <- curvefit(y, t, tout, methods[1:5]),
#     fit2 <- curvefit(y, t, tout, methods[1:5], use.cpp = TRUE),
#     replications = 10
# )

# test_that("curvefit works", {
#     expect_silent(suppressWarnings(
#         fit <- curvefit(y, t, tout, methods)
#     ))
# })

# test_that("print.curvefit works", {
#     expect_silent(print(fit))
# })
