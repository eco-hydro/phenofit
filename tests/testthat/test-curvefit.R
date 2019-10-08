context("test-curvefit")

# simulate vegetation time-series
fFUN = doubleLog.Beck
par  = c(
    mn  = 0.1,
    mx  = 0.7,
    sos = 50,
    rsp = 0.1,
    eos = 250,
    rau = 0.1)
t    <- seq(1, 365, 8)
tout <- seq(1, 365, 1)
y <- fFUN(par, t)

methods <- c("AG", "Beck", "Elmore", "Gu", "Zhang", "Klos")
suppressWarnings(fit <- curvefit(y, t, tout, methods))

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

# Klos is not convergent.
test_that("plot.fFITs works", {
    expect_silent(plot.fFITs(fit))
})

test_that("get_GOF works", {
    expect_silent({info <- get_GOF.fFITs(fit)})
    print(info)
})
