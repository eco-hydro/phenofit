# library(magrittr)
# library(purrr)

## 1. Test for `optim_pheno`
# simulate vegetation time-series
FUN = doubleLog_Beck
par  = c( mn  = 0.1 , mx  = 0.7 , sos = 50 , rsp = 0.1 , eos = 250, rau = 0.1)
par0 = c( mn  = 0.15, mx  = 0.65, sos = 100, rsp = 0.12, eos = 200, rau = 0.12)

t    <- seq(1, 365, 8)
tout <- seq(1, 365, 1)
y    <- FUN(par, t)

methods = c("BFGS", "ucminf", "nlm", "nlminb")
opt1 <- I_optim(par0, doubleLog_Beck, y, t, methods) # "BFGS", "ucminf", "nlm",
# opt2 <- I_optimx(prior, fFUN, y, t, tout, )

sFUN   = "doubleLog.Beck" # doubleLog.Beck
r <- optim_pheno(par0, sFUN, y, t, tout, method = methods[4],
                 nptperyear = 46, iters = 2, wFUN = wTSM, verbose = FALSE, use.julia = FALSE)

microbenchmark::microbenchmark(
    r <- optim_pheno(par0, sFUN, y, t, tout, method = methods[4],
                     nptperyear = 46, iters = 2, wFUN = wTSM, verbose = FALSE, use.julia = TRUE),
    r <- optim_pheno(par0, sFUN, y, t, tout, method = methods[4],
                     nptperyear = 46, iters = 2, wFUN = wTSM, verbose = FALSE, use.julia = FALSE)
)

profvis::profvis(
    foreach(i = 1:100) %do% {
        r <- optim_pheno(par0, sFUN, y, t, tout, method = methods[4],
                         nptperyear = 46, iters = 2, wFUN = wTSM, verbose = FALSE, use.julia = TRUE)
    }
)

## 2. Test for `curvefits`

source('tests/testthat/helper_MOD13A1.R')
wFUN = wTSM # wBisquare #

# The `maxExtendMonth` in season_mov and curvefits is different
# lambda   <- init_lambda(INPUT$y) # lambda for whittaker
brks2 <- season_mov(INPUT,
    rFUN = "smooth_wWHIT", wFUN = wFUN,
    plotdat = d, IsPlot = IsPlot, IsPlot.OnlyBad = F, print = F)

param <- list(
    INPUT, brks2,
    methods = c("AG", "Zhang", "Beck", "Elmore", "Gu")[-5], #,"klos",
    verbose = F,
    wFUN = wFUN,
    nextend = 2, maxExtendMonth = 2, minExtendMonth = 1,
    minPercValid = 0.2,
    print = FALSE
)

julia_init()
fit_R      = do.call(curvefits, c(param, list(use.julia = FALSE)))
fit_julia  = do.call(curvefits, c(param, list(use.julia = TRUE)))

bench::mark(
    fit_R     <- do.call(curvefits, c(param, list(use.julia = FALSE))),
    fit_Julia <- do.call(curvefits, c(param, list(use.julia = TRUE))),
    min_time = 60
)

profvis::profvis(
    for (i in 1:10) {
        fit_R     = do.call(curvefits, c(param, list(use.julia = FALSE)))
        fit_Julia <- do.call(curvefits, c(param, list(use.julia = TRUE)))
    }
)


## 3. Test for parallel performance
# method = "BFGS"
# par0 <- c(t0 = 100, par0)
# rbenchmark::benchmark(
#   r1 = optim_pheno(par0, sFUN, y, t, tout, methods[4],
#               nptperyear = 46, iters = 2, wFUN = wTSM, verbose = FALSE)
# )
