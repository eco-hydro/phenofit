# library(magrittr)
# library(purrr)

# simulate vegetation time-series
FUN = doubleLog_Beck
par  = c( mn  = 0.1 , mx  = 0.7 , sos = 50 , rsp = 0.1 , eos = 250, rau = 0.1)
par0 = c( mn  = 0.15, mx  = 0.65, sos = 100, rsp = 0.12, eos = 200, rau = 0.12)

t    <- seq(1, 365, 8)
tout <- seq(1, 365, 1)
y <- FUN(par, t)

methods = c("BFGS", "ucminf", "nlm", "nlminb")
opt1 <- I_optim(par0, doubleLog_Beck, y, t, methods) # "BFGS", "ucminf", "nlm",
# opt2 <- I_optimx(prior, fFUN, y, t, tout, )

sFUN   = "doubleLog.Beck" # doubleLog.Beck
r <- optim_pheno(par0, sFUN, y, t, tout, method = methods[4],
                 nptperyear = 46, iters = 2, wFUN = wTSM, verbose = FALSE)
# method = "BFGS"
# par0 <- c(t0 = 100, par0)
# rbenchmark::benchmark(
#   r1 = optim_pheno(par0, sFUN, y, t, tout, methods[4],
#               nptperyear = 46, iters = 2, wFUN = wTSM, verbose = FALSE)
# )
