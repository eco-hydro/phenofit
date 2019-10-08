# library(magrittr)
# library(purrr)

# simulate vegetation time-series
FUN = doubleLog.Beck
par  = c( mn  = 0.1 , mx  = 0.7 , sos = 50 , rsp = 0.1 , eos = 250, rau = 0.1)
par0 = c( mn  = 0.15, mx  = 0.65, sos = 100, rsp = 0.12, eos = 200, rau = 0.12)

t    <- seq(1, 365, 8)
tout <- seq(1, 365, 1)
y <- FUN(par, t)

opt1 <- I_optim(par0, FUN, y, t, tout, c("nlminb")) # "BFGS", "ucminf", "nlm",
# opt2 <- I_optimx(prior, fFUN, y, t, tout, )

methods = c("BFGS", "ucminf", "nlm", "nlminb")[4]
r <- optim_pheno(par0, "doubleLog.Beck", y, t, tout, method = methods,
                 nptperyear = 46, iters = 2, wFUN = wTSM,
                 use.cpp = FALSE, verbose = FALSE)
# method = "BFGS"
# sFUN   = "doubleLog.Zhang" # doubleLog.Beck
# par0 <- c(t0 = 100, par0)
# rbenchmark::benchmark(
#   r1 = optim_pheno(par0, sFUN, y, t, tout, method ,
#               nptperyear = 46, iters = 2, wFUN = wTSM,
#               use.cpp = FALSE, verbose = FALSE),
#   r2 = optim_pheno(par0, sFUN, y, t, tout, method,
#               nptperyear = 46, iters = 2, wFUN = wTSM,
#               use.cpp = TRUE, verbose = FALSE)
# )
