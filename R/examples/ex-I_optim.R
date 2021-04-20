# simulate vegetation time-series
FUN  = doubleLog_Beck
par  = c( mn  = 0.1 , mx  = 0.7 , sos = 50 , rsp = 0.1 , eos = 250, rau = 0.1)
par0 = c( mn  = 0.15, mx  = 0.65, sos = 100, rsp = 0.12, eos = 200, rau = 0.12)

t <- seq(1, 365, 8)
y <- FUN(par, t)

methods = c("BFGS", "ucminf", "nlm", "nlminb")
opt1 <- I_optim (par0, FUN, y, t, methods)
opt2 <- I_optimx(par0, FUN, y, t, methods)

# \dontrun{
# microbenchmark::microbenchmark(
#     opt1 = I_optim (par0, FUN, y, t, methods),
#     opt2 = I_optimx(par0, FUN, y, t, methods),
#     times = 2
# )
# }
