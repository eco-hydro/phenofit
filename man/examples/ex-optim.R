# test the computing performance of `I_optim` function
# library(ggplot2)
library(magrittr)
library(purrr)

# simulate vegetation time-series
fFUN = doubleLog.Beck
par = c(
  mn  = 0.1,
  mx  = 0.7,
  sos = 50,
  rsp = 0.1,
  eos = 250,
  rau = 0.1)
t    <- seq(1, 365, 8)
tout <- seq(1, 365, 1)
y <- fFUN(par, t)

# initial parameter
par0 <- c(
  mn  = 0.15,
  mx  = 0.65,
  sos = 100,
  rsp = 0.12,
  eos = 200,
  rau = 0.12)

objective <- f_goal # goal function
prior <- as.matrix(par0) %>% t() %>% rbind(., .)

opt1 <- I_optim(prior, fFUN, y, t, tout, c("nlminb")) # "BFGS", "ucminf", "nlm",
# opt2 <- I_optimx(prior, fFUN, y, t, tout, c("BFGS", "ucminf", "nlm", "nlminb"))

## Not run:
# mean: 130 ms
methods = c("BFGS", "ucminf", "nlm", "nlminb")
microbenchmark::microbenchmark(
  I_optim(prior, fFUN, y, t, tout, methods, fn = phenofit::f_goal),
  I_optim(prior, fFUN, y, t, tout, methods, fn = phenofit::f_goal_r),
  # I_optimx(prior, fFUN, y, t, tout, c("BFGS", "ucminf", "nlm", "nlminb")),
  times = 20
)
