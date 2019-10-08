par  = c( mn  = 0.1 , mx  = 0.7 , sos = 50 , rsp = 0.1 , eos = 250, rau = 0.1)
par0 = c( mn  = 0.15, mx  = 0.65, sos = 100, rsp = 0.12, eos = 200, rau = 0.12)

t    = seq(1, 365, 8)
tout = seq(1, 365, 1)
y    = doubleLog.Beck(pars[-1], t)


# methods = c("nlminb")
methods = c("BFGS", "ucminf", "nlm", "nlminb")

{
  # beck
  prior = par0[-1]
  opt1 <- I_optim(prior, doubleLog.Beck, y, t, tout, methods)

  # 1.6 times faster for Beck
  rbenchmark::benchmark(
    I_optim(prior, doubleLog_Beck, y, t, tout, methods, fn = f_goal),
    I_optim(prior, doubleLog_Beck, y, t, tout, methods, fn = f_goal),
    # I_optim(prior, FUN = doubleLog_Beck, y, t, tout, methods, fn = phenofit::f_goal_r),
    # I_optimx(prior, fFUN, y, t, tout, c("BFGS", "ucminf", "nlm", "nlminb")),
    replications = 2
  )
}

# {
#   # beck
#   prior = par0[-1]
#   opt1 <- I_optim(prior, doubleLog.Zhang, y, t, tout, methods)

#   # 1.6 times faster for Beck
#   rbenchmark::benchmark(
#     I_optim(prior, doubleLog.Zhang, y, t, tout, methods, fn = f_goal),
#     I_optim(prior, doubleLog_Zhang, y, t, tout, methods, fn = f_goal),
#     # I_optim(prior, FUN = doubleLog_Beck, y, t, tout, methods, fn = phenofit::f_goal_r),
#     # I_optimx(prior, fFUN, y, t, tout, c("BFGS", "ucminf", "nlm", "nlminb")),
#     replications = 10
#   )
# }

# test_that("Rcpp double logistics works", {
#     expect_equal(2 * 2, 4)
# })
