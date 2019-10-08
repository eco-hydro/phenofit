context("test-optim")

library(ggplot2)
library(magrittr)
library(purrr)

# simulate vegetation time-series
fFUN = doubleLog_Beck
par  = c( mn  = 0.1 , mx  = 0.7 , sos = 50 , rsp = 0.1 , eos = 250, rau = 0.1)
par0 = c( mn  = 0.15, mx  = 0.65, sos = 100, rsp = 0.12, eos = 200, rau = 0.12)

t    <- seq(1, 365, 8)
tout <- seq(1, 365, 1)
y    <- fFUN(par, t)
ypred<- t*0

optFUNs <- c("opt_ucminf", "opt_nlminb", "opt_nlm", "opt_optim") # %>% set_names(., .)

optim_fFUN <- function(optFUN, objective){
    # print(optFUN)
    optFUN <- get(optFUN)

    opt <- suppressMessages( optFUN(par0, objective, y = y, t = t, pred = ypred, fun = doubleLog_Beck) )
    opt$ysim <- fFUN(opt$par, t)
    opt
}

test_that("optFUNs works", {
    # 1. test for error situation
    # f_goal # goal function
    # opts <- lapply(optFUNs, optim_fFUN, objective = NULL)
    # convcode <- map_dbl(opts, "convcode")
    # expect_equal(convcode, rep(9999, length(optFUNs)))

    # 2. test for normal situation
    opts <- lapply(optFUNs, optim_fFUN, f_goal)
    convcode <- map_dbl(opts, "convcode")
    expect_true(all(convcode <= 1))
})


test_that("I_optim works", {
    prior <- as.matrix(par0) %>% t() %>% rbind(., .)

    opt1 <- I_optim(prior, fFUN, y, t, c("BFGS", "ucminf", "nlm", "nlminb"))
    opt2 <- I_optimx(prior, fFUN, y, t, c("BFGS", "ucminf", "nlm", "nlminb"))

    expect_equal(nrow(opt1), 8)
    expect_equal(nrow(opt2), 8)
})
