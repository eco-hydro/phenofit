library(phenofit)
library(ggplot2)
library(magrittr)
library(purrr)
library(data.table)

# simulate vegetation time-series
fFUN = doubleLog_Beck
par  = c( mn  = 0.1 , mx  = 0.7 , sos = 50 , rsp = 0.1 , eos = 250, rau = 0.1)
par0 = c( mn  = 0.15, mx  = 0.65, sos = 100, rsp = 0.12, eos = 200, rau = 0.12)

t    <- seq(1, 365, 8)
tout <- seq(1, 365, 1)
y <- fFUN(par, t)

optFUNs <- c("opt_ucminf", "opt_nlminb", "opt_nlm", "opt_optim") %>% set_names(., .)
opts <- lapply(optFUNs, function(optFUN){
    optFUN <- get(optFUN)
    opt    <- optFUN(par0, f_goal, y = y, t = t, fun = fFUN)
    opt$ysim <- fFUN(opt$par, t)
    opt
})

# visualization
df   <- map(opts, "ysim") %>% as.data.table() %>% cbind(t, y, .)
pdat <- data.table::melt(df, c("t", "y"), variable.name = "optFUN")

ggplot(pdat) + 
    geom_point(data = data.frame(t, y), aes(t, y), size = 2) + 
    geom_line(aes(t, value, color = optFUN), linewidth = 0.9)
