library(phenofit)
# simulate vegetation time-series
FUN = doubleLog.Beck
par  = c( mn  = 0.1, mx  = 0.7, sos = 50, rsp = 0.1, eos = 250, rau = 0.1)
t    <- seq(1, 365, 8)
tout <- seq(1, 365, 1)
y <- FUN(par, t)
methods <- c("AG", "Beck", "Elmore", "Gu", "Zhang") # "Klos" too slow
fit <- curvefit(y, t, tout, methods)
