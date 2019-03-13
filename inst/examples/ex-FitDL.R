library(phenofit)
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
methods <- c("AG", "Beck", "Elmore", "Gu", "Zhang")

r <- FitAG(y, t, tout)
plot(t, y)
lines(tout, r$zs$iter2, col = "red")
legend('topright', c('Original time-series', 'AG smoothed'), 
    lty = c(0, 1), pch = c(16, NA), col = c("black", "red"))
