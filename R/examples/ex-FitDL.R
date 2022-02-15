# simulate vegetation time-series
t    <- seq(1, 365, 8)
par  <- c(mn = 0.1, mx = 0.7, sos = 50, rsp = 0.1, eos = 250, rau = 0.1)
y    <- doubleLog.Beck(par, t)
data <- data.frame(t, y)
# methods <- c("AG", "Beck", "Elmore", "Gu", "Zhang")
tout <- seq(1, 365, 1)
r <- FitDL.Elmore(y, t, tout)

plot(r, data)
get_GOF(r, data)
get_param(r)
