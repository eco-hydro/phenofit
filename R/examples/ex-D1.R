# doubleLog.Beck
t    <- seq(1, 365, 8)
tout <- seq(1, 365, 1)
par  = c(mn  = 0.1, mx  = 0.7, sos = 50, rsp = 0.1, eos = 250, rau = 0.1)
y <- doubleLog.Beck(par, t)

methods <- c("AG", "Beck", "Elmore", "Gu", "Zhang") # "Klos" too slow
fit <- curvefit(y, t, tout, methods)
x  <- fit$model$AG
d1 <- D1(x)
d2 <- D2(x)
d_k <- curvature(x)
