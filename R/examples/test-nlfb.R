# 1. nlfr package
start <- set_names(prior[1, ], attr(FUN, "par"))
# lower = lower, upper = upper
f_nlf <- nlfb(start, resfn = function(par, t, y, fun, ...){
    fun(par, t) - y
}, trace = FALSE, y = y, t=t, w = w, fun = FUN, ...)
yfit  <- FUN(f_nlf$coefficients, t)
plot(t, y, type = "b")
lines(t, yfit)
print(f_nlf)
