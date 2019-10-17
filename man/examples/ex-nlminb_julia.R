\dontrun{

t    = seq(1.0, 366, 8)
fun  = doubleLog_Beck
par  = c(0.1 , 0.7, 50, 0.1, 250, 0.1)
par0 = c(0.05, 0.6 , 45, 0.1, 200, 0.2)

ypred = t*0
y     = fun(par, t)

julia_init()
r_julia  <- opt_nlminb_julia(par0, "doubleLog_Beck", y, t)
r_R <- opt_nlminb(par0, f_goal, fun = fun, y = y, t = t, pred = ypred)

list(julia = r_julia, R = r_R) %>%
    map(~c(.$par, .$objective, .$value)) %>%
    do.call(rbind, .)# %>%

n <- length(t)
w <- rep(0.2, n)
# julia is 5 times faster
{
    # microbenchmark::microbenchmark : 18.939826 ms in R
    info <- rbenchmark::benchmark(
        r1 <- opt_nlminb_julia(par0, "doubleLog_Beck", y, t, w),
        r2 <- opt_nlminb(par0, f_goal, fun = fun, y = y, t = t, pred = ypred),
        replications = 500
    )
    print(info)
}

}
