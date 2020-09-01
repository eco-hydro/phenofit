set.seed(1000)

n = 1e5
y = rnorm(1e5)
system.time({
    for (i in 1:1e3) {
        z = rcpp_SG(y, 7, 2)
    }
})
