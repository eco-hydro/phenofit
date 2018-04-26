f <- function(x,y,groups,...){
    s2  <- smooth.spline(x, y, df = freedom[panel.number()])
    # s2  <- smooth.spline(x, y, lambda = freedom[panel.number()])

    panel.xyplot(x, y, groups,...)
    panel.lines(s2)
    xpos <- range(x) %>% {.[1] + diff(.)/10}
    ypos <- range(y) %>% {.[1] + diff(.)/10}
    panel.text(xpos, ypos, adj = c(0,1), sprintf("cv.crit = %.1e\npen.crit = %.3e", s2$cv.crit, s2$pen.crit))
}

library(lattice)

freedom <- seq(2, min(length(yi), 40), length.out = 20) %>% floor
df <- data.frame(yi, ti, freedom = as.factor(rep(freedom, each =length(yi))))

xyplot(yi~ti|freedom, data = df, type=c('p'), pch=19, groups=freedom,
       panel = f, as.table = T)
