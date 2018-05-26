
plot_qc <- function(t, y, w, ...){
    # factor w
    # wf <- findInterval(w, c(0, 0.5, 1, Inf), include.lowest = T, right = F)
    wf <- 4 - findInterval(w, c(0, 0.5, 1, Inf))
    
    colors <- c("grey60", "#00BFC4", "#F8766D", "#C77CFF", "blue", "red", "black")
    pch <- c(19, 15, 4)
    
    plot(t, y, type = "l", ...)
    Ids <- unique(wf)
    for (i in 1:3){
        I = wf == i
        add <- ifelse(i == 1, F, T)
        points(t[I], y[I], pch = pch[i], col = colors[i], cex = 0.8)
    }
    
    abline(v = t[seq(1, length(y), nptperyear)], col = "grey60", lty = 3)
}
