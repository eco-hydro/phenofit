#' Plot INPUT returned by check_input
#' 
#' @inheritParams season
#' @param ... other parameter will be ignored.
#' @export
plotdata <- function(INPUT, wmin = 0.1, ...){
    t <- INPUT$t
    y <- INPUT$y
    w <- INPUT$w
    # nptperyear <- INPUT$nptperyear
    
    npt <- length(y)
    # show grid lines
    par(mgp = c(1.5, 0.5, 0)) #oma = c(1, 2, 3, 1)
    # at <- t[seq(1, npt, nptperyear)]
    # fmt <- ifelse(yday(at[1]) == 1, "%Y", "%Y/%m/%d")
    # axis(side=1, at = at, labels = format(at, fmt))
    
    wf <- 4 - findInterval(w, c(-Inf, wmin, 0.5, 1), left.open = T)
    
    colors <- c("grey60", "#00BFC4", "#F8766D", "#C77CFF", "blue", "red", "black")
    pch    <- c(19, 15, 4)
    
    years    <- year(t)

    if (length(unique(years)) < 3){
        plot(t, y, type = "l", xaxt="n", ann = FALSE, ...)
        axis.Date(1, at=seq(min(t), max(t), by="month"), format="%Y-%m")
    } else {
        plot(t, y, type = "l", ann = FALSE, ...)
    }
    
    Ids <- unique(wf)
    for (i in 1:3){
        I = wf == i
        add <- ifelse(i == 1, F, T)
        points(t[I], y[I], pch = pch[i], col = colors[i], cex = 0.6)
    }

    # ylab = expression(paste('GPP ( gC ', m^-2, d^-1, ')'))
    date_beg <- ymd( min(years) *1e4 + 0101 )
    date_end <- ymd( max(years) *1e4 + 0101 )
    
    t_grids  <- seq.Date(date_beg, date_end, by = "year")
    abline(v = t_grids, col = "grey60", lty = 3)
    grid(nx = NA, NULL)
    ylu <- INPUT$ylu
    if (!is.null(ylu)) abline(h=ylu, col="red", lty=2) # show ylims
}
