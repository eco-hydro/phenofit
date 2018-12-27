#' plot_season
#' 
#' @param brks A list object returned by \code{season_3y}.
#' @inheritParams check_fit
#' @rdname season
#' 
#' @export
plot_season <- function(INPUT, brks, plotdat, ylu, IsOnlyPlotbad = FALSE){
    stat <- stat_season(INPUT, brks)
    stat_txt  <- stat[c("R2", "NSE", "sim.cv", "obs.cv")] %>% unlist() %>% 
        {paste(names(.), round(., 3), sep = "=", collapse = ", ")}

    # if (NSE < 0 | (cv < 0.1 & NSE < 0.1)) {}
    # if(IsPlot && (NSE < 0 && cv < 0.2)){
    if (IsOnlyPlotbad && stat['NSE'] < 0.3) return()

    t  <- brks$whit$t
    dt <- brks$dt
    zs <- dplyr::select(brks$whit, dplyr::matches("ziter*"))
    ypred <- last(zs)

    # if (missing(xlim))
    xlim  <- c(first(brks$dt$beg), last(brks$dt$end))

    ## PLOT CURVE FITTING TIME-SERIES
    #  need to plot outside, because y, w have been changed.
    # plotdat   <- INPUT[c("t", "y", "w", "ylu")]
    # if (!is.null(INPUT$y0)) plotdat$y <- INPUT$y0
    plotdata(plotdat)

    colors <- c("red", "blue", "green")
    iters  <- ncol(zs)
    if (iters < 3) colors <- c("red", "green")

    for (i in 1:iters){
        lines(t, zs[[i]], col = colors[i], lwd = 2)
    }

    # 7.2 plot break points
    points(dt$peak, dt$y_peak, pch=20, cex = 1.8, col="red")
    points(dt$beg , dt$y_beg , pch=20, cex = 1.8, col="blue")
    points(dt$end , dt$y_end , pch=20, cex = 1.8, col="blue")

    if (!missing(ylu)) abline(h=ylu, col="red", lty=2) # show ylims
    legend('topleft', stat_txt, adj = c(0.05, 0.8), bty='n', text.col = "red")
}
