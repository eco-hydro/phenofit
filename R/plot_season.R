#' plot_season
#'
#' Plot growing season divding result.
#'
#' @inheritParams season
#' @inheritParams wHANTS
#' @param brks A list object returned by `season` or `season_mov`.
#' @param show.legend Whether to show legend?
#'
#' @importFrom grid viewport pushViewport grid.draw
#' @export
plot_season <- function(INPUT, brks, plotdat, ylu,
                        IsPlot.OnlyBad = FALSE, show.legend = TRUE){
    if (missing(plotdat)) {
        plotdat <- INPUT
    }
    
    stat <- stat_season(INPUT, brks)
    stat_txt  <- stat[c("R2", "NSE", "sim.cv", "obs.cv")] %>% unlist() %>%
        set_names(c("R2", "NSE", "CV_sim", "CV_obs")) %>%
        # set_names(c("R2", "NSE", expression(CV[sim]), expression(CV[obs]))) %>%
        {paste(names(.), round(., 3), sep = "=", collapse = ", ")}

    # if (NSE < 0 | (cv < 0.1 & NSE < 0.1)) {}
    # if(IsPlot && (NSE < 0 && cv < 0.2)){
    if (IsPlot.OnlyBad && stat['NSE'] < 0.3) return()

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
    # par.old <- par()
    par(mar = c(3.5, 3, 1, 1), mgp = c(1.2, 0.6, 0))
    plot_input(plotdat)

    # colors <- c("blue", "red", "green")
    NITER  <- ncol(zs)
    # if (NITER < 3) colors <- c("blue", "red")
    lines_colors <- iter_colors(NITER)

    for (i in 1:NITER){
        lines(t, zs[[i]], col = lines_colors[i], lwd = 2)
    }

    # 7.2 plot break points
    points(dt$peak, dt$y_peak, pch=20, cex = 1.8, col="red")
    points(dt$beg , dt$y_beg , pch=20, cex = 1.8, col="blue")
    points(dt$end , dt$y_end , pch=20, cex = 1.8, col="blue")

    if (!missing(ylu)) abline(h=ylu, col="red", lty=2) # show ylims
    legend('topleft', stat_txt, adj = c(0.05, 0.8), bty='n', text.col = "red")

    if (show.legend){
        lgd <- make_legend_nmax(paste0("iter", 1:NITER), lines_colors, plotdat$QC_flag)
        # fix in the futher
        pos.y <- 0.04
        pushViewport(viewport(x = unit(0.5, "npc"), y = unit(pos.y, "npc"),
                              width = unit(0.8, "npc"), height = unit(pos.y*2, "npc")))
        # grid.rect()
        grid.draw(lgd)
    }
}
