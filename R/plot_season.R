#' plot_season
#'
#' Plot growing season divding result.
#'
#' @inheritParams season
#' @inheritParams smooth_wHANTS
#' @param brks A list object returned by `season` or `season_mov`.
#' @param show.legend Whether to show legend?
#' @param title The main title (on top)
#'
#' @importFrom grid viewport pushViewport grid.draw
#' @export
plot_season <- function(
    INPUT, brks, plotdat, ylu,
    IsPlot.OnlyBad = FALSE, show.legend = TRUE,
    ylab = "VI", title = NULL)
{
    if (missing(plotdat)) {
        plotdat <- INPUT
    }
    if (is.data.frame(brks$dt[[1]])) {
        brks$dt %<>% do.call(rbind, .)
    }
    ylu = INPUT$ylu
    A = diff(ylu)
    # ylu2 = c(ylu[1], ylu[2] + 0.1*A)

    stat <- stat_season(INPUT, brks)
    stat_txt <- eval(substitute(
        expression(
            R^2 == R2 * ", NSE = " * NSE * ", " * CV[sim] ~ "= " * sim.cv *
                ", " * CV[obs] ~ "= " * obs.cv),
        stat[c("R2", "NSE", "sim.cv", "obs.cv")] %>% round(3) %>% as.list()
    ))
    if (IsPlot.OnlyBad && stat['NSE'] < 0.3) return()

    t  <- brks$fit$t
    dt <- brks$dt
    zs <- select_vars(brks$fit, "ziter")
    ypred <- last(zs)

    # if (missing(xlim))
    xlim  <- c(first(brks$dt$beg), last(brks$dt$end))

    ## PLOT CURVE FITTING TIME-SERIES
    #  need to plot outside, because y, w have been changed.
    # plotdat   <- INPUT[c("t", "y", "w", "ylu")]
    # if (!is.null(INPUT$y0)) plotdat$y <- INPUT$y0
    # par.old <- par()
    bottom = ifelse(show.legend, 3.2, 1) + 0.5
    par.old <- par(mar = c(bottom, 3, 1, 1), mgp = c(1.2, 0.6, 0))
    on.exit(par(par.old))
    
    ymin = min(INPUT$y0, na.rm = TRUE)
    ymax = max(INPUT$y0, na.rm = TRUE)
    ylim = c(ymin, (ymax - ymin)*0.08 + ymax)
    plot_input(plotdat, ylab = ylab, ylim = ylim)
    plot_season_boundary(dt)

    NITER  <- ncol(zs)
    lines_colors <- iter_colors(NITER)
    if (NITER == 1) lines_colors <- c("red")

    for (i in 1:NITER){
        lines(t, zs[[i]], col = lines_colors[i], lwd = 2)
    }

    # 7.2 plot break points
    points(dt$peak, dt$y_peak, pch=20, cex = 1.8, col="red")
    points(dt$beg , dt$y_beg , pch=20, cex = 1.8, col="blue")
    points(dt$end , dt$y_end , pch=20, cex = 1.8, col="blue")

    if (!missing(ylu)) abline(h=ylu, col="red", lty=2) # show ylims
    legend('topleft', stat_txt, adj = c(0.05, 0.4), bty='n', text.col = "red")

    if (!is.null(title)) title(title, line = -1)
    if (show.legend){
        lgd <- make_legend_nmax(paste0("iter", 1:NITER), lines_colors, plotdat$QC_flag)
        # fix in the futher
        pos.y <- 0.032
        pushViewport(viewport(x = unit(0.5, "npc"), y = unit(pos.y, "npc"),
                              width = unit(0.8, "npc"), height = unit(pos.y*2, "npc")))
        # grid.rect()
        grid.draw(lgd)
    }
}

plot_season_boundary <- function(dt) {
    ylim = par("usr")[3:4]
    A    = diff(ylim)
    ylim = ylim + c(1, -1)*A/200

    nGS = nrow(dt)
    if (!is.null(nGS) && nGS > 0) {
        xs = lapply(1:nrow(dt), function(i){
            c(dt$beg[i], dt$end[i]) %>% c(., rev(.), NA)
        }) %>% do.call(c, .)
        ys = rep(ylim, each = 2) %>% c(NA) %>% rep(nGS)
        polygon(xs, ys,
                # density = 2, angle = 30,
                col = alpha("grey", 0.2),
                border = alpha("grey", 0.4))
    }
}
