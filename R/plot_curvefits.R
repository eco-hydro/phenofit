#' plot_curvefits
#'
#' @param d_fit data.frame of curve fittings returned by [get_fitting()].
#' @param seasons Growing season dividing object returned by [season()]
#' and [season_mov()].
#' @param d_obs data.frame of original vegetation time series, with the columns
#' of `t`, `y` and `QC_flag`. If not specified, it will be determined from `d_fit`.
#' @param title String, title of figure.
#' @param xlab,ylab String, title of `xlab` and `ylab`.
#' @param yticks ticks of y axis
#' @param font.size Font size of axis.text
#' @param show.legend Boolean
#' @param theme ggplot theme to be applied
#' @param shape the shape of input VI observation? `line` or `point`
#' @param cex point size for VI observation.
#' @param angle `text.x` angle
#'
#' @example inst/examples/ex-check_input.R
#' @example inst/examples/ex-visual.R
#'
#' @export
plot_curvefits <- function(
    d_fit,
    seasons,
    d_obs = NULL,
    title = NULL,
    xlab = "Time", ylab = "Vegetation Index",
    yticks = NULL,
    font.size = 14,
    theme = NULL,
    cex = 2,
    shape = "point", angle = 30,
    show.legend = TRUE)
{
    methods <- d_fit$meth %>% unique() %>% rm_empty() # in case of NA
    nmethod <- length(methods) # how many curve fitting methods?

    if (is.null(d_obs)) {
        d_obs <- d_fit[meth == methods[1]] # t, y
        d_obs$meth <- NULL
    }

    last_iter_rough <- colnames(seasons$fit) %>% last()
    iters_name_fine <- colnames(d_fit) %>% .[grep("ziter", .)] %>% sort() # "ziter1", "ziter2"
    lines_colors <- iter_colors(length(iters_name_fine)) %>% set_names(iters_name_fine) # only for smoothed time-series

    nyear <- diff(range(year(d_obs$t), na.rm = TRUE)) + 1
    nyear_lean <- 8 # more than `nyear_lean`, then lean axis.text.x 30deg
    # pdat    <- get_fitting(fit)
    p <- ggplot(d_fit, aes_string("t", "y", color = "iters")) +
        geom_vline(data = seasons$dt, aes(xintercept = beg), size = 0.4, linetype = 2, color = "blue") +
        geom_vline(data = seasons$dt, aes(xintercept = end), size = 0.4, linetype = 2, color = "red") +
        facet_grid(meth ~ .) +
        # scale_x_date(breaks = seasons$dt$beg, date_labels = "%Y/%m") +
        ggtitle(title) +
        # theme_gray(base_size = font.size) +
        theme_bw(base_size = font.size) +
        theme(
            axis.title = element_text(size = font.size),
            # axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
            axis.text = element_text(size = font.size - 2),
            plot.margin = margin(0.2, 0.2, -0.3, 0.2, "lines")
        ) +
        labs(x = xlab, y = ylab)
    if (nyear >= nyear_lean) {
        p <- p + theme(axis.text.x = element_text(angle = angle, hjust = 1, vjust = 1))
    }

    if ("QC_flag" %in% colnames(d_obs)) {
        p <- p + geom_point(
            data = d_obs,
            aes_string("t", "y", shape = "QC_flag", color = "QC_flag", fill = "QC_flag"),
            size = cex, alpha = 0.7
        )
    } else {
        p <- if (shape == "point") {
            p + geom_point(data = d_obs, aes_string("t", "y"), size = cex, alpha = 0.6, color = "grey60")
        } else if (shape == "line") {
            p + geom_line(data = d_obs, aes_string("t", "y"), size = 0.8, alpha = 0.8, color = "grey60")
        }
    }

    # iterations of smoothed time-series
    p <- p + geom_line(data = seasons$fit, aes_string("t", last_iter_rough), color = "black", size = 0.8)
    for (i in seq_along(iters_name_fine)) {
        iter_name <- iters_name_fine[i]
        p <- p + geom_line(aes_string(y = iter_name), size = 0.8, alpha = 0.7, color = lines_colors[i])
    }

    xlim <- d_fit$t %>% range() %>% add(ddays(c(1, -1) * 90))
    p <- p + scale_color_manual(values = c(qc_colors, lines_colors), drop = F) +
        scale_fill_manual(values = qc_colors, drop = F) +
        scale_shape_manual(values = qc_shapes, drop = F) +
        coord_cartesian(xlim = xlim)

    if (!is.null(theme)) p <- p + theme
    if (!is.null(yticks)) p <- p + scale_y_continuous(breaks = yticks)
    if  (is.null(title)) p <- p + theme(plot.title = element_blank())

    if (show.legend) {
        iters_name_fine = c("Rough fitting", "Fine fitting")
        lines_colors = c("black", "red")
        lgd <- make_legend_nmax(iters_name_fine, lines_colors, d_obs$QC_flag)
        p <- p + theme(legend.position = "none")
        p <- arrangeGrob(p, lgd,
            nrow = 2, heights = c(min(5 * nmethod, 15), 1),
            padding = unit(1, "line")
        )
    }
    return(p)
}
# ggnewscale::new_scale_color
