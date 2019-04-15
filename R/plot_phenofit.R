# update to support MOD09A1, also suit for previous MOD13A1

#' plot_phenofit
#'
#' @param d_fit data.frame of curve fittings returned by \code{\link{get_fitting}}.
#' @param seasons Growing season dividing object returned by \code{\link{season}}
#' and \code{\link{season_mov}}.
#' @param title String, title of figure.
#' @param font.size Font size of axis.text
#' @param show.legend Boolean
#'
#' @importFrom dplyr left_join
#'
#' @example inst/examples/ex-check_input.R
#' @example inst/examples/ex-visual.R
#'
#' @export
plot_phenofit <- function(d_fit,
                          seasons,
                          title = NULL,
                          title.ylab = "Vegetation Index",
                          font.size = 14,
                          show.legend = TRUE)
{
    methods <- d_fit$meth %>% table() %>% names()
    nmethod <- length(methods) # how many curve fitting methods?

    d_obs <- d_fit[meth == meth[1]] # t, y
    d_obs$meth <- NULL

    last_iter_rough <- colnames(seasons$whit) %>% last()

    iters_name_fine <- colnames(d_obs) %>% .[grep("ziter", .)] %>% sort() # "ziter1", "ziter2"
    lines_colors     <- iter_colors(length(iters_name_fine)) %>% set_names(iters_name_fine)# only for smoothed time-series

    # browser()
    nyear <- diff(range(year(d_obs$t), na.rm = TRUE)) + 1
    nyear_lean <- 7 # more than `nyear_lean`, then lean axis.text.x 30deg

    # pdat    <- get_fitting(fit)
    p <- ggplot(d_fit, aes_string("t", "y", color = "iters")) +
        geom_line (data = seasons$whit, aes_string("t", last_iter_rough), color = "black", size = 0.8) + # show in front
        geom_vline(data = seasons$dt, aes(xintercept = as.numeric(beg)), size = 0.4, linetype=2, color = "blue") +
        geom_vline(data = seasons$dt, aes(xintercept = as.numeric(end)), size = 0.4, linetype=2, color = "red") +
        # geom_point(data = seasons$dt, aes_string(peak, y_peak), color= "red") +
        # geom_point(data = seasons$dt, aes_string(beg , y_beg ), color= "blue") +
        # geom_line(size = 0.4) +
        facet_grid(meth~.) +
        scale_x_date(breaks = seasons$dt$beg, date_labels = "%Y/%m") +
        ggtitle(title) +
        theme_gray(base_size = font.size) +
        theme(axis.title = element_text(size = font.size),
            # axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
            axis.text = element_text(size = font.size - 2)) +
        labs(x = 'Time', y = title.ylab)

    if (nyear >= nyear_lean) {
        p <- p + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
    }

    if ('QC_flag' %in% colnames(d_obs)){
        #     guides(shape = guide_legend(override.aes_string = list(size = 2)))
        p  <- p + geom_point(data = d_obs, aes_string("t", "y", shape="QC_flag",
                                                  color = "QC_flag", fill = "QC_flag"),
                             size = 2, alpha = 0.7)
    } else {
        p <- p + geom_point(aes_string("t", "y"), size = 2, alpha = 0.5, color = "grey60")
            # geom_line (data = seasons$whit, aes_string(t, ziter2), color = "black", size = 0.8) + # show in front
    }

    # iterations of smoothed time-series
    for (i in seq_along(iters_name_fine)) {
        iter_name <- iters_name_fine[i]
        p <- p + geom_line(aes_string(y = iter_name), size = 0.8, alpha = 0.7, color = lines_colors[i])
    }

    p <- p + scale_color_manual(values = c(qc_colors, lines_colors), drop = F) +
        scale_fill_manual(values = qc_colors, drop = F) +
        scale_shape_manual(values = qc_shapes, drop = F)

    # geom_vline(data = pdat2, aes_string(xintercept=date, linetype = pmeth, color = pmeth), size = 0.4, alpha = 1) +
    # scale_linetype_manual(values=c(2, 3, 1, 1, 1, 4))
    # p + facet_grid(meth~pmeth)

    # if (plotly){
    #     plotly::ggplotly(p)
    # }else{

    if (show.legend){
        lgd <- make_legend_nmax(iters_name_fine, lines_colors, d_obs$QC_flag)

        p <- p + theme(legend.position="none")
        p <- arrangeGrob(p, lgd, nrow = 2, heights = c(min(5*nmethod, 15), 1),
            padding = unit(1, "line"))
    }
    return(p)
}
