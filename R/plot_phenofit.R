# update to support MOD09A1, also suit for previous MOD13A1

#' plot_phenofit
#'
#' @param df_fit data.frame of curve fittings returned by \code{\link{get_fitting}}.
#' @param seasons Growing season dividing object returned by \code{\link{season}}
#' and \code{\link{season_3y}}.
#' @param title String, title of figure.
#' @param font.size Font size of axis.text
#' @param show.legend Boolean
#'
#' @importFrom dplyr left_join
#' @export
plot_phenofit <- function(df_fit,
                          seasons,
                          title = NULL,
                          font.size = 14,
                          show.legend = TRUE)
{
    methods <- df_fit$meth %>% table() %>% names()
    nmethod <- length(methods) # how many curve fitting methods?
    
    d_obs <- df_fit[meth == meth[1]] # t, y
    d_obs$meth <- NULL
    
    # pdat    <- get_fitting(fit)
    # seasons <- fit$seasons
    p <- ggplot(df_fit, aes_string("t", "y", color = "iters")) +
        geom_line (data = seasons$whit, aes_string("t", "ziter2"), color = "black", size = 0.8) + # show in front
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
            axis.text = element_text(size = font.size - 2))

    if ('QC_flag' %in% colnames(d_obs)){
        # p <- p + geom_point(data = x, aes_string(date, y, shape = SummaryQA, color = SummaryQA), size = 1, alpha = 0.7) +
        #     scale_shape_manual(values = c(21,22, 24:25)) +
        #     scale_fill_manual(values = c("grey40", "#7CAE00", "#F8766D", "#C77CFF")) +
        #     guides(shape = guide_legend(override.aes_string = list(size = 2)))
        p  <- p + geom_point(data = d_obs, aes_string("t", "y", shape="QC_flag",
                                                  color = "QC_flag", fill = "QC_flag"),
                             size = 2, alpha = 0.7) +
            # geom_line (data = seasons$whit, aes_string(t, ziter2), color = "black", size = 0.8) + # show in front
            geom_line(aes_string(y = "ziter1"), size = 0.8, alpha = 0.7, color = "blue") +
            geom_line(aes_string(y = "ziter2"), size = 0.8, alpha = 0.7, color = "red") +
            scale_color_manual(values = c(qc_colors,"iter1" = "blue", "iter2" = "red"), drop = F) +
            scale_fill_manual(values = qc_colors, drop = F) +
            scale_shape_manual(values = qc_shapes, drop = F) +
            ylab('Vegetation Index')
    }else{
        p <- p + geom_point(aes_string("t", "y"), size = 2, alpha = 0.5, color = "grey60") +
            # geom_line (data = seasons$whit, aes_string(t, ziter2), color = "black", size = 0.8) + # show in front
            geom_line(aes_string(color = "iters"), size = 1)
    }

    # geom_vline(data = pdat2, aes_string(xintercept=date, linetype = pmeth, color = pmeth), size = 0.4, alpha = 1) +
    # scale_linetype_manual(values=c(2, 3, 1, 1, 1, 4))
    # p + facet_grid(meth~pmeth)

    # if (plotly){
    #     plotly::ggplotly(p)
    # }else{

    if (show.legend){
        p <- p + theme(legend.position="none")
        p <- arrangeGrob(p, lgd, nrow = 2, heights = c(min(5*nmethod, 15), 1),
            padding = unit(1, "line"))
    }
    return(p)
}

# make_legend
make_legend <- function(linename = c("iter1", "iter2", "whit"),
        linecolor = c("blue", "red", "black")){
    npoints   <- length(qc_levels)

    labels <- c(qc_levels,linename)
    colors <- c(qc_colors, linecolor)

    # labels <- c(" good", " margin", " snow/ice", " cloud", linename)
    # colors <- c("grey60", "#00BFC4", "#F8766D", "#C77CFF", linecolor)
    nline <- length(linename)
    pch <- c(qc_shapes, rep(NA, nline))

    lty <- rep(0, npoints);  lty[3] <- 1
    lty <- c(lty, rep(1, nline))
    lwd <- c(rep(1, npoints), rep(3, nline))

    I   <- 1:length(colors)
    lgd <- grid::legendGrob(labels[I], pch = pch[I], nrow = 1,
                       # do.lines = TRUE,
                       gp=grid::gpar(lty = lty[I], lwd = lwd[I],
                               cex = 1.4,
                               col = colors[I], fill = colors[I]))
    lgd$children[[5]]$children[[1]]$children %<>% .[2] # fix cross point type
    return(lgd)
}
lgd <- make_legend()
