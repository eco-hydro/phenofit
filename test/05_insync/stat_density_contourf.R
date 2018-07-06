# main functions
# source("test/Figs/stat_density_contourf.R")
source("test/Figs/stat_prop.R")
fontsize = 14

brks <- c(0.8, 0.5, 0.2)
linewidth <- 1.1
nbrk <- length(brks)
colors <- scales::hue_pal()(nbrk)
colors <- set_names(colors, brks)

# constrained d region manual
fix_value <- function(x, value = 100){
    x[x > value] <- NA#value
    return(x)
}
# d_all$GPP_avg %<>% fix_value(100)
# d_all$EVI_avg %<>% fix_value(100)

p_density_naked <- function(d, title = "(a)", ...){
    linewidth = 1.1
    p_main <- ggplot(d, aes(GPP_avg, EVI_avg)) +
        geom_point(alpha = 0.4)  + #aes(color = index),
        geom_abline(slope = 1, col = "red", size = linewidth) +
        # geom_smooth(method = "lm", formula = y ~ x, color = "black") +
        # geom_contour(data = dc, aes(x, y, z = prob, color=..level..), size = 1.2, breaks = prob) +
        # stat_prob_2d(aes(color=..level..), breaks = prob) +
        # scale_color_continuous(breaks = prob) +
        stat_prob_2d(aes(color=as.factor(..level..)), breaks = brks, size = linewidth, show.legend = F) +
        # scale_color_discrete(labels = sprintf("%d%%", brks*100)) +
        # scale_color_continuous(limits=c(0, 0.8), breaks=prob) +
        # lims(x = c(-100, 100), y = c(-100, 100)) +
        # facet_wrap(~gof, scales = "free") +
        # geom_rug(alpha = 0.3, size = 0.5) +
        labs(color = "density", title = title) +
        scale_x_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
        scale_y_continuous(sec.axis = dup_axis(labels = NULL, name = NULL)) +
        theme_bw(base_size = fontsize) +
        theme(legend.position="none",
            plot.margin=unit(c(0,0,0,0),"points"),
            # panel.grid.major = element_blank(),
            # panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            plot.title = element_text(vjust = -4)
            # axis.ticks.length=unit(-0.2, "cm"),
            # axis.ticks.margin=unit(0.4, "cm")
        )

    if (length(grep("Bias", title)) > 0){
        p_main <- p_main +
            geom_vline(xintercept = 0, col = "black", linetype = 2, size = 0.6) +
            geom_hline(yintercept = 0, col = "black", linetype = 2, size = 0.6)
    }

    ggMarginal(p_main, type = "density")
}

g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

# https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
titleGrob <- function(label, margin=unit(c(b=2, l=1, t=0, r=0), "line"), ..., debug=FALSE){
    library(gtable)
    tg <- textGrob(label, ...)
    w <- grobWidth(tg)
    h <- grobHeight(tg)

    g <- gtable("title",
                widths = unit.c(margin[2], w, margin[4]),
                heights = unit.c(margin[3], h, margin[1]), respect = FALSE)
    if(debug){
        rg <- rectGrob()
        pos <- expand.grid(x=1:3, y=1:3)[-5,]
        g <- gtable_add_grob(g, list(rg, rg, rg, rg, rg, rg, rg, rg), t=pos$y, l=pos$x)
    }
    gtable_add_grob(g, tg, t=2, l = 2)
}
################################################################################

## 1. construct a legend
p_lgd <- ggplot(data.frame(x = 1:100, y = 1:100), aes(x, y)) +
    # geom_point(alpha = 0.4)  +
    # geom_abline(slope = 1, col = "red", size = linewidth) +
    # geom_smooth(method = "lm", formula = y ~ x, color = "black") +
    stat_prob_2d(aes(color=as.factor(..level..)), breaks = brks, size = linewidth, show.legend = T) +
    labs(color = "density") +
    scale_color_discrete(labels = sprintf("%d%%", c(0.2, 0.5, 0.8)*100)) +
    theme_bw(base_size = 16) +
    theme(
        plot.title = element_text(size = fontsize, face = "bold"),
        legend.position="bottom",
        legend.text = element_text(size = fontsize),
        legend.key.width = unit(2,"cm")) +
    guides(colour = guide_legend(NULL, nrow = 1))
lgd <- g_legend(p_lgd)
