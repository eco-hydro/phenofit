#' ggplot_dual_axis
#'
#' Takes 2 ggplot plots and makes a dual y-axis plot
#' function takes 2 compulsory arguments and 1 optional argument
#'
#' @param lhs is the ggplot whose y-axis is to be displayed on the left
#' @param rhs is the ggplot whose y-axis is to be displayed on the right
#' @param axis.title.y.rhs takes value "rotate" to rotate right y-axis label
#'
#' The function does as little as possible, namely:
#' display the lhs plot without minor grid lines and with a
#' transparent background to allow grid lines to show
#' display the rhs plot without minor grid lines and with a secondary y axis,
#' a rotated axis label, without minor grid lines.
#' Justify the y-axis label by setting 'hjust = 0' in 'axis.text.y'
#' rotate the right plot 'axis.title.y' by 270 degrees, for symmetry
#' rotation can be turned off with 'axis.title.y.rhs' option
#'
#' @references
#' [1]. http://stackoverflow.com/questions/18989001/how-can-i-put-a-transformed-scale-on-the-right-side-of-a-ggplot2
#' [2]. https://stackoverflow.com/questions/28186861/ggplot2-when-overlapping-two-plots-to-get-axes-on-the-right-legend-from-second
#'
ggplot_dual_axis <- function(lhs, rhs, g3,  add_yaxis_r = T, axis.title.y.rhs = "rotate") {

    g1 <- lhs
    g2 <- rhs

    theme_bg <- theme(
            panel.grid = element_blank(),
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA)
        )

    if ("ggplot" %in% class(lhs)){
        # 3a. Use only major grid lines for the left axis
        lhs <- lhs + theme(panel.grid.minor = element_blank())
        g1 <- ggplotGrob (lhs)
    }

    if ("ggplot" %in% class(rhs)){
        # 1. Fix the right y-axis label justification
        rhs <- rhs + theme(axis.text.y = element_text(hjust = 0))
        # 2. Rotate the right y-axis label by 270 degrees by default
        if (missing(axis.title.y.rhs) | axis.title.y.rhs %in% c("rotate", "rotated")) {
            rhs <- rhs + theme(axis.title.y = element_text(angle = 270))
        }
        # 3b. Use only major grid lines for the right axis
        #     force transparency of the backgrounds to allow grid lines to show
        rhs <- rhs + theme_bg
        g2 <- ggplotGrob (rhs)
    }
    # Process gtable objects
    # 4. Extract gtable
    # 5. Overlap the panel of the rhs plot on that of the lhs plot
    pp <- c(subset(g1$layout, name == "panel", se = t:r))
    g  <- gtable_add_grob(g1,
                          g2$grobs[[which(g2$layout$name == "panel")]], pp$t, pp$l, pp$b, pp$l)

    if (!missing(g3)){
        g3 %<>% `+`(theme_bg)
        g3 <- ggplotGrob(g3)
        g  <- gtable_add_grob(g,
                          g3$grobs[[which(g3$layout$name == "panel")]], pp$t, pp$l, pp$b, pp$l)
    }

    if (add_yaxis_r){
        # Tweak axis position and labels
        I <- which(g2$layout$name == "axis-l")
        ga <- g2$grobs[[I]]
        ax <- ga$children[["axis"]]  # ga$children[[2]]
        ax$widths <- rev(ax$widths)
        ax$grobs  <- rev(ax$grobs)
        ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
        g <- gtable_add_cols(g, g2$widths[g2$layout[I,]$l], length(g$widths) - 1)
        g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b, name = "yaxis-r2")
        # g <- gtable_add_grob(g, g2$grobs[[7]], pp$t, length(g$widths), pp$b)

        ## add dual ylab
        I    <- which(g2$layout$name == "ylab-l")
        ncol <- length(g$widths)
        g = gtable_add_grob(g, g2$grob[[I]], pp$t, ncol, pp$b, name='ylab-r2')
        g$widths[ncol] <- sum(g2$grobs[[I]]$widths) + unit(2, "pt")
        # grid.draw(g_temp)
    }

    # add legend on top
    if ("guide-box" %in% g1$layout$name) {
        dimGB1 <- c(subset(g1$layout, name == "guide-box", se = t:r))
        g <- gtable_add_grob(g,
                             g1$grobs[[which(g1$layout$name == "guide-box")]],
                             dimGB1$t, dimGB1$l, dimGB1$b, dimGB1$l, z = -Inf)
    }
    # Display plot with arrangeGrob wrapper arrangeGrob(g)
    # g <- arrangeGrob(g)
    grid.newpage()
    grid.draw(g)
    return(g)
}

library("gridExtra")
library("grid")
library("gtable") # loads the grid package
#
# p1 <- ggplot(x,aes(date, whit_gee)) +geom_point() + geom_line() #+
#     # scale_y_continuous(sec.axis = sec_axis(~., name = "GPP")); p1
# p2 <- ggplot(x,aes(date, GPP_NT)) +geom_point() + geom_line(); p2

# ggplot_dual_axis(p1, p2)
