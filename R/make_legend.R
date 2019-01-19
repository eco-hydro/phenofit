# make_legend
make_legend <- function(linename = c("iter1", "iter2", "whit"),
        linecolor = c("blue", "red", "black"), cex = 1.2){
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
                               cex = cex,
                               col = colors[I], fill = colors[I]))
    lgd$children[[5]]$children[[1]]$children %<>% .[2] # fix cross point type
    return(lgd)
}
lgd <- make_legend()
