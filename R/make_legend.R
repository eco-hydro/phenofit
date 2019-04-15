# make_legend
make_legend <- function(
    linename = c("iter1", "iter2", "whit"),
    linecolor = c("blue", "red", "black"), cex = 1.2,
    nmax_points = 6, nrow = 1)
{
    nline <- length(linename)
    nmax_points <- min(nmax_points, length(qc_levels))

    I_sel <-  1:nmax_points
    labels <- c(qc_levels[I_sel], linename)
    colors <- c(qc_colors[I_sel], linecolor)[1:(nline+nmax_points)]

    # labels <- c(" good", " margin", " snow/ice", " cloud", linename)
    # colors <- c("grey60", "#00BFC4", "#F8766D", "#C77CFF", linecolor)
    pch <- c(qc_shapes[I_sel], rep(NA, nline))

    lty <- rep(0, nmax_points);  lty[3] <- 1
    lty <- c(lty, rep(1, nline))
    lwd <- c(rep(1, nmax_points), rep(3, nline))

    I   <- 1:length(colors)
    lgd <- grid::legendGrob(labels[I], pch = pch[I],
        nrow = nrow, byrow = T,
        # do.lines = TRUE,
        gp=grid::gpar(lty = lty[I], lwd = lwd[I],
               cex = cex,
               col = colors[I], fill = colors[I]))
    lgd$children[[5]]$children[[1]]$children %<>% .[2] # fix cross point type (3th)
    return(lgd)
}


# legend of MOD13A1 or MOD09A1
# make sure QC_flag is factor
make_legend_nmax <- function(linename, linecolor, QC_flag, ...){
    nmax_points <- 4
    if (!is.null(QC_flag)) {
        nmax_points <- factor(QC_flag) %>% as.numeric() %>% max(na.rm = T)
        nmax_points <- ifelse(nmax_points <= 4, 4, 6)
    } 

    make_legend(linename , linecolor, nmax_points = nmax_points, ...)
}


lgd <- make_legend()
lgd_short <- make_legend(nmax_points = 4)
