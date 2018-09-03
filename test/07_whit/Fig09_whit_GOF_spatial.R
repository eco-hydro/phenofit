library(lattice)
library(latticeExtra)
library(grid)
library(gridExtra)
library(rgdal)
windowsFonts(Times = windowsFont("Times New Roman"),
             Arial = windowsFont("Arial"))

poly <- readOGR("F:/ArcGIS/continent.shp")

d <- merge(df[meth == "wWH" & iter == "iter1", .(site, R2, Bias, RMSE, Roughtness = Rg)],
           st[, .(site, lon, lat)]) %>% .[, 2:ncol(.)]

sp <- df2sp(d)

################################################################################
title = eval(substitute(expression(bold("("*lab*") "* text)),
                                  list(lab = letters[1], text = quote(R^2))))
titles <- list(title, "(b) Bias", "(c) RMSE", "(d) Roughtness")




plot_sppoint <- function(sp, i, brks, sp.layout, cols){
    brks   <- c(-Inf, brks, Inf)
    ncolor <- length(brks) - 1

    fmt <- ifelse(diff(brks[2:3]) < 0.01, "%.3f", "%.2f")
    labels <- sprintf(fmt, brks[-c(1, ncolor+1)])
    at <- 1:(ncolor-1) + 0.5

    data <- sp@data[, i]
    sp@data <- data.table(x = cut(data, brks, include.lowest = T) )

    # browser()
    spplot(sp, as.table = T, cex = 0.4,
           # scales = list(at = brks),
           sp.layout = sp.layout,
           col.regions = cols,
           # edge.col = alpha("grey70", 0.9),
           xlim = c(-182, 182), ylim = c(-60, 87),
           # cuts = 1:ncolor,
           colorkey = list(
               right = list( # see ?levelplot in package trellis, argument colorkey:
                   fun = draw.colorkey,
                   args = list(
                       key = list(
                           at = 0:ncolor+0.5, # colour breaks
                           labels = list(labels = labels, at = at), #at = brks
                           # tck = 1,
                           col = cols, # colours
                           height = 0.95, width = 1.3,
                           axis.line = list(col = 'black')
                       )
                   )
               )
           ),
           par.settings = theme.nopadding)
}

{
    lst_brks <- list(
        seq(0.2, 0.8, 0.2),
        c(seq(0, .05, .01)),
        seq(0, 0.1, 0.02),
        seq(0.01, 0.03, 0.005)
    )

    sp_poly <- list("sp.polygons", poly, col = "grey60", fill = "grey85")
    sp_text <- list("sp.text", -170, 80, title,
                    fontfamily = "Times", cex = 1.3,
                    font = 2, fontface = "bold", adj = 0)
    sp.layout <- c(sp_poly, sp_text)

    ps <- list()
    for (i in 1:4){
        # print(i)
        brks <- lst_brks[[i]]
        ncolor <- length(brks) + 1
        sp_text <- list("sp.text", c(-170, 80), titles[[i]],
                        fontfamily = "Times", cex = 1.3,
                        font = 2, adj = 0)
        sp.layout <- list(sp_poly, sp_text)

        print(brks)
        cols <- bpy.colors(11)
        if (i == 1) {
            red <- RColorBrewer::brewer.pal(11, "Spectral")[c(2)]
            cols <- rev(RColorBrewer::brewer.pal(ncolor, "Spectral"))
            cols[ncolor] <- red
        } else {
            cols <- RColorBrewer::brewer.pal(ncolor, "BrBG") #YlOrRd
        }
        ps[[i]] <- plot_sppoint(sp, i, brks, sp.layout, cols)
    }

    file <- "Fig9_wWH_spatial_performance.jpg"
    # tiff(file, 15, 6, units = "in", res = 300, compression = "lzw")
    CairoPNG(file, 15, 6, units = "in", dpi = 300)
    p <- arrangeGrob(grobs = ps, nrow = 2, ncol = 2, padding = unit(0, "line"))
    grid.newpage()
    grid.draw(p)
    dev.off()
    file.show(file)
}


theme.nopadding <-
    list(axis.line = list(col = "transparent"),
         layout.heights =
             list(top.padding = 0,
                  main.key.padding = 0,
                  key.axis.padding = 0,
                  axis.xlab.padding = 0,
                  xlab.key.padding = 0,
                  key.sub.padding = 0,
                  bottom.padding = 0),
         layout.widths =
             list(left.padding = 0,
                  key.ylab.padding = 0,
                  ylab.axis.padding = 0,
                  axis.key.padding = 0,
                  right.padding = 0))
