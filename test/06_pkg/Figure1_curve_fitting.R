library(lubridate)

ratio    <- 1.8
lambda   <- 2
sites_sm <- sites[c(1, 3, 4, 6, 9)]
df_sm    <- df[site %in% sites_sm, ]

lambda = 2
ps   <- list()
fits <- list()

for (i in seq_along(sites)){
    runningId(i)
    sitename <- sites[i]

    d <- df[site == sitename , ]
    # Check input data and initial parameters for phenofit
    INPUT <- check_input(d$date, d$y, d$w, trim = T, maxgap = nptperyear / 4, alpha = 0.02)
    # The detailed information of those parameters can be seen in `season`.
    brks <- season(INPUT, lambda, nptperyear, iters = 3, wFUN = wFUN, IsPlot = F,
                   south = d$lat[1] < 0,
                   Aymin_less = 0.6, ymax_min = ymax_min,
                   max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5) #, ...
    fit  <- curvefits(INPUT, brks, lambda =lambda, IsPlot = T,
                      methods = c("AG", "zhang", "beck", "elmore", 'Gu')[c(4)], #,"klos"
                      nptperyear = nptperyear, debug = F,
                      wFUN = wFUN,
                      qc = d$SummaryQA,
                      extend_month = 2,
                      south = d$lat[1] < 0)
    fit$INPUT   <- INPUT
    fit$seasons <- brks
    fits[[i]]   <- fit

    titlename <- sprintf('(%s) %s, %s', letters[i], d$site[1], d$IGBP[1])
    p <- plot_phenofit(fit, d, titlename, show.legend = F) +
        scale_x_date(breaks = ymd(paste0(2001:2016, "0101")),
                       labels = 2001:2016) + ylab("EVI") +
        theme_light() +
        theme(legend.position="none",
              plot.title = element_text(vjust = -6, hjust = 0.01,
                                        margin = margin(1, 1, 1, 1, "cm")*0),
              plot.margin = unit(c(0,1,1,1)*0.3, "cm"),
              panel.grid.minor = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.major.x = element_line(size = 0.5),
              strip.text = element_blank()
              ) +
        coord_cartesian(xlim = c(ymd(20000101), ymd(20161231)))
    if (i < length(sites)){
        p <- p + theme(
            axis.title.x = element_blank(),
            axis.text.x  = element_blank()
        )
    }
    # temp <- ExtractPheno(fit$fits$ELMORE, IsPlot = T)
    ps[[i]] <- p
}

names(fits) <- sites

n <- length(sites)
library(gridExtra)
file <- "Figure4_phenofit_curvefitting.pdf"
cairo_pdf(file, 10, 20)
fig <- do.call(arrangeGrob, c(ps, list(ncol = 1)))
grid.arrange(fig, lgd, ncol = 1, heights = c(n * 5, 1))

dev.off()
file.show(file)

ratio = 1.15
cairo_pdf("Figure5_Phenology_Extraction.pdf", 8*ratio, 6*ratio)
temp <- ExtractPheno(fit$fits$ELMORE[1:6], IsPlot = T, TRS = 0.5)
dev.off()
file.show("Figure5_Phenology_Extraction.pdf")
# grid.arrange(n)
# Cairo::CairoPDF("phenofit_5st_v4.pdf", 10*ratio, 6*ratio)
# fits <- dlply(df_sm, .(site), function(d){
#     fit <- NA
#     tryCatch({
#
#         print(p)
#     }, error = function(e){
#         message(sprintf("[e] %s:%s", e$message, d$site[1]))
#     })
#     return(fit)
# })
# dev.off()

p + theme_light() +
    theme(
            # legend.position="top",
          plot.title = element_text(vjust = -6, hjust = 0.01,
                                    margin = margin(1, 1, 1, 1, "cm")*0),
          plot.margin = unit(c(0,1,1,1)*0.3, "cm"),
          panel.grid.minor = element_blank()
    ) +
#
#     scale_shape(guide = FALSE) +
#     scale_linetype(guide = TRUE) +
    guides(color = guide_legend(title = "lgd"),
           shape = guide_legend(title = "lgd"),
           fill  = guide_legend("lgd", nrow = 1, order = 1)
           ) +
    labs(shape="Male/Female", colour="Male/Female")

grid.newpage()

grid.draw(lgd1)

grid.newpage()
I <- 5:7
lgd2 <- legendGrob(labels[I],
                   pch = c(NA, NA, NA),  nrow = 1,
                   do.lines = T,
                   gp=gpar(lty = c( rep(1, 3)), lwd = c(2, 2, 2),
                           col = colors[I], fill = colors[I]))
grid.draw(lgd2)
p <-  grid.arrange(lgd1, lgd2, nrow = 2, padding = unit(0, "line"))

