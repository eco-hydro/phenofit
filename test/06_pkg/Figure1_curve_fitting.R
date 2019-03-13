library(lubridate)
library(data.table)
library(magrittr)
library(phenofit)
library(plyr)
library(Cairo)
library(gridExtra)

load("Y:/R/phenofit/data/phenofit_MultipleINPUT_flux212.rda")
stations212 <- fread("C:/Users/kon055/Google Drive/Github/data/phenology/station/flux-212.txt")

# shpfile <- "st_20.shp"
# st <- readShapePoints(shpfile) %>% 
#     { data.table( set_colnames(coordinates(.), c("lon", "lat")), .@data)}

df <- lst$MOD13A1 %>%
    merge(stations212[, .(site, lat, IGBP)], by = "site") %>%
    reorder_name(c("site", "IGBP", "lat", "long"))
df <- df[date >= ymd(20010101), ]

sites   <- unique(df$site)
# rename variable name as y
varname <- "EVI_500m"
I_var   <- which(colnames(df) == varname[1])
colnames(df)[I_var] <- "y"

# remap SummaryQA factor level, plot_phenofit also use this variable
if ('SummaryQA' %in% colnames(df)){
    values <- c("0", "1", "2", "3")
    levels <- c(" good", " margin", " snow&ice", " cloud")
    df$SummaryQA %<>% factor() %>% mapvalues(values, levels)
}
## -----------------------------------------------------------------------------
ratio    <- 1.8
lambda   <- 2
sites_sm <- sites[c(1, 3, 4, 6, 9)]
df_sm    <- df[site %in% sites_sm, ]

lambda <- 2
ps   <- list()
fits <- list()

nptperyear <- 23

# test seasons end --------------------------------------------------------
for (i in seq_along(sites)){
    runningId(i)
    sitename <- sites[i]

    d <- df[site == sitename , ]
    dnew <- add_HeadTail(d)
    # Check input data and initial parameters for phenofit
    INPUT <- check_input(dnew$t, dnew$y, dnew$w, trim = T, maxgap = nptperyear/4, alpha = 0.02)
    INPUT$y0 <- dnew$y

    # 2. The detailed information of those parameters can be seen in `season`.
    lambda <- init_lambda(INPUT$y)#*2
    brks <- season(INPUT, nptperyear,
                   wFUN = wFUN, iters = 3,
                   FUN = whitsmw2, lambda = lambda,
                   IsPlot = F, south = d$lat[1] < 0,
                   rymin_less = 0.6, ymax_min = ymax_min,
                   max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5) #, ...
    titlename <- sprintf('(%d) %s, %s', i, d$site[1], d$IGBP[1])

    brks2 <- whit_brks(dnew, nptperyear = 23, FUN = whitsmw2, IsPlot = T, partial = F)

    # 3. curve fitting
    fit  <- curvefits(INPUT, brks2, lambda =lambda, IsPlot = T,
                      methods = c("AG", "zhang", "beck", "elmore", 'Gu'), #,"klos"
                      nptperyear = nptperyear, debug = F, wFUN = wTSM,
                      ymax_min = ymax_min,
                      extend_month = 2,
                      qc = as.numeric(dnew$SummaryQA),
                      south = d$lat[1] < 0)
    fit$INPUT   <- INPUT
    fit$seasons <- brks2
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
              strip.text = element_blank()) +
        coord_cartesian(xlim = c(ymd(20000101), ymd(20161231)))
    # if (i < length(sites)){
    #     p <- p + theme(
    #         axis.title.x = element_blank(),
    #         axis.text.x  = element_blank()
    #     )
    # }
    # temp <- get_pheno(fit$fits$ELMORE, IsPlot = T)
    ps[[i]] <- p
}
names(fits) <- sites

## check curve fitting
n <- length(sites)
file <- "Figure4_phenofit_curvefitting.pdf"
cairo_pdf(file, 10, 20)
fig <- do.call(arrangeGrob, c(ps, list(ncol = 1)))
grid.arrange(fig, lgd, ncol = 1, heights = c(n * 5, 1))

dev.off()
file.show(file)

ratio = 1.15
cairo_pdf("Figure5_Phenology_Extraction.pdf", 8*ratio, 6*ratio)
temp <- get_pheno(fit$fits$ELMORE[1:6], IsPlot = T, TRS = 0.5)
dev.off()
file.show("Figure5_Phenology_Extraction.pdf")
# grid.arrange(n)
# Cairo::CairoPDF("phenofit_5st_v4.pdf", 10*ratio, 6*ratio)
# fits <- dlply(df_sm, .(site), function(d){
#     fit <- NA
#     tryCatch({
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
          panel.grid.minor = element_blank() ) +
#     scale_shape(guide = FALSE) +
#     scale_linetype(guide = TRUE) +
    guides(color = guide_legend(title = "lgd"),
           shape = guide_legend(title = "lgd"),
           fill  = guide_legend("lgd", nrow = 1, order = 1) ) +
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
