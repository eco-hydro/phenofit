source('test/stable/load_pkgs.R')
source("test/07_whit/main_phenofit_test.R")
source("test/07_whit/dat_flux&cam_phenofit.R")
source('test/stable/ggplot/geom_boxplot_no_outlier.R')

load("data_test/whit_lambda/MOD13A1_st_1e3_20180731.rda")
methods <- c('AG', 'BECK', 'ELMORE', 'ZHANG', 'whit_R', 'wHANTS', 'wSG', 'wWH', "wWH2")
methods2 <- c('AG', 'ZHANG', 'wHANTS', 'wSG', 'wWH', "wWH2")

quantile_envelope <- function(x, alpha){
    names <- "ymean"
    # if (alpha != 0.5) {
        alpha <-  c(alpha, 1-alpha)
        names <- c("ymin", "ymax")
    # }
    res <- quantile(x, alpha, na.rm = T)
    set_names(res, names)
}

fix_name <- function(x) {
    names(x) %<>% str_extract(".*(?=_)")
    # names(x)
    x
}

tidy_info <- function(file){
    x <- readRDS(file)
    a <- llply(x, fix_name) %>% purrr::transpose() %>%
        llply(function(l) {
            info <- map(l, "info") %>% do.call(rbind, .)
            rough <- map(l, "rough") %>% do.call(rbind, .)
            list(info = info, rough = rough)
        })
    names <- names(a)

    for (i in 1:length(a)){
        name <- names[i]
        d_info  <- a[[i]]$info
        d_rough <- a[[i]]$rough

        if (name != "phenofit"){
            d_info$meth <- name
            d_rough$meth <- name
        }
        d_info  %<>% reorder_name(c("site", "meth"))
        d_rough %<>% reorder_name(c("site", "meth"))

        a[[i]] <- merge(d_info, d_rough) %>%
            reorder_name(c("site", "meth", "type", "iter"))#list(info = d_info, rough = d_rough)
    }
    a %<>% do.call(rbind, .) #transpose() %>% map(~
    # a <- fix_name(a)
    # d <- .[NSE > 0, ] %>%
    #     melt(id.vars = c("site", "meth", "type"), variable.name = "index") %>%
    #     .[index %in% c("R2", "NSE", "RMSE")]

    # d$meth %<>% factor(methods)
    a$Rg %<>% unlist()
    return(a)
}

###############################################################################
indir <- "V:/result/val_info/info_0/"
files <- dir(indir, full.names = T, recursive = T)

df <- llply(files, tidy_info, .progress = "text") %>%
    # set_names(c("p10%", "p30%", "p50%")) %>%
    # melt_list("perc")
    do.call(rbind, .)
df$meth %<>% factor(methods)

# df <- melt_list(lst, "meth")
# df$Rg %<>% unlist()

st$site %<>% as.character()
st$IGBPname %<>% factor(IGBPnames_006)

indice <- c("R^2", "Bias", "RMSE", "Roughness")
labels <- c("bold((a)*' '*R^2)", "bold((b)*' '*Bias)",
            "bold((c)*' '*RMSE)", "bold((d)*' '*Roughness)")

d <- df[, .(site, meth, type, iter, RMSE, `R^2` = R2, Bias, Roughness = Rg)] %>%
    melt(id.vars = c("site", "meth", "type", "iter"), variable.name = "index") #, "perc"
d <- merge(st[, .(site, IGBPname)], d)

d$index %<>% factor(indice)# <- indice# factor(indices)
d$index %<>% mapvalues(indice, labels)

d1 <- d[meth %in% c("wWH", "wWH2") & iter == "iter1"] %>%
    dcast(site+IGBPname+index~meth, value.var = "value")
# wWH vs wWH2 -------------------------------------------------------------

dmin <- 0.1
d1[, `:=`(diff = wWH - wWH2, kind = 0)]
d1[diff >  dmin, kind := 1]
d1[diff < -dmin, kind := -1]
table(d1[index == "R^2", ]$kind)

source("test/stable/ggplot/stat_prop.R")
brks <- c(0.8, 0.5, 0.2)

d2 <- melt(d1, c("site", "IGBPname", "index"), variable.name = "meth")
ggplot(d2, aes(meth, value)) +
    geom_boxplot2() +
    facet_wrap(~index, scale = "free", labeller = label_parsed)
p1 <- ggplot(d1, aes(wWH2, wWH)) +
    geom_point(alpha = 0.1) +
    geom_abline(slope = 1, color = "red", size = 0.8) +
    # stat_prob_2d(aes(fill=as.factor(..level..), alpha=1 - ..level..), geom = "polygon",
    #              breaks = brks, size = 1, show.legend = T) +
    # stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black')+
    facet_wrap(~index, scale = "free", labeller = label_parsed) +
    theme_grey(base_size = 12) + #, base_family = "Arial"
    theme(strip.text = element_text(color = "black", face="bold", size = 12,
                                    margin = margin(1, 0, 1, 0, "pt")),
          axis.text = element_text(size = 11),
          axis.title = element_text(size = 12, face = "bold"))
    # coord_equal()
    # geom_hex(bins = 30) +
    # scale_fill_continuous(low="green",high="red")


file = "Fig8_wWH_vs_wWH2.pdf"
write_fig(p1, file, 9, 7, T)
write_fig(p1, gsub(".pdf", ".tif", file), 9, 7, T)

file.show(file)



# vs other methods --------------------------------------------------------
d <- df[, .(site, meth, type, iter, RMSE, `R^2` = R2, Bias, Roughness = Rg)] %>%
    melt(id.vars = c("site", "meth", "type", "iter"), variable.name = "index") #, "perc"
d <- merge(st[, .(site, IGBPname)], d)

d$index %<>% factor(indice)# <- indice# factor(indices)
# d$index %<>% mapvalues(indice, labels)

pdat <- d[meth %in% methods2 & iter == "iter1"]
# pdat$IGBPname %<>% as.numeric()
# pdat$facet <- factor(p$index, levels = labels)
d_lab <- data.table(index = factor(indice, indice),
                    label = sprintf("(%s)", letters[1:length(indice)]))
grid_x <- geom_vline(xintercept = (1:(length(IGBPnames_006) - 3)) + 0.5,
                     linetype = 2, size = 0.3, color = "grey70")

colors <- scales::hue_pal()(length(methods2))
colors[length(colors) - 1] <- "black"
# important lesson, numeric and factor can't be mixed together
p <- ggplot(pdat, aes(IGBPname, value, color = meth)) +
    # geom_boxplot() +
    facet_grid(index~., scale = "free", labeller=label_parsed, switch = "y") + #, nrow = 5, strip.position = "left"
    # geom_boxplot() +
    # scale_x_discrete(limits = IGBPnames_006[1:16]) +
    # scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
    labs(x = "IGBP", y = NULL) +
    scale_color_manual(values = colors) +
    guides(color = guide_legend(nrow = 1)) + #, keywidth=0.6
    theme_light(base_size = fontsize) + #, base_family = "Arial"
    theme(legend.position = "bottom",
          legend.spacing.x = unit(0.1, 'cm'),
          legend.margin = margin(-0.4, unit = "cm"),
          legend.title = element_blank(),
          plot.margin = margin(0.2, 0.1, 0.2, 0, unit = "cm"),
          # panel.spacing = unit(0.2, "lines"),
          # panel.grid.major = element_line(linetype = 2),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(size = 0.3, linetype = 2),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(color = "black", face="bold", size = 14,
                                    margin = margin(1, 0, 1, 0, "pt")),
          strip.placement = "outside",
          axis.text = element_text(color = "black"))

boxplot2 <- function(p, width = 0.95, size = 0.7){
    # width  <- 0.95
    width2 <- width - 0.15
    dodge <- position_dodge(width = width)

    p +
        stat_summary(fun.data = box_qtl,
                     position = dodge, size = size,
                     geom = "errorbar", width = width2) +
        geom_boxplot2(coef = 0,
                      width = width2,
                      lwd = size - 0.2,
                      notch = F, outlier.shape = NA, position=dodge) +
        grid_x +
        geom_text(data = d_lab, aes(x = "ENF",
                                    y = Inf, color = NULL, label = label),
                  vjust = 1.5, hjust = 1.1, fontface = "bold", size =5, show.legend = F)

}

p2 <- boxplot2(p, width = 0.88, size = 0.62); #p2

file <- "Fig8_compare_with_other_methods_iter1.pdf"
file_tiff <-  gsub(".pdf", ".tif", file)
write_fig(p2, file_tiff, 11, 9, T, res = 200)

write_fig(p2, file, 11, 9, T)

# write_fig(p1, , 9, 7, T)
# write_tiff(p2, gsub(".pdf", ".tiff", file), 11, 9)



colors <- scales::hue_pal()(3)
colors <- c(colors[2], "grey50", colors[1])

ggplot(d2[meth == "wWH2"], aes(wWH, value)) +
    geom_point(aes(color = kind), alpha = 0.2) +
    # facet_wrap(~label, scales = "free") +
    facet_wrap(~meth+index, labeller = label_value, scales = "free", ncol = 2) +
    geom_abline(slope = 1, color ="red", size = 0.1) +
    scale_color_manual(values = colors)
    # geom_text(data = d_dominant, aes(label = Bigger), fontface = "bold",
    #           x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5, show.legend = F, color = colors[1]) +
    # geom_text(data = d_dominant, label = " vs ", fontface = "bold",
    #           x = -Inf, y = Inf, hjust = -0.6, vjust = 1.5, show.legend = F) +
    # geom_text(data = d_dominant, aes(label = Smaller), fontface = "bold",
    #           x = -Inf, y = Inf, hjust = -1.3, vjust = 1.5, show.legend = F, color = colors[3])



###############################################################################
## 2. acf

boxplot <- function(p, width = 0.95, size = 0.7){
    # width  <- 0.95
    width2 <- width - 0.15
    dodge <- position_dodge(width = width)

    p + stat_summary(fun.data = box_qtl,
                     position = dodge, size = size,
                     geom = "errorbar", width = width2) +
        geom_boxplot2(coef = 0,
                      width = width2,
                      lwd = size,
                      notch = F, outlier.shape = NA, position=dodge) +
        theme_light(base_size = fontsize, base_family = "Arial") +
        theme(legend.position = c(1-0.01, 0.01), legend.justification = c(1, 0),
              panel.grid.major = element_line(linetype = 2),
              panel.grid.minor = element_blank(),
              legend.title=element_blank(),
              axis.text = element_text(color = "black"))
}

d_acf <- df[, .(site, meth, type, iter, acf)] %>%
    melt(id.vars = c("site", "meth", "type", "iter"), variable.name = "index")

acf <- df$acf %>% do.call(rbind, .) %>%
    set_colnames(paste0("tag", 1:10)) %>% data.table()
d_acf <- cbind(df[, .(site, meth, type, iter)], acf) %>%
    melt(c("site", "meth", "type", "iter"), variable.name = "tag")
d_acf$tag %<>% mapvalues(levels(.), gsub("tag", "", levels(.)))

d_acf <- d_acf[meth %in% methods2, ]
# d[, median(value, na.rm = T), .(meth, type, index)] %>%
#     dcast(meth+type~index, value.var = "V1") %>%
#     .[order(type)]

p <- ggplot(d_acf[meth != "whit_R"], aes(tag, value, color = meth)) +
    labs(x = "lag (images)", y = "ACF of residual")
    # geom_boxplot2()
p2 <- boxplot(p, width = 0.8) +
    guides(color = guide_legend(nrow = 1)) +
    scale_color_manual(values = colors) +
    theme(
        legend.position = c(1.04, 0.98)*0.95, legend.justification = c(1, 1),
        legend.margin = margin(-0.4, unit = "cm"),
        legend.spacing.x = unit(0.1, 'cm')
        # legend.spacing.x = unit(0.1, 'cm')
    )

file <- "Fig9_autocorrelation.pdf"
file_tiff <-  gsub(".pdf", ".tif", file)
write_fig(p2, file_tiff, 9, 5, T)

write_fig(p2, file, 10, 5.5, T)


###############################################################################
## 03. Fig10. The influence of good values percentage
levels <- seq(0, 1, 0.01)
xmid   <- c((levels[-1] + levels[-length(levels)])/2, 1)

df[, grp_perc := findInterval(perc_good, levels, rightmost.closed = T)]
df[, xmid := xmid[grp_perc]]

d <- df[, .(site, meth, type, iter, RMSE, R2, Bias, Rg, grp_perc, xmid)] %>% #, "NSE"
    melt(id.vars = c("site", "meth", "type", "iter", "grp_perc", "xmid"), variable.name = "index")

ggplot(d[index == "R2"], aes(xmid,value)) + geom_point() + geom_density2d() +
    facet_wrap(~meth)
# d <- df[meth == "wWH" & iter == "iter1"]

alphas <- c(.05, .1, .25, .5) %>% set_names(., .)
d_envelope <- llply(alphas, function(alpha){
    # print(alpha)
    expr <- substitute(quote(res <- ddply_dt(d, .(quantile_envelope(value, alpha)), .(meth, index, iter, grp_perc))),
               list(alpha = alpha))
    # print(eval(expr))
    eval(eval(expr))
})

d_enve <- d_envelope %>% melt_list("alpha")
d_enve[, xmid := xmid[grp_perc]]

# hue_pal()(4) %>% show_col()
cols <- hue_pal()(4)
colors <- c(trans_col(cols[1], 0.5), cols[2:3], "white")

pdat <- d_enve[meth == "wWH" & iter == "iter1"]

indice <- c("R2", "Bias", "RMSE", "Rg")
labels <- c("bold((a)*' '*R^2)", "bold((b)*' '*Bias)",
            "bold((c)*' '*RMSE)", "bold((d)*' '*Roughness)")
# pdat <- d[meth != "whit_R" & iter == "iter2" & index %in% indice]
pdat$index %<>% factor(indice, labels)

# sub_id <- geom_text(data = data.table(index = factor(labels, labels),
#                                        label = sprintf("(%s)", letters[1:length(labels)])),
#                     aes(label = label),
#                     x = -Inf, y = Inf, vjust = 1, hjust = 0)
d_vline <- data.table(index = factor(labels[1:3], labels),
                      x = c(.3, .25, .25))
v_line <- geom_vline(data = d_vline ,
                     aes(xintercept = x), color = "black", linetype = 2, size = 1)
v_line_lab <- geom_text(data = d_vline, hjust = -0.1, vjust = -0.4,
                        aes(x, y = -Inf,label = sprintf("%2d%%", x*100)))

p <- ggplot(pdat, aes(xmid, ymin)) +
    # geom_point() + geom_density2d() +
    geom_ribbon(aes(x = xmid, ymin = ymin, ymax = ymax, fill = alpha)) +
    facet_wrap(~index, scales = "free", labeller=label_parsed) +
    # sub_id +
    geom_line(data = pdat[alpha == 0.5, ], color ="white") +
    v_line + v_line_lab +
    # geom_vline(xintercept = c(0.25), color= "blue", size = 1, linetype = 2) +
    scale_fill_manual(values = colors, labels = c("  5% and 95%", "10% and 90%", "25% and 75%", "50%")) +
    guides(fill = guide_legend(override.aes = list(color = "black"))) +
    labs(x = "The percentage of good values (%)") +
    scale_x_continuous(labels = function(x) sprintf("%d", x*100)) +
    scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
    theme_gray(base_size = 14) +
    theme(
        legend.title = element_blank(),
        legend.position = c(0.98, 0.985),
        legend.justification = c(1, 1),
        legend.spacing.x = unit(0.1, 'cm'),
        strip.text = element_text(face = "bold", size = 13, margin = margin(2, 0, 1, 0, "pt")),  # not work for meth expr
        axis.text = element_text(size = 12),
        axis.title.x = element_text(face = "bold", size = 13)
    )+
    ylab(NULL)

file <- "Fig10_good_values_percentage_impact.pdf"
cairo_pdf(file, 9.5, 6)
print(p)
dev.off(); file.show(file)

file <- gsub(".pdf", ".tiff", file)
write_tiff(p, file, 9.5, 6)

## 4. GOF spatial distribution

d <- merge(df[meth == "wWH", .(site, meth, type, iter, RMSE, R2, Bias, Roughtness = Rg)],
           st[, .(site, lon, lat)]) %>%
    melt(c("site", "meth", "type", "iter", "lon", "lat"), variable.name = "index")
d
ggplot(d, aes(lon, lat)) + geom_point(aes(color = R2))

indice <- c("R2", "Bias", "RMSE", "Roughtness")
ps <- list()
for (i in 1:4){
    pdat <- d[index == indice[i] & iter == "iter1", ]
    p <- ggplot(pdat) +
        geom_polygon(data = d_poly, aes(long, lat, group = group), fill = "grey85", colour = "black") +
        coord_fixed(xlim = c(-180, 180), ylim = c(-55, 85), ratio = 1) +
        geom_point(aes(lon, lat, color = value), size = 0.8) +
        scale_colour_gradientn(colors = colors) +
        theme_void() +
        guides(color = guide_colorbar(barheight  = 6)) +
        theme(plot.margin = margin(-50, -10, -50, -10, "pt"),
              axis.text = element_blank(),
              axis.title = element_blank()
              )
    ps[[i]] <- p
}
ps[[1]]

tiff("a.png",10, 15, units = "in", res = 300, dpi = 300)
CairoPNG()
p <- arrangeGrob(grobs = ps, nrow = 4, ncol = 1)
grid.newpage()
grid.draw(p)
dev.off()
# file.show("a.png")

library(rgdal)
d_poly <- readOGR("F:/ArcGIS/continent.shp") %>%
    tidy(.[, "CONTINENT"], region = "CONTINENT")

colors <- rev(RColorBrewer::brewer.pal(11, "Spectral"))


    # coord_map()
