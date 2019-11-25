source("test/load_pkgs.R")
source("test/07_whit/main_phenofit_test.R")
source("test/07_whit/main_Whittaker_Figures.R")
source("test/07_whit/dat_flux&cam_phenofit.R")
source('test/stable/ggplot/geom_boxplot_no_outlier.R')

load("data_test/whit_lambda/MOD13A1_st_1e3_20180731.rda")
methods <- c('AG', 'BECK', 'ELMORE', 'ZHANG', 'whit_R', 'smooth_wHANTS', 'smooth_wSG', 'wWH', "wWH2")
methods2 <- c('AG', 'ZHANG', 'smooth_wHANTS', 'smooth_wSG', 'wWH', "wWH2")

st$site %<>% as.character()
st$IGBPname %<>% factor(IGBPnames_006)

indice       <- c("R2", "Bias", "RMSE", "Rg")
indice_label <- c("bold((a)*' '*R^2)", "bold((b)*' '*Bias)",
                  "bold((c)*' '*RMSE)", "bold((d)*' '*Roughness)")

###############################################################################
indir <- "V:/result/val_info/info_0_v2/"
files <- dir(indir, full.names = T, recursive = T)

df <- llply(files, tidy_info, .progress = "text") %>%
    # set_names(c("p10%", "p30%", "p50%")) %>%
    # melt_list("perc")
    do.call(rbind, .)


# df$meth %<>% factor(methods2)

# df <- melt_list(lst, "meth")
# df$Rg %<>% unlist()

d <- df[, .(site, meth, type, iters, RMSE, R2 = R2, Bias, Rg, Rg_norm_by_obs, Rg_norm_by_pred)] %>%
    melt(id.vars = c("site", "meth", "type", "iters"), variable.name = "index") #, "perc"
d <- merge(st[, .(site, IGBPname)], d)
d$index %<>% factor(indice, indice_label)# <- indice# factor(indices)

d1 <- d[meth %in% c("wWH", "wWH2") & iters == "iter1"] %>%
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
d <- df[, .(site, meth, type, iters, RMSE, `R^2` = R2, Bias, Roughness = Rg)] %>%
    melt(id.vars = c("site", "meth", "type", "iters"), variable.name = "index") #, "perc"
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

d_acf <- df[, .(site, meth, type, iters, acf)] %>%
    melt(id.vars = c("site", "meth", "type", "iters"), variable.name = "index")

acf <- df$acf %>% do.call(rbind, .) %>%
    set_colnames(paste0("tag", 1:10)) %>% data.table()
d_acf <- cbind(df[, .(site, meth, type, iters)], acf) %>%
    melt(c("site", "meth", "type", "iters"), variable.name = "tag")
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



## 4. GOF spatial distribution
