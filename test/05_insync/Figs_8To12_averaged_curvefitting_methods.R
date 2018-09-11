source("test/Figs/stat_density_contourf.R")
source('F:/Github/PML_v2/fluxsites_tidy/R/pkg_vis.R')
library(ggExtra)
library(grid)

# 4 Figures main function: only for averaged curve fitting methods
p_boxplot <- function(type = "SOS", real = "TRS1",
    df_list = df_lst_nometh[c("GPP_mod","GPP_vpm", "MOD13A1_EVI", "MOD13Q1_EVI" )])
{
    title <- paste(type[1], real, sep = '_')

    lst <- llply(df_list,
        gof_prod_phenofit, df_obs = df_obs_nometh, type=type, real=real)
    df_gof <- map(lst, 'gof') %>% melt_list('prod') %>%
        gather(gof, val, -site, -index, -meth, -prod) %>% data.table()

    df_gof <- subset(df_gof, gof %in% c("Bias", "RMSE"))
    df_gof$gof   %<>% factor() %>% mapvalues(c("Bias", "RMSE"), c("Bias (d)", "RMSE (d)"))
    df_gof$index %<>% fix_level()

    # 2. get mean and sd
    # 2.1 for every sites
    stat_d <- df_gof[, .(median = median(val, na.rm = T),
           sd = sd(val, na.rm = T)), .(prod, gof, index)] %>%
        plyr::mutate(label = sprintf("%.1f±%.1f", median, sd)) %>%
        .[, c('median', 'sd') := NULL] %>%
        spread(prod, label) %>%
        .[order(gof), ]
    # 2.2 for every products
    stat_prod <- df_gof[, .(median = median(val, na.rm = T),
                           sd = sd(val, na.rm = T)), .(prod, gof)] %>%
        plyr::mutate(label = sprintf("%.1f±%.1f", median, sd)) %>%
        .[, c('median', 'sd') := NULL] %>%
        spread(prod, label) %>%
        .[order(gof), ]

    plot <- ggplot(df_gof, aes(index, val)) + #, color = index
        stat_summary(fun.data = box_qtl, geom = "errorbar", width = 0.5) +
        geom_boxplot2(notch = TRUE, outlier.shape = NA, coef = 0, width = 0.8) +
        geom_hline(data = data.frame(gof = c("Bias (d)", "Bias (d)", "RMSE (d)"),
                                     yintercept = c(-15, 15, 15)),
                   aes(yintercept = yintercept), color = "red", linetype = 2) +
        geom_hline(data = data.frame(gof = c("Bias (d)"), yintercept = 0),
                   aes(yintercept = yintercept), color = "blue", linetype = 2) +
        facet_grid(gof~prod, scales = "free_y") + #, labeller = label_parsed
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              strip.text  = element_text(margin = margin(0,0,0,0, "cm"), lineheight = 0)) +
        # ggtitle(title) +
        labs(x = NULL, y = NULL, color = NULL, title = title) + #'Phenophase'
        guides(color=FALSE)
    return(list(stat = stat_d, stat_prod = stat_prod, plot = plot))
}

################################################################################
## 1. aggregated diff curve fitting methods
#  Caution: GPP_mod, GPP_vpm need to fix doy bias. Data copied from `Fig2_scatter_plot.rmd`.
df_lst %<>% set_names(c("GPP[mod]", "GPP[obs]", "GPP[vpm]",
                        "MOD13A1_EVI", "MOD13A1_NDVI", "MOD13Q1_EVI", "MOD13Q1_NDVI",
                        "GPP[avg]", "EVI[avg]"))
df_obs <- merge(df_lst$`GPP[obs]`, stations)
vars   <- names(df_obs)[4:(4+18)]

df_obs_nometh      <- df_obs[meth != "ZHANG", lapply(.SD, median, na.rm = T),
    .(site, IGBP, lat, flag, origin),.SDcols = vars]
df_obs_nometh$meth <- "avg"

df_lst_nometh <- llply(df_lst, function(df_sim){
    df <- df_sim[meth != "ZHANG", lapply(.SD, median, na.rm = T),
        .(site, flag, origin),.SDcols = vars]
    df$meth <- "avg"
    return(df)
})

## debug
df_obs  <- df_obs_nometh
df_list <- df_lst_nometh[c("GPP[mod]", "GPP[vpm]", "MOD13A1_EVI", "MOD13Q1_EVI" )]

# df_sim <- df_lst_nometh[[i]] #GPP_mod
pars <- expand.grid(type = c('SOS', 'EOS'), real = c('all', 'TRS1')) %>%
    set_rownames(paste(.$type, .$real, sep = "_"))

# for (i in 1:4){
#     runningId(i)
#     p_boxplot(type = pars$type[i], real = pars$real[i], df_list = df_list)
# }

# selected 4 prods
df_list <- df_lst_nometh[c("GPP[mod]", "GPP[vpm]", "MOD13A1_EVI", "MOD13Q1_EVI" )]
ps <- mlply(pars, p_boxplot, df_list = df_list, .progress = "text") %>% set_names(rownames(pars))
stat_df   <- map(ps, 'stat')
stat_prod <- map(ps, 'stat_prod') %>% melt_list("prod")

writelist_ToXlsx(stat_df, "t02_average_curve_fitting_methods.xlsx")
file = "Fig_8To12_averaged_curve_methods.pdf"

for (i in seq_along(ps)){
    cairo_pdf(sprintf("Fig_8To12_averaged_curve_methods_%d.pdf", i), 9, 5)
    if (i %in% c(2, 4)){
        p <- ps[[i]]$plot + geom_vline(xintercept = c(2, 5) + 0.5, linetype = 1, size = 0.3, color = "grey70")
    }else{
        p <- ps[[i]]$plot + geom_vline(xintercept = c(4, 7) + 0.5, linetype = 1, size = 0.3, color = "grey70")
    }
    print(p)
    dev.off()
}

## tables


# file.show(file)
# averaged GPP and EVI products
df_list <- df_lst_nometh[c("GPP[avg]","EVI[avg]")]
ps_avg  <- mlply(pars, p_boxplot, df_list = df_list, .progress = "text") %>% set_names(rownames(pars))

data     <- map(ps, "data")
data_avg <- map(ps_avg, "data")
# 2:GPP_obs



## 2. density contourf plot of GPP_avg and EVI_avg
# GPP_avg: TRS1.SOS, TRS1.EOS
# EVI_avg: TRS1.SOS, TRS2.EOS
sel_column <- . %>% {
    dt = .[prod == "GPP_avg" & str_detect(index, "TRS1") |
      prod == "EVI_avg" & str_detect(index, "TRS2")]
    dt[, index := NULL]
}

d_avg <- list(SOS = data_avg$SOS_TRS1,
          EOS = data_avg$EOS_TRS1) %>% map(sel_column) %>%
    melt_list("phase") %>% spread(prod, val) %>% na.omit()
d_avg$phase %<>% factor(c("SOS", "EOS"))

d_all <- list(SOS = data_avg$SOS_all,
              EOS = data_avg$EOS_all) %>%
    melt_list("phase") %>% spread(prod, val) %>% na.omit()
d_all$phase %<>% factor(c("SOS", "EOS"))

## 2. main plot
df <- d_avg
upperval <- 100
df    <- merge(stations[!(site %in% c("US-Me6", "IT-SRo"))], df, by = "site") %>% #, "DE-Akm"
    .[abs(GPP_avg) < upperval & abs(EVI_avg) < upperval, ]

dl   <- dlply(df, .(phase, gof), . %>% {.})
pars <- expand.grid(gof = c('Bias', 'RMSE'), phase = c('SOS', 'EOS')) %>%
    add_column(title = sprintf("(%s) %s %s", letters[1:nrow(.)], .$phase, .$gof))
pars$d <- dl

# ps <- mlply(pars, p_density_naked)
ps <- list()
for (i in seq_along(pars)){
    ps[[i]] <- p_density_naked(dl[[i]], pars$title[i])
}
# ps <- mlply(pars, .(phase, gof), p_density_naked)
p  <- do.call(grid.arrange,
              c(ps,list(nrow = 2, ncol = 2,
                        bottom = textGrob("GPP_avg", gp=gpar(fontsize=fontsize, fontface = "bold")),
                        left   = textGrob("EVI_avg", gp=gpar(fontsize=fontsize, fontface = "bold"), rot = 90))))

file <- 'Fig9_density_contour2.pdf'
cairo_pdf(file, 10, 7.5)
grid.arrange(p, lgd, nrow = 2, heights = c(15, 1))#add legend and export to pdf
dev.off()
file.show(file)

## 3.
stations[IGBP == "WSA", IGBP := "SAV"]

#' @param ... other parameters that pass to ggplot, e.g. color = IGBP, shape = IGBP
p_scatter <- function(df, file = "test.pdf", upper_bias = 80, upper_rmse = 100, showtext = TRUE, ...){
    p_common <- function(p){
        d       <- p$data
        # constrain axis equal
        by      <-  quote(.(gof, phase)) # or .(gof, phase)
        df_lims <- get_facetlims(d, GPP_avg, EVI_avg, .(gof, phase))

        p <- p +
            geom_blank(data = df_lims, aes(GPP_avg, EVI_avg, color = NULL, shape = NULL)) +
            geom_abline(slope = 1, color = "red", size = 1.1) +
            geom_point(size = 2, alpha = 1) +
            facet_grid(gof~phase, scales = "free") +
            coord_equal(ratio = 1) +
            scale_shape_manual(values = shapevalues) +
            stat_prob_2d(aes(color=as.factor(..level..), linetype = as.factor(..level..), shape = NULL),
                         breaks = brks, size = linewidth, show.legend = F)
        if (showtext){
            data_txt <- d[(gof == "Bias (d)" & (abs(GPP_avg) > upper_bias | abs(EVI_avg) > upper_bias)) |
                          (gof == "RMSE (d)" & (abs(GPP_avg) > upper_rmse | abs(EVI_avg) > upper_rmse))]
            # print(data_txt)
            p <- p + geom_text_repel(data = data_txt, aes(label = site), show.legend = F)
        }
        p
    }

    # upper <- 80
    upperval <- 250
    df    <- merge(stations[!(site %in% c("US-Me6", "IT-SRo"))], df, by = "site") %>% #, "DE-Akm"
        .[abs(GPP_avg) < upperval & abs(EVI_avg) < upperval, ]

    df1   <- df[gof == "Bias (d)"]
    df2   <- df[gof == "RMSE (d)"]

    shapevalues <- c(16, 17, 15, 3, 7, 8, 4, 6, 10)
    p1 <- ggplot(df1, aes(GPP_avg, EVI_avg, ...)) + #color = IGBP, shape = IGBP
        geom_hline(data = data.frame(gof = c("Bias (d)"), yintercept = 0),
                   aes(yintercept = yintercept), color = "black", linetype = 2) +
        geom_vline(data = data.frame(gof = c("Bias (d)"), xintercept = 0),
                   aes(xintercept = xintercept), color = "black", linetype = 2) +
        guides(color = FALSE, shape = FALSE)
    p1 %<>% p_common()

    p2 <- ggplot(df2, aes(GPP_avg, EVI_avg, ...)) +
        theme(legend.position = "bottom")
    p2 %<>% p_common()

    cairo_pdf(file, 8, 9)
    grid.arrange(p1, p2, nrow = 2, heights = c(10, 13))
    dev.off()
}

p_scatter(d_avg, file = 'F11_avg.pdf', upper_bias = 50, upper_rmse = 100, showtext = TRUE)
# p_scatter(d_all, file = 'F12_all.pdf', upper_bias = 80, upper_rmse = 120, showtext = FALSE)


add_alpha <- function(colorname, alpha = 0.6){
    {col2rgb(colorname)/255} %>% {rgb(.[1], .[2], .[3], alpha = alpha)}
}

# d_gpp <- df_area[prod %in% c("GPP_mod", "GPP_vpm")]

df_quantile$prod %<>% factor(c("GPP_obs", "GPP_vpm", "GPP_mod", "MOD13A1_EVI", "MOD13Q1_EVI"))
d_obs <- df_quantile[prod == "GPP_obs"]
d_obs[, prod := NULL]

d_gpp <- df_quantile[prod %in% c("GPP_mod", "GPP_vpm")] #, "GPP_vpm"
d_evi <- df_quantile[prod %in% c("MOD13A1_EVI", "MOD13Q1_EVI")] #, "GPP_vpm"

cols <- scales::hue_pal()(3)[c(1, 3)]
color <- "black"
colors <- c(cols, "black")
# colors %<>% map_chr(add_alpha)
alpha_top  <- 0.6
alpha_back <- 0.4

## 1. gpp
p1 <- ggplot(d_gpp, aes(doy, mean, fill = prod, color = prod)) +
    # geom_ribbon(data = d_obs, aes(ymin = `10%`, ymax = `90%`), color = color, fill = color, alpha = 0.2) +
    # geom_ribbon(data = d_obs, aes(ymin = `25%`, ymax = `75%`), color = color, fill = color, alpha = 0.4) +
    geom_ribbon(data = d_obs[alpha == 0.5, ], aes(x = doy, ymin = ymin, ymax = ymax, fill = NULL, color = NULL),
                color = "black", alpha = alpha_back) +
    # geom_ribbon(data = d_gpp[alpha == 0.5 & prod == "GPP_obs", ], aes(ymin = ymin, ymax = ymax), alpha = alpha_top) +

    # geom_ribbon(data = d_gpp[alpha == 0.2 & prod == "GPP_mod", ], aes(ymin = ymin, ymax = ymax), alpha = alpha_back) +
    geom_ribbon(data = d_gpp[alpha == 0.5, ], aes(ymin = ymin, ymax = ymax), alpha = alpha_top) +

    # geom_ribbon(aes(ymin = `25%`, ymax = `75%`), alpha = 0.5) +
    # geom_errorbar(aes(ymin = `25%`, ymax = `75%`), alpha = 1) +
    # geom
    # geom_line(data = d_obs, color = "black", size = 1.1) +
    geom_line(size = 1.1, alpha = 1) +
    geom_line(data = d_obs[alpha == 0.5], aes(y = mean, fill = NULL, color = NULL), linetype = 1, size = 1, alpha = 0.6) +
    # geom_line(data = d_obs[alpha == 0.5], aes(y = ymax, fill = NULL, color = NULL), linetype = 2, size = 1, alpha = 0.6) +
    # geom_line(data = d_gpp[alpha == 0.5], aes(y = ymin), linetype = 2, size = 1, alpha = 0.6) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    facet_wrap(~prod, ncol = 1)

## 2. evi
a <- 32
b <- 0.22

# for gpp
fix_scale <- function(x){
    x/a + b
}
d_obs <- df_quantile[prod == "GPP_obs"]
d_obs[, prod := NULL]
vars <- c("mean", "ymin", "ymax")

d_obs[, (vars) := lapply(.SD, fix_scale), .SDcols = vars]

# scale_y_continuous(sec.axis = sec_axis(~(.+b)/a), name = "GPP_obs") +

df <- d_evi
p2 <- ggplot(df, aes(doy, mean, fill = prod, color = prod, ymin = ymin, ymax = ymax)) +
    # geom_ribbon(data = d_obs, aes(ymin = `10%`, ymax = `90%`), color = color, fill = color, alpha = 0.2) +
    # geom_ribbon(data = d_obs, aes(ymin = `25%`, ymax = `75%`), color = color, fill = color, alpha = 0.4) +
    geom_ribbon(data = d_obs[alpha == 0.5, ], aes(x = doy, fill = NULL, color = NULL),
                color = "black", alpha = alpha_back) +
    # geom_ribbon(data = d_gpp[alpha == 0.5 & prod == "GPP_obs", ], alpha = alpha_top) +

    # geom_ribbon(data = d_gpp[alpha == 0.2 & prod == "GPP_mod", ], alpha = alpha_back) +
    geom_errorbar(data = df[alpha == 0.5, ], alpha = alpha_top) +

    # geom_ribbon(aes(ymin = `25%`, ymax = `75%`), alpha = 0.5) +
    # geom_errorbar(aes(ymin = `25%`, ymax = `75%`), alpha = 1) +
    # geom
    # geom_line(data = d_obs, color = "black", size = 1.1) +
    geom_line(size = 1.1, alpha = 1) +
    geom_line(data = d_obs[alpha == 0.5], aes(y = mean, fill = NULL, color = NULL), linetype = 1, size = 1, alpha = 0.6) +
    # geom_line(data = d_obs[alpha == 0.5], aes(y = ymax, fill = NULL, color = NULL), linetype = 2, size = 1, alpha = 0.6) +
    # geom_line(data = d_gpp[alpha == 0.5], aes(y = ymin), linetype = 2, size = 1, alpha = 0.6) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    facet_wrap(~prod, ncol = 1) +
    scale_y_continuous(sec.axis = sec_axis(~(.-b)*a, name = "GPP_obs"))

cairo_pdf("F13_reason.pdf")
grid.arrange(p1, p2)
dev.off()

