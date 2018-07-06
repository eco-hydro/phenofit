windowsFonts(Times = windowsFont("Times New Roman"),
             Arial = windowsFont("Arial"))
fontsize <- 14
mytheme_grey <- theme_grey(base_size = fontsize, base_family = "Arial") +
    theme(
        # legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
        # legend.position="bottom",
        text = element_text(colour = "black", size = fontsize),
        # axis.title.x = element_blank(),
        # plot.title = element_text(hjust = 0),
        axis.title = element_text(face = "bold", size = fontsize),
        axis.text = element_text(colour = "black", size = fontsize),
        strip.text = element_text(colour = "black", size = fontsize, face = "bold"),
        # axis.text.x = element_text(angle = 0, size = fontsize),
        # legend.margin = unit(0.2, "cm"),
        legend.text = element_text(size = fontsize, face = "bold"),
        legend.title = element_blank(),
        # panel.grid.minor.x = element_blank(),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(size = 0.2),
        plot.margin = unit(c(1,3,1,1)*0.2, "cm"))

# theme_grey, theme_gray, theme_light
mytheme_light  <- theme_light(base_size = fontsize, base_family = "Times") +
    theme(
        # legend.position = c(0.02, 0.02), legend.justification = c(0, 0), 
        # legend.position="bottom", 
        text = element_text(colour = "black", size = fontsize), 
        # axis.title.x = element_blank(),
        # plot.title = element_text(hjust = 0),
        axis.title = element_text(face = "bold", size = fontsize), 
        axis.text  = element_text(colour = "black", size = fontsize), 
        strip.text = element_text(colour = "black", size = fontsize, face = "bold"),
        # axis.text.x = element_text(angle = 0, size = fontsize), 
        # legend.margin = unit(0.2, "cm"), 
        legend.text  = element_text(size = fontsize, face = "bold"), 
        legend.title = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        # panel.grid.major = element_line(colour = "white"), 
        panel.grid.minor = element_line(size = 0.2),
        plot.margin = unit(c(1,3,1,1)*0.2, "cm"))
theme_set(mytheme_light)

##
pbar <- function(p){
    p +
        # geom_violin() +
        stat_summary(fun.data = box_qtl, geom = "errorbar", width = 0.5) +
        geom_boxplot2(width = width, notch = T, coef =0) +
        guides(color = F)
}

p_season_prod <- function(info){
    info.meth <- info[, .(RMSE = median(RMSE, na.rm = T)),.(site, meth, prod, phase)]

    fmt <- "%.1f±%.1f" #median ± sd
    mean.meth <- info.meth[, .(median = median(RMSE, na.rm = T),
                               sd = sd(RMSE, na.rm = T),
                               RMSE = quantile(RMSE, 0.53, na.rm = T)), .(meth, prod, phase)] %>%
        mutate(label = sprintf(fmt, median, sd))
    # , color = prod
    p2 <- ggplot(info.meth, aes(meth, RMSE)) %>% pbar +
        geom_blank(data = mean.meth) +
        geom_text(data = mean.meth, aes(meth, median, label = label), vjust=-0.8,
                  show.legend = F, size = 3.4) +
        labs(x = "Curve fitting methods", y = "RMSE (d)") +
        facet_grid(prod~phase, scales = "free")
    p2
    # print(p1)
    # print(p2)
}

p_season <- function(info, title1, title2){
    info.meth <- info[, .(RMSE = median(RMSE, na.rm = T)),.(site, meth)]
    info.prod <- info[, .(RMSE = median(RMSE, na.rm = T)),.(site, prod)]

    fmt <- "%.1f±%.1f" #median ± sd
    mean.meth <- info.meth[, .(median = median(RMSE, na.rm = T),
                               sd = sd(RMSE, na.rm = T),
                               RMSE = quantile(RMSE, 0.53, na.rm = T)), .(meth)] %>%
        mutate(label = sprintf(fmt, median, sd))
    mean.prod <- info.prod[, .(median = median(RMSE, na.rm = T),
                               sd = sd(RMSE, na.rm = T),
                               RMSE = quantile(RMSE, 0.53, na.rm = T)), .(prod)] %>%
        mutate(label = sprintf(fmt, median, sd))

    # , color = meth
    p1 <- ggplot(info.meth, aes(meth, RMSE)) %>% pbar +
        geom_blank(data = mean.meth) +
        geom_text(data = mean.meth, aes(meth, median, label = label), vjust=-1.3,
                  show.legend = F, size = 3) +
        annotate("text", x = -Inf, y = Inf, label = title1, hjust = 0, vjust = 1) +
        labs(x = "Curve fitting methods", y = "RMSE (d)") #, title =title1

    # , color = prod
    p2 <- ggplot(info.prod, aes(prod, RMSE)) %>% pbar +
        geom_blank(data = mean.prod) +
        geom_text(data = mean.prod, aes(prod, median, label = label), vjust=-1.2,
                  show.legend = F, size = 3) +
        theme(axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1)) +
        annotate("text", x = -Inf, y = Inf, label = title2, hjust = 0, vjust = 1) +
        labs(x = "Products", y = "RMSE (d)") #, title =title2
    list(p1, p2)
    # print(p1)
    # print(p2)
}

#1. median of diff curve fitting methods was select
#2.
MEAN_rmse <- function(df_sim){
    # df_sim <- df_lst$MOD13A1_EVI
    I <- match(with(df_obs, paste(meth, site, flag)), with(df_sim, paste(meth, site, flag)))
    df_sim <- df_sim[I]

    vars <- intersect(names(df_sim), names(df_obs))
    d <- list(obs = df_obs[, ..vars], sim = df_sim[, ..vars]) %>%
        melt_list("var") %>%
        gather(index, val, -site, -var, -meth, -flag, -origin) %>% data.table() %>%
        {.[!is.na(val), ]}

    # aggregate phenology metrics by different curve fitting methods
    d  <- d[meth != "ZHANG", .(val = median(val, na.rm = T)), .(site, flag, origin, index, var)]
    ds <- spread(d, var, val) %>% mutate(RE = sim - obs) %>% data.table()
    ds <- ds[!is.na(RE), ]
    info_all <- ddply(ds, .(site, index), function(d){ GOF2(d$RE)})
    return(info_all)
}



