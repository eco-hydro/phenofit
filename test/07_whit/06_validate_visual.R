source('test/stable/load_pkgs.R')
library(ggrepel)

lst_tidy <- lst %>% tidy_rough_fitting()
d <- merge(lst_tidy$melt, df, by = c("site", "t"))


get_GOF <- function(d, is_valid = 1, by_site = T){
    grp <- if(by_site) .(meth, site) else .(meth)
    info <- ddply(d[I_valid == is_valid & iters == "iter1"], grp,function(d){
        with(d, GOF(y0, value))
    }) %>% data.table()
}


info <- list(cal = get_GOF(d, 0),
             val = get_GOF(d, 1))
info_df <- melt_list(info, "type")[order(site)]

dmin = 0.2


meth <- c("wHANTS", "wSG", "wWH2")[1]# "wWH",
pdat <- info_df %>% melt(id.vars = c("site", "meth", "type"), variable.name = "index") %>%
    .[, name := sprintf("%s ~ %s", type, index)] %>%
    .[index %in% c("R2", "NSE", "RMSE", "AI")] %>%
    dcast(site+index+type+name~meth, value.var = "value")

eval(parse(text = sprintf("pdat[, `:=`(diff = wWH - %s, kind = 0)]", meth)))
pdat[diff >  dmin, kind := 1]
pdat[diff < -dmin, kind := -1]

pdat$kind %<>% factor(levels = c(-1, 0, 1), labels = c("bad", "mid", "good"))
pdat[, table(kind), .(index)]
# table(pdat$kind)


p1 <- ggplot(pdat[kind != "mid",  ],
       aes_string(meth, "wWH", color = "kind")) +
    geom_abline(slope = 1, col = "red", size = 1) +
    geom_point(data = pdat, aes(color = NULL), color = "black") +
    geom_point() +
    # geom_text_repel(aes(label = site), show.legend = F) +
    facet_wrap(~name, scale = "free", nrow = 2) +
    theme(aspect.ratio=1) +
    ggtitle(basename(outfile))

p2 <- ggplot(pdat[kind != "mid" & index == "NSE",  ],
       aes_string(meth, "wWH", color = "kind")) +
    geom_abline(slope = 1, col = "red", size = 1) +
    geom_point(data = pdat, aes(color = NULL), color = "black") +
    geom_point() +
    geom_text_repel(aes(label = site), show.legend = F) +
    # facet_wrap(~index, scale = "free") +
    theme(aspect.ratio=1) +
    coord_cartesian(xlim = c(-0.2, 1), ylim = c(-0.2, 1))

print(p2)
print(p1)


# +  +
#     geom_point(data = , color = "green") +
#      +
#     geom_point(data = pdat[wWH - wWH2 < -dmin, ], color = "red") +
#     geom_text_repel(data = pdat[wWH - wWH2 < -dmin, ], aes(label = site), color = "red")

# meth       RMSE       NSE        R2        MAE        AI         Bias    Bias_perc         R pvalue n_sim
# 1 wHANTS 0.05392916 0.8887526 0.8951324 0.03717092 0.9690407 -0.011877982 -0.036886719 0.9461144      0 13650
# 2    wSG 0.08038025 0.7528611 0.7549066 0.05399808 0.9290617 -0.001550526 -0.004815114 0.8688536      0 13650
# 3    wWH 0.06589713 0.8338977 0.8366325 0.04496132 0.9510410 -0.005110335 -0.015869993 0.9146762      0 13650
# 4   wWH2 0.06503765 0.8382024 0.8423709 0.04347108 0.9529469 -0.008935293 -0.027748285 0.9178077      0 13650

# brks   <- season(INPUT, nptperyear,
#                FUN = whitsmw2, wFUN = wFUN, iters = 2,
#                lambda = lambda,
#                IsPlot = IsPlot, plotdat = d,
#                south = d$lat[1] < 0,
#                rymin_less = 0.6, ymax_min = ymax_min,
#                max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5) #, ...
# get growing season breaks in a 3-year moving window

# d_rough_flux <- tidy_rough_fitting(lst)
# d[I_valid == 1, GOF(y0, value), .(meth)]


tidy_rough_fitting <- function(file){
    if (is.character(file)){
        x <- readRDS(file)
    } else{
        x <- file
    }

    d_rough <- x  %>% map(function(l){
        map(l, "whit") %>% melt_list("site")
    }) %>% melt_list("meth")

    d_melt <- d_rough %>% .[, .(site, t, iter1 = ziter1, iter2 = ziter2, meth)] %>%
        melt(id.vars = c("site", "t", "meth"),
             measure.vars = c("iter1", "iter2"), variable.name = "iters") %>%
        .[, .(site, t, iters, value, meth)]

    list(rough = d_rough, melt = d_melt)
}


#
# d_rough_flux <- tidy_rough_fitting("data_test/phenoflux166_rough_val.RDS")
# d_rough_cam  <- tidy_rough_fitting("data_test/phenocam133_rough_val.RDS")

# save(d_rough_cam, d_rough_flux, file = "data_test/phenofit_rough.rda")
# d2 <- readRDS("data_test/phenocam133_rough.RDS")
