rm(list = ls())
source("test/load_pkgs.R")
source("test/07_whit/main_phenofit_test.R")
source("test/07_whit/dat_flux&cam_phenofit.R")
source('test/stable/ggplot/geom_boxplot_no_outlier.R')
# source("R/plot_phenofit.R")

# stations212 <- fread("F:/Github/MATLAB/PML/data/flux-212.txt")

infile     <- file_cam
dir_gdrive <- "D:/Document/GoogleDrive/"

file <- "data_test/phenofit.rda"
# GEE Whittaker result is also read.
if (!file.exists(file)){
    # 1. gee whit ---------------------------------------------------------
    # Modified 20180807, multiple Whittakers version were tested.
    #
    # "gcesapelo" and "coaloilpoint" are located at the seaside.
    # It's normal sometimes nodata.

    ## rename meth
    sub_levels <- . %>% {
        level <- levels(.)
        mapvalues(., level, gsub("phenoflux166_|phenocam133_|_v105", "", level))
    }

    dir_gdrive   <- "D:/Document/GoogleDrive/whit" #data/gee_phenofit/v2/

    files    <- dir(dir_gdrive, "*.geojson", full.names = T)
    patterns <- str_extract(basename(files), ".*(?=_\\d{4}_)") %>% unique()

    df_whit <- llply(patterns, function(pattern) readwhitMAT(dir_gdrive, pattern),
                .progress = "text") %>% set_names(patterns) %>% melt_list("meth")
    df_whit$meth %<>% as.factor()

    # whit_flux <- df[grep("phenoflux", meth), ] %>% data.frame() %>% data.table()
    # whit_cam  <- df[grep("phenocam", meth), ] %>% data.frame() %>% data.table()

    # whit_flux$meth %<>% sub_levels
    # whit_cam$meth  %<>% sub_levels

    # fwrite(df_flux, "data_test/gee_whit_phenoflux166.csv")
    # fwrite(df_cam , "data_test/gee_whit_phenocam133.csv")

    # 2. merge others ---------------------------------------------------------
    lst_cam  <- get_phenofit_fitting(file_cam, "data_test/phenofit_fitting_phenocam133.rda")
    lst_flux <- get_phenofit_fitting(file_flux, "data_test/phenofit_fitting_phenoflux166.rda")

    df_cam  <- get_phenofit_update_whit(lst_cam , df_whit)
    df_flux <- get_phenofit_update_whit(lst_flux, df_whit)

    # load
    load("rought_fitting.rda")
    # load("data_test/phenofit_rough.rda")
    df_cam$fits_merge %<>% rbind(d_rough)
    df_cam %<>% get_phenofit_result_merge()

    df_flux$fits_merge %<>% rbind(d_rough)
    df_flux %<>% get_phenofit_result_merge()

    st_cam  <- fread(file_st_cam)
    st_cam <- st_cam[!(site %in% c("cedarcreek", "rosemount"))]
    st_flux <- fread(file_st_flux)
    save(df_cam, df_flux, st, st_cam, st_flux, file = file)
}else{
    load(file)
}

# FR-Pue, AU-Dry ""

sites_sel <- c("RU-Fyo", "GF-Guy", "US-UMB", "CN-Cha", "US-KS2", "US-Whs",
               "AU-How", "ZA-Kru", "CN-HaM", "US-Los", "DE-Geb")

{
    library(maptools)
    st_flux[match(sites_sel, site), `:=`(
        selected = 1,
        label    = sprintf("(%s), %s", letters[1:.N], IGBPname)
    )]

    # df2sp(st_flux) %>% writePointsShape("F:/Github/GEE/gee_whittaker/ArcGIS/shp/st_flux166.shp")
    # df2sp(st_cam)  %>% writePointsShape("st_cam133.shp")
}

# load("D:/Documents/GoogleDrive/phenofit.rda")
methods <- c('AG', 'BECK', 'ELMORE', 'ZHANG', 'whit_R', 'whit_gee')[-5]

##
source("test/07_whit/main_phenofit_test.R")
i <- 1

get_over_perform <- function(i = 1){
    prefix <- c("phenoflux", "phenocam")[i]
    df <- if(i == 1) df_flux else df_cam
    st <- if(i == 1) st_flux else st_cam

    df <- df$all
    st$IGBPname %<>% factor(IGBPnames_006)
    st <- st[order(IGBPname, site), ] # reorder according to IGBP
    st[, ID := 1:.N]

    st[site %in% sites_sel, titlestr := sprintf("(%s) %s, %s", letters[1:.N],site, IGBPname)]

    sites    <- st$site    # sites <- unique(df$site)
    sitename <- sites[100]

    # visualization -----------------------------------------------------------
    # make sure different curve fitting methods have the same length fitting
    formula <- if(i == 1) formula(site+date+t+y+GPP_NT+GPP_DT+SummaryQA+iters~meth) else
        formula(site+date+t+y+gcc+vci+SummaryQA+iters~meth)
    IGBP.all = T
    ylim2 <- if(i == 1) c(35, 100) else c(64, 100)

    xlab <- st[, .N, IGBPname] %>% rbind(data.table(IGBPname = "all", N = nrow(st)))
    xlab[, label:=sprintf("%s\n(n=%2d)", IGBPname, N)]
    outfile <- sprintf("Fig8_valid_%s.pdf", prefix)

    info_flux <- over_perform(df[iters == "iter2"], st, i, formula, prefix, xlab, ylim2, IGBP.all = IGBP.all, outfile)
    info_flux
}
info_flux <- get_over_perform(1)
table(info_flux$meth)


info_cam  <- get_over_perform(2)


info <- info_flux
d <- dcast(info, site+ID+IGBPname~meth, value.var = "R")
d <- d[IGBPname != "all"]
d <- melt(d[, c(1:8, 10:11)], id.vars = c("site", "ID", "IGBPname", "whit_fluxcam_wWH"),
          variable.name = "meth")
colnames(d)[4] <- "wWH"
d <- d[meth != "raw" & meth != "wWH2"]

d[,`:=`(diff = wWH - value, kind=0)]
dmax <- 0.05

d[diff > dmax,  kind := 1]
d[diff < -dmax, kind := -1]
d$kind %<>% factor(levels = c(1, 0, -1),
                   labels = c("Bigger", "Similar", "Smaller"))
d_dominant <- ddply_dt(d, .(table(kind)), .(meth))
d_dominant[, `:=`(
    text_big   = sprintf("bold('%s' * phantom(' vs ' * '%s'))", Bigger, Smaller),
    text_vs    = sprintf("bold(phantom('%s') * ' vs ' * phantom('%s'))", Bigger, Smaller),
    text_small = sprintf("bold(phantom('%s vs ') * '%s')", Bigger, Smaller))]


# p <- ggplot(d[meth != "raw"], aes(wWH, value, color = kind)) +
#     geom_point() + geom_abline(slope = 1, color = "red", size = 1) +
#     facet_wrap(~meth)
#
# x <- d_dominant %>% melt("meth", variable.name = "kind")



# ggplot(df_trim, aes(y, GPP_NT, color = SummaryQA)) +
#     geom_point(alpha = 0.1) +
#     geom_smooth(aes(color = NULL), method = "lm")
#
# ggplot(df_trim, aes(whit_fluxcam_wWH, GPP_NT, color = SummaryQA)) +
#     geom_point(alpha = 0.1) +
#     geom_smooth(aes(color = NULL), method = "lm")


# d <- list(info_flux, info_cam) %>% set_names(types) %>% melt_list("type")
# d <- d[meth != "raw" & type == "(a) FLUXNET"]
# levels <- levels(d$meth)
# levels_new <- sprintf("(%s) %s", letters[seq_along(levels[-1])], levels[-1])
# levels_new %<>% c("raw", .)
# d$meth %<>% mapvalues(levels, levels_new)
#
# d_dominant <- ddply_dt(d, .(table(kind)), .(meth, type))
# d_dominant[, label := sprintf("%2d vs %2d", Better, Worse)]


# x <- d_dominant %>% melt(c("meth", "type"), variable.name = "kind")

# library(lattice)
# library(latticeExtra)
# xyplot(value~wWH|meth*type, d)
# ggplot(x[kind != "middle"], aes(kind, value)) + facet_grid(meth~type) + geom_point()

## phenocam can't illustrate the merits of wWH, so exclude it.
colors <- scales::hue_pal()(3)
colors <- c(colors[2], "grey", colors[1])
text_size <- 4
p <- ggplot(d, aes(wWH, value, color = kind)) +
    geom_abline(slope = 1, color = "red", size = 0.8) +
    geom_point() +
    facet_wrap(~ meth) +
    geom_text(data = d_dominant, aes(label = text_big, color = NULL), parse = T,
              color = colors[1], fontface = "bold", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, show.legend = F, size = text_size) +
    geom_text(data = d_dominant, aes(label = text_vs, color = NULL), parse = T,
              color = "black", fontface = "bold", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, show.legend = F, size = text_size) +
    geom_text(data = d_dominant, aes(label = text_small, color = NULL), parse = T,
              color = colors[3], fontface = "bold", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, show.legend = F, size = text_size) +
    # geom_text(data = d_dominant, aes(label = label, color = NULL), fontface = "bold",
    #           x = -Inf, y = Inf, hjust = -0.3, vjust = 1.5, show.legend = F) +
    ylab(expression(R^2)) +
    # coord_equal() +
    # theme_gray(base_size = 14) +
    theme_light(base_size = fontsize, base_family = "Arial") +
    theme(
        legend.position = "bottom",
        legend.margin = margin(-15, unit = "pt"),
        legend.title=element_blank(),
        # plot.margin = margin(),
        # axis.title.x = element_text(margin = margin()),
        # legend.position = c(1-0.01, 0.01), legend.justification = c(1, 0),
          panel.grid.major = element_line(linetype = 2),
          panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold", margin = margin(4, 0, 4, 0, unit = "pt")),
          axis.text = element_text(color = "black")) +
    scale_color_manual(values = colors)

p + ggrepel::geom_text_repel(data = d[diff > 0.4 | diff < -0.3], aes(label = site))


##
file = "Fig11_valid_fluxnet.pdf"
file_tiff <- gsub(".pdf", ".tif", file)
width = 7; height = 5.5
write_fig(p, file_tiff, width, height, T, res = 300)
write_fig(p, file, width, height, T)




# site figure data input
df_trim <- dcast(df, formula, value.var = "value", fun.aggregate = mean) # %>% na.omit()
df_trim$SummaryQA %<>% factor(qc_levels)
# df_trim <- melt(df_trim, measure.vars = methods, variable.name = "meth")

# vars <- c("get_range", "save_pdf", "lgd_vci", "lgd_gpp",
#           "qc_levels", "qc_colors", "qc_shapes", "methods")
# cl <- cluster_Init(pkgs = c("data.table", "ggplot2", "magrittr"))
# clusterExport(cl, vars)
# res <- parLapplyLB(cl, sites, plot_whit,
#                    df_trim, st, prefix_fig = paste0("whit_", prefix))

## 1. show all sites
check_Whittaker(sites, df_trim, "gee_whit_flux166_v2.pdf") # all sites
# merge_pdf('../whit_phenoflux166.pdf', indir = "./")
# merge_pdf('whit_phenoflux166.pdf', indir = "Figure/")
merge_pdf('whit_phenocam133.pdf', indir = "Figure", "whit_phenocam.*.pdf", del = F)
merge_pdf('whit_phenoflux166.pdf', indir = "Figure", "whit_phenoflux.*.pdf", del = F)

## 2. select representative points for fluxnet
check_Whittaker(sites_sel, df_trim, "gee_whit_flux11.pdf", "whit_fluxcam_wWH") # all sites



################################################################################
# colnames(df)[4] <- "raw"

## check impact of lambda and weigth updating
names(table(df$meth)) %>% .[grep("*WH*", .)] %>%
    paste(collapse = "', '") %>% paste0("'", ., "'") %>% cat
methods <- c('wWH', 'WH_p15', 'WH_p2', 'wWH_p15', 'wWH_p2')

df[, `:=`(
    re = value - y,
    doy = yday(date)
)]



## 1. autocorrelation shows that weighted function works better, and wWH_p2 is the
# best in all methods.
# Initial lambda is useless here.
df_acf <- ddply(df[, ], .(meth, site), function(d){
    with(d, {
        acf(-re, lag.max = 10, plot = F, na.action = na.pass)$acf[,,1][-1]
    })
}) %>% data.table() %>% melt(c("meth", "site"), variable.name = "lag")
df_acf$lag %<>% factor() %>% as.numeric() %>% as.factor()

ggplot(df_acf, aes(lag, value, color = meth)) + geom_boxplot()


## 2. residual analysis
# Result show that no significant difference
ggplot(df[], aes(as.factor(doy), re)) +
    # geom_hline(yintercept = c(0, 0.05, 0.1)
    #            # color = c("black", "blue", "red")
    #            # linetype = c(1, 2, 1)
    #            ) +
    geom_boxplot2(aes(color = meth)) +
    geom_hline(yintercept = 0, color = "black") +
    geom_hline(yintercept = 0.025/2, color = "blue", linetype =2) +
    geom_hline(yintercept = 0.05/2, color = "red") +
    scale_y_continuous(limits = c(-1, 1)*0.025)+
    facet_wrap(~meth)

## 3.


# grid.newpage(); grid.draw(lgd)
# arrangeGrob(p, lgd, nrow = 2, heights = c(15, 1), padding = unit(1, "line")) #return,

# p0 <- ggplot(d[meth == "raw"], aes(date, value, shape = SummaryQA, color = SummaryQA)) +
#     geom_point() +
#     theme(legend.position = "none") +
#     scale_color_manual(values = c("good" = "grey60", "margin" = "#00BFC4",
#                                   "snow/ice" = "#F8766D", "cloud" = "#C77CFF"), drop = F) +
#     scale_shape_manual(values = c(19, 15, 4, 17), drop = F) +
#     scale_y_continuous(lim = lim_raw)
# p3 <- ggplot_dual_axis(p1, p2) #%>% as.ggplot()
# p <- ggplot_dual_axis(p3, p0, add_yaxis_r = F)
