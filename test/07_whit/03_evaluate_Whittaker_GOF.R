source('test/stable/load_pkgs.R')
source("test/07_whit/main_phenofit_test.R")
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

    df <- llply(patterns, function(pattern) readwhitMAT(dir_gdrive, pattern),
                .progress = "text") %>% set_names(patterns) %>% melt_list("meth")
    df$meth %<>% as.factor()


    whit_flux <- df[grep("phenoflux", meth), ] %>% data.frame() %>% data.table()
    whit_cam  <- df[grep("phenocam", meth), ] %>% data.frame() %>% data.table()

    whit_flux$meth %<>% sub_levels
    whit_cam$meth  %<>% sub_levels

    # fwrite(df_flux, "data_test/gee_whit_phenoflux166.csv")
    # fwrite(df_cam , "data_test/gee_whit_phenocam133.csv")

    # 2. merge others ---------------------------------------------------------
    lst_cam  <- get_phenofit_fitting(file_cam, "data_test/phenofit_fitting_phenocam133.rda")
    lst_flux <- get_phenofit_fitting(file_flux, "data_test/phenofit_fitting_phenoflux166.rda")

    df_cam  <- get_phenofit_update_whit(lst_cam , whit_cam)
    df_flux <- get_phenofit_update_whit(lst_flux, whit_flux)

    # load
    load("data_test/phenofit_rough.rda")
    df_cam$fits_merge %<>% rbind(d_rough_cam$melt)
    df_cam %<>% get_phenofit_result_merge()

    df_flux$fits_merge %<>% rbind(d_rough_flux$melt)
    df_flux %<>% get_phenofit_result_merge()

    st_cam  <- fread(file_st_cam)
    st_cam <- st_cam[!(site %in% c("cedarcreek", "rosemount"))]
    st_flux <- fread(file_st_flux)
    save(df_cam, df_flux, st_cam, st_flux, file = file)
}else{
    load(file)
}

# FR-Pue, AU-Dry ""

sites_sel <- c("RU-Fyo", "GF-Guy", "US-UMB", "CN-Cha", "US-KS2", "US-Whs",
               "AU-How", "ZA-Kru", "CN-HaM", "US-Los", "DE-Geb")

{
    st_flux[site %in% sites_sel, selected := 1]

    df2sp(st_flux) %>% writePointsShape("st_flux166.shp")
    df2sp(st_cam)  %>% writePointsShape("st_cam133.shp")
}


# load("D:/Documents/GoogleDrive/phenofit.rda")
methods <- c('AG', 'BECK', 'ELMORE', 'ZHANG', 'whit_R', 'whit_gee')[-5]

##

i <- 1
prefix <- c("phenoflux", "phenocam")[i]
df <- if(i == 1) df_flux else df_cam
st <- if(i == 1) st_flux else st_cam

df <- df$all
st$IGBPname %<>% factor(IGBPnames_006)
st <- st[order(IGBPname,site), ] # reorder according to IGBP
st[site %in% sites_sel, titlestr := sprintf("(%s) %s, %s", letters[1:.N],site, IGBPname)]

# make sure different curve fitting methods have the same length fitting
formula <- if(i == 1) formula(site+date+t+y+GPP_NT+GPP_DT+SummaryQA+iters~meth) else
    formula(site+date+t+y+gcc+vci+SummaryQA+iters~meth)
IGBP.all = T
ylim2 <- if(i == 1) c(35, 100) else c(64, 100)

xlab <- st[, .N, IGBPname] %>% rbind(data.table(IGBPname = "all", N = nrow(st)))
xlab[, label:=sprintf("%s\n(n=%2d)", IGBPname, N)]
outfile <- sprintf("Fig8_valid_%s.pdf", prefix)

source("test/07_whit/main_phenofit_test.R")
over_perform(df[iters == "iter2"], formula, prefix, xlab, ylim2, IGBP.all = IGBP.all, outfile)




# site figure data input
df_trim <- dcast(df, formula, value.var = "value", fun.aggregate = mean) # %>% na.omit()
df_trim$SummaryQA %<>% factor(qc_levels)
# df_trim <- melt(df_trim, measure.vars = methods, variable.name = "meth")

sites    <- st$site
# sites <- unique(df$site)
sitename <- sites[100]

# vars <- c("get_range", "save_pdf", "lgd_vci", "lgd_gpp",
#           "qc_levels", "qc_colors", "qc_shapes", "methods")
# cl <- cluster_Init(pkgs = c("data.table", "ggplot2", "magrittr"))
# clusterExport(cl, vars)
# res <- parLapplyLB(cl, sites, plot_whit,
#                    df_trim, st, prefix_fig = paste0("whit_", prefix))

plot_whits <- function(sites, df_trim, file){
    methods <- c('WH_p2',  "WH_p15",'wWH_p2', "wWH_p15",'wWH') #'WH_p15', 'wWH_p15', 'wWH_v13'
    # methods <- "wWH_v13"
    nmeth   <- length(methods)
    # Cairo::CairoPDF("gee_whit_flux166_all_meth_v2.pdf", 9, nmeth*1.8)

    ps <- list()
    for (i in seq_along(sites)){
        runningId(i)
        sitename <- sites[i]
        # for single method, ggplot obj return
        p <- plot_methods(sitename, df_trim, st, prefix_fig = paste0("whit_", prefix),
                          methods = methods, show.legend = T)
        # print(p)
        #
        # ps[[i]] <- p
        # for all methods, grob obj return
        # p <- plot_methods(sitename, df_trim, st, prefix_fig = paste0("whit_", prefix), methods)
    }
    # dev.off()
    ylab_r  <- expression("GPP ( gC "*mm^-1*d^-1*" )")
    # ylab_r  <- "VCI"
    # FigsToPages(ps, lgd_gpp, ylab_r, file, width = 10)
}

st[, ID := 1:.N]
plot_whits(sites, df_trim, "gee_whit_flux166_v2.pdf") # all sites


## 2. select representative points for fluxnet
plot_whits(sites_sel, df_trim, "gee_whit_flux11.pdf") # all sites

# merge_pdf('../whit_phenoflux166.pdf', indir = "./")
# merge_pdf('whit_phenoflux166.pdf', indir = "Figure/")
merge_pdf('whit_phenocam133.pdf', indir = "Figure", "whit_phenocam.*.pdf", del = F)
merge_pdf('whit_phenoflux166.pdf', indir = "Figure", "whit_phenoflux.*.pdf", del = F)

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
