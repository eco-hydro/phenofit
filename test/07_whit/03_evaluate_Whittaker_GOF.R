source('test/stable/load_pkgs.R')
source("test/07_whit/main_phenofit_test.R")
source('test/stable/ggplot/geom_boxplot_no_outlier.R')
# source("R/plot_phenofit.R")

# stations212 <- fread("F:/Github/MATLAB/PML/data/flux-212.txt")

infile     <- file_cam
dir_gdrive <- "D:/Document/GoogleDrive/"

file <- "data_test/phenofit.rda"
if (!file.exists(file)){
    df_cam  <- get_phenofit_result(file_cam)
    df_flux <- get_phenofit_result(file_flux)

    st_cam  <- fread(file_st_cam)
    st_flux <- fread(file_st_flux)
    save(df_cam, df_flux, st_cam, st_flux, file = file)
}else{
    load(file)
}

# load("D:/Documents/GoogleDrive/phenofit.rda")

methods <- c('AG', 'BECK', 'ELMORE', 'ZHANG', 'whit_R', 'whit_gee')[-5]

##
i <- 1
prefix <- c("phenoflux", "phenocam")[i]
df <- if(i == 1) df_flux else df_cam
st <- if(i == 1) st_flux else st_cam

st$IGBPname %<>% factor(IGBPnames_006)
st <- st[order(IGBPname,site), ] # reorder according to IGBP
st[site %in% sites_sel, titlestr := sprintf("(%s) %s, %s", letters[1:.N],site, IGBPname)]

# make sure different curve fitting methods have the same length fitting
formula <- if(i == 1) formula(site+date+t+y+GPP_NT+GPP_DT+SummaryQA+iters~meth) else
    formula(site+date+t+y+gcc+vci+SummaryQA+iters~meth)
IGBP.all = F
if (i == 1){
    over_perform(df[iters == "iter2"], formula, prefix, IGBP.all = IGBP.all)
} else{
    over_perform(df[iters == "iter2"], formula, prefix, ylim2 = c(64, 100), IGBP.all = IGBP.all)
}

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
    ps <- list()
    for (i in seq_along(sites)){
        runningId(i)
        sitename <- sites[i]
        # for single method, ggplot obj return
        ps[[i]] <- plot_methods(sitename, df_trim, st, prefix_fig = paste0("whit_", prefix), methods = "whit_gee", show.legend = F)
        # print(p)
        # for all methods, grob obj return
        # p <- plot_methods(sitename, df_trim, st, prefix_fig = paste0("whit_", prefix), methods)
    }

    ylab_r  <- expression("GPP ( gC "*mm^-1*d^-1*" )")
    # ylab_r  <- "VCI"
    FigsToPages(ps, lgd_gpp, ylab_r, file, width = 10)
}
plot_whits(sites, df_trim, "gee_whit_flux166.pdf") # all sites

## 2. select representative points for fluxnet
sites_sel <- c("RU-Fyo", "FR-Pue", "US-UMB", "BE-Vie", "US-KS2", "US-Whs",
    "AU-How", "AU-Dry", "CH-Fru", "US-Los", "DE-Geb")
plot_whits(sites_sel, df_trim, "gee_whit_flux11.pdf") # all sites

# merge_pdf('../whit_phenoflux166.pdf', indir = "./")
# merge_pdf('whit_phenoflux166.pdf', indir = "Figure/")
merge_pdf('whit_phenocam133.pdf', indir = "Figure", "whit_phenocam.*.pdf", del = F)
merge_pdf('whit_phenoflux166.pdf', indir = "Figure", "whit_phenoflux.*.pdf", del = F)

################################################################################
# colnames(df)[4] <- "raw"

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
