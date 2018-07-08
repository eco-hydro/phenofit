source('test/stable/load_pkgs.R')
source("test/07_whit/main_phenofit_test.R")
source('test/stable/ggplot/geom_boxplot_no_outlier.R')
# source("R/plot_phenofit.R")

# stations212 <- fread("F:/Github/MATLAB/PML/data/flux-212.txt")

infile     <- file_cam
dir_gdrive <- "D:/Document/GoogleDrive/"

# df_cam  <- get_phenofit_result(file_cam)
# df_flux <- get_phenofit_result(file_flux)
#
# st_cam  <- fread(file_st_cam)
# st_flux <- fread(file_st_flux)
load("D:/Documents/GoogleDrive/phenofit.rda")

methods <- c('AG', 'BECK', 'ELMORE', 'ZHANG', 'whit_R', 'whit_gee')[-5]

##
i <- 2
prefix <- c("phenoflux", "phenocam")[i]
df <- if(i == 1) df_flux else df_cam
st <- if(i == 1) st_flux else st_cam

df      <- df[iters == "iter2"]
st$IGBPname %<>% factor(IGBPnames_006)

# make sure different curve fitting methods have the same length fitting
formula <- if(i == 1) formula(site+date+t+y+GPP_NT+GPP_DT+SummaryQA~meth) else
    formula(site+date+t+y+gcc+vci+SummaryQA~meth)

over_perform(df, formula, prefix)

# site figure data input
qc_levels <- c("good", "margin", "snow/ice", "cloud")
qc_colors <- c("grey60", "#00BFC4", "#F8766D", "#C77CFF") %>% set_names(qc_levels)
qc_shapes <- c(19, 15, 4, 17) %>% set_names(qc_levels)

df_trim <- dcast(df, formula, value.var = "value", fun.aggregate = mean) # %>% na.omit()
df_trim$SummaryQA %<>% factor(qc_levels)
# df_trim <- melt(df_trim, measure.vars = methods, variable.name = "meth")

sites <- unique(df$site)
sitename <- sites[100]

vars <- c("get_range", "save_pdf", "lgd_vci", "lgd_gpp",
          "qc_levels", "qc_colors", "qc_shapes", "methods")
cl <- cluster_Init(pkgs = c("data.table", "ggplot2", "magrittr"))
clusterExport(cl, vars)
res <- parLapplyLB(cl, sites, plot_whit,
                   df_trim, st, prefix_fig = paste0("whit_", prefix))
for (i in seq_along(sites)){
    runningId(i)
    sitename <- sites[i]
    plot_whit(sitename, df_trim, st, prefix_fig = paste0("whit_", prefix))
}

# merge_pdf('../whit_phenoflux166.pdf', indir = "./")
# merge_pdf('whit_phenoflux166.pdf', indir = "Figure/")
merge_pdf('whit_phenocam133.pdf', indir = "Figure", del = F)
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
