source('test/stable/load_pkgs.R')
source('test/stable/s1_MCD12Q2.R')
stations <- fread("test/Figs/phenofit_st99.csv")

# main --------------------------------------------------------------------
dirs <- dir("Y:/R/phenofit/OUTPUT/", full.names = T) %>%
    set_names(basename(.))

lst    <- llply(dirs, get_slurm_out, IsSave = F)
df_lst <- llply(lst, function(x) tidy_pheno(rm_empty(x)))

# save(lst, df_lst, file = "phenofit_OUTPUT.rda")
load("test/Figs/phenofit_OUTPUT.rda")
df_obs <- df_lst$GPP_obs
df_obs <- merge(df_obs, stations)
## 1. GOF figures for PMLv2

plot_cmp <- function(type, real){
    gof_prod_phenofit(df_sim, df_obs, type, real, trim = TRUE)
}

products <- names(df_lst)[-2]
pars <- expand.grid(type = c('SOS', 'EOS'), real = c('TRS1', 'all')) %>%
    set_rownames(paste(.$type, .$real, sep = "_"))

type = 'SOS'; real = "all"

res <- list()
for (i in 1:length(products)){
    name <- products[i]
    df_sim <- df_lst[[name]]

    figname <- sprintf('%d. %s_phenofit_gof.pdf', i, name)
    CairoPDF(figname, width = 11.5, height = 6.5, pointsize = 8)
    res[[i]] <- mlply(pars, plot_cmp, .progress = "text") %>% set_names(rownames(pars))
    dev.off()
}
res %<>% set_names(products)

# MCD12Q2 -----------------------------------------------------------------
MCD12Q2 <- fread("file:///Y:/R/phenofit/data/MCD12Q2_flux212_2001-2014.csv")
MCD12Q2$date %<>% ymd()
df <- merge(MCD12Q2, df_obs, by = c("site", "flag"))

CairoPDF("0. MCD12Q2_phenofit_GOF.pdf", width = 11.5, height = 6.5)
MCD12Q2.sos <- gof_prod_MCD12Q2(df, type = 'SOS', trim = F)
MCD12Q2.eos <- gof_prod_MCD12Q2(df, type = 'EOS', trim = F)
# a <- gof_prod_MCD12Q2(df, type = 'SOS', trim = T)
# b <- gof_prod_MCD12Q2(df, type = 'EOS', trim = T)
dev.off()



# ggplot(subset(pdat, GOFIndex %in% c("Bias", "RMSE")), aes(index, val, color = index)) +
#     stat_summary(fun.data = box_qtl, geom = "errorbar", width = 0.5) +
#     geom_boxplot_noOutliers(notch = TRUE, coef = 0) +
# # geom_boxplot() +
# # geom_jitter(width = 0.5) +
# facet_grid(GOFIndex~meth, scales = "free_y") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#     ggtitle(title) +
#     labs(x = 'Phenophase', y = NULL, color = 'Phenophase') +
#     # coord_flip() +
#     guides(color=FALSE)
