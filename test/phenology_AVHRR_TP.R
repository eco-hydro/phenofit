library(foreach)
library(iterators)
library(Cairo)
# setwd("..")
source("test/load_pkgs.R")

# 1. rep26
# file_json <- "test/setting_AVHRR_TP_rep26.json"
# options   <- setting.read(file_json)
#
# InitCluster(6)
# system.time(r <- phenofit_process(options, nsite=-1, .parallel = TRUE))
# save(r, file = "phenofit_AVHRR_rep26.rda")

# 2. TP TSF
{
    killCluster()

    file_json <- "test/setting_AVHRR_TP_TSF.json"
    options   <- setting.read(file_json)

    InitCluster(8)
    system.time(r <- phenofit_TSM.avhrr(options, nsite=-1, outdir = "test/AVHRR_TP/",
                                        .parallel = TRUE, dateRange = NULL))
    save(r, file = "phenofit_AVHRR_TP.rda")
}

# foreach(x = 1:8, .packages = "phenofit") %dopar% {
#     phenofit_finefit
# }

# load("phenofit_AVHRR_rep26.rda")
Fig1_rep = FALSE
if (Fig1_rep) {
    title.ylab <- "NDVI"
    NITER <- r[[1]]$fit[[1]]$fFIT[[1]]$zs %>% length()
    lines_colors <- iter_colors(NITER)
    lgd <- make_legend(paste0("iter", 1:NITER), lines_colors, nmax_points = 4)

    # 1. plot all figures
    Fig1 <- FALSE
    if (Fig1){
        methods <- NULL
        CairoPDF("phenofit_AVHRR_rep26.pdf", 11, 7)
        foreach(p = r, i = icount(2), name = names(r)) %do% {
            newpage <- ifelse(i == 1, FALSE, TRUE)
            phenofit_plot(p, "fitting", methods, title = name, title.ylab, newpage = newpage)
        }
        dev.off()
    }

    # 2. plot only one method
    methods <- c("Beck")
    ps <- foreach(p = r, i = icount(), name = names(r)) %do% {
        newpage <- ifelse(i == 1, FALSE, TRUE)
        phenofit_plot(p, "fitting", methods, title = name, title.ylab,
                      Isplot = FALSE, show.legend = FALSE, newpage = newpage)
    }
    FigsToPages(ps, lgd, "NDVI", file = "b.pdf", width = 10, nrow=6)
}
