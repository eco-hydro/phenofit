library(foreach)
library(iterators)
library(Cairo)
# setwd("..")
source("test/load_pkgs.R")
load("D:/Documents/OneDrive - mail2.sysu.edu.cn/phenology/data/00basement_TP.rda")


# 2. TP phenology ---------------------------------------------------------
{
    killCluster()

    outdir <- "D:/Documents/OneDrive - mail2.sysu.edu.cn/phenology/AVHRR_TP_v0.1.2/"
    # outdir <- "F:/phenology/AVHRR_TP_v0.1.2/"
    if (!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE) }

    file_json <- "test/setting_AVHRR_TP_nosnow.json"
    options   <- setting.read(file_json)

    InitCluster(8)
    I_all <- seq_along(gridclip)
    I_lst <- chunk(I_all, 2)
    I_part <- I_lst[[1]] #%>% rev()
    system.time(r <- phenofit_TS.avhrr(options,
        I_part = I_all,
        dateRange = NULL,
        outdir = outdir,
        exportType = "pheno",
        overwrite = FALSE,
       .parallel = TRUE))
    # save(r, file = "phenofit_AVHRR_TP.rda")
}

## post-process
# pheno_mean <- function(file){
#     x <- read_rds(file)
#     # get the mean value of those values
#     d <- x$pheno$doy %>% melt_list("meth")
#     d[, lapply(.SD, sd, na.rm = TRUE), .(meth), .SDcols = colnames(d)[-c(1:2, ncol(d))]]
# }

# 1. rep26 ----------------------------------------------------------------
s1_rep26 = FALSE
# s1_rep26 = TRUE
if (s1_rep26) {
    file_json <- "test/setting_AVHRR_TP_rep26.json"
    options   <- setting.read(file_json)

    InitCluster(6)
    system.time(r <- phenofit_process(options, nsite=-1, .parallel = TRUE))
    save(r, file = "phenofit_AVHRR_rep26.rda")
}

# load("phenofit_AVHRR_rep26.rda")
Fig1_rep = FALSE
# Fig1_rep = TRUE
if (Fig1_rep) {
    title.ylab <- NULL #"NDVI"
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
    FigsToPages(ps, lgd, "NDVI", file = "phenofit_TP_rep26.pdf", width = 10, nrow=6)
}
