library(foreach)
library(iterators)
library(Cairo)
# setwd("..")
source("test/load_pkgs.R")
load("../phenology_TP/data/00basement_TP.rda")
# 1. rep26 ----------------------------------------------------------------

poly <- matrix(c(1,1,1,2,2,2,2,1,1,1),ncol=2, byrow=TRUE) %>%
    list() %>%
    st_polygon()
p1 = st_sample(poly, 6)
st_crs(4326)

s1_rep26 = FALSE
s1_rep26 = TRUE
if (s1_rep26) {
    file_json <- "test/setting_AVHRR_TP_rep26.json"
    options   <- setting.read(file_json)

    InitCluster(6)
    system.time(r <- phenofit_process(options, nsite=-1, .parallel = TRUE))
    save(r, file = "phenofit_AVHRR_rep26.rda")
}

# load("phenofit_AVHRR_rep26.rda")
Fig1_rep = FALSE
Fig1_rep = TRUE
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


# 2. TP phenology ---------------------------------------------------------
{
    killCluster()

    outdir <- "F:/phenology/AVHRR_TP_v0.1.2/"
    if (!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE) }

    file_json <- "test/setting_AVHRR_TP_nosnow.json"
    options   <- setting.read(file_json)

    InitCluster(8)
    system.time(r <- phenofit_TS.avhrr(options, nsite=-1,
                                       dateRange = NULL,
                                       outdir = outdir,
                                       exportType = "pheno",
                                       overwrite = TRUE,
                                       .parallel = TRUE))
    # save(r, file = "phenofit_AVHRR_TP.rda")
}


## post-process
pheno_mean <- function(file){
    x <- read_rds(file)
    # get the mean value of those values
    d <- x$pheno$doy %>% melt_list("meth")
    d[, lapply(.SD, sd, na.rm = TRUE), .(meth), .SDcols = colnames(d)[-c(1:2, ncol(d))]]
}

r <- microbenchmark::microbenchmark(
    pheno_mean(file), times = 1000
)
file <- "D:/Documents/OneDrive - mail2.sysu.edu.cn/phenofit_20267.RDS"



files <- dir("G:/Github/phenology/China_Phenology/data/txt", full.names = T)[-4]
llply(files, function(file){
     fwrite(fread(file), file)
})
