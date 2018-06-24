library(phenofit)
library(Ipaper)

library(data.table)
library(plyr)
library(lubridate)
library(magrittr)
library(purrr)

source("test/GEE/pkg_seasonality.R")

get_phenofit <- function(i, df){
    d     <- df[Id == i, ] # get the first site data
    titlestr <- with(d[1, ], sprintf('[%03d]%s, lat = %5.2f, lon = %6.2f',
                                     Id, IGBPname, lat, lon))
    tryCatch({
        dnew  <- add_HeadTail(d)
        # 1. Check input data and initial parameters for phenofit
        INPUT <- check_input(dnew$t, dnew$y, dnew$w, maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
        INPUT$y0 <- dnew$y
        IsPlot   <- FALSE # for brks
        # 2. The detailed information of those parameters can be seen in `season`.
        lambda <- init_lambda(INPUT$y)#*2
        # brks   <- season(INPUT, nptperyear,
        #                FUN = whitsmw2, wFUN = wFUN, iters = 2,
        #                lambda = lambda,
        #                IsPlot = IsPlot, plotdat = d,
        #                south = d$lat[1] < 0,
        #                rymin_less = 0.6, ymax_min = ymax_min,
        #                max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5) #, ...
        # get growing season breaks in a 3-year moving window
        brks2 <- season_3y(INPUT, nptperyear, south = d$lat[1] < 0, FUN = whitsmw2,
                           IsPlot = IsPlot, print = print, partial = F)

        # 3. curve fitting
        fit  <- curvefits(INPUT, brks2, lambda =lambda,
                          methods = c("AG", "zhang", "beck", "elmore"), #,"klos",, 'Gu'
                          nptperyear = nptperyear, debug = F, wFUN = wTSM,
                          nextent = 2, maxExtendMonth = 3, minExtendMonth = 1,
                          qc = as.numeric(dnew$SummaryQA), minPercValid = 0.2,
                          print = print)
        fit$INPUT   <- INPUT
        fit$seasons <- brks2

        # svg("Figure1_phenofit_curve_fitting.svg", 11, 7)
        file_pdf = sprintf('Figure/phenofit_v3_%03d.pdf', d$Id[1])
        Cairo::CairoPDF(file_pdf, 11, 6) #
        # grid::grid.newpage()
        plot_phenofit(fit, d, titlestr) %>% grid::grid.draw()# plot to check the curve fitting
        dev.off()

        # temp <- ExtractPheno(fit$fits$ELMORE[1:5], IsPlot = T) # check extracted phenology
        ## 3.2 Get GOF information
        stat  <- ldply(fit$fits, function(fits_meth){
            ldply(fits_meth, statistic.phenofit, .id = "flag")
        }, .id = "meth")
        fit$stat <- stat

        # 4. extract phenology, check extracted phenology extraction for one method.
        # ratio = 1.15
        # file <- "Figure5_Phenology_Extraction_temp.pdf"
        # cairo_pdf(file, 8*ratio, 6*ratio)
        # temp <- ExtractPheno(fit$fits$ELMORE[2:6], IsPlot = T)
        # dev.off()
        # file.show(file)

        # pheno: list(p_date, p_doy)
        p <- lapply(fit$fits, ExtractPheno)

        pheno  <- map(p, tidyFitPheno, origin = INPUT$t[1]) %>% purrr::transpose()

        fit$pheno  <- pheno
        return(fit)
    }, error = function(e){
        message(sprintf('[e]: %s, %s', titlestr, e$message))
    })
}

################################################################################
# lambda     <- 5    # Whittaker parameter
ypeak_min    <- 0.1  # the maximum ymax shoud be greater than `ymax_min`
rytrough_max <- 0.8  # trough < ymin + A*rymin_less
nptperyear   <- 23   # How many points for a single year
wFUN         <- wTSM # Weights updating function, could be one of `wTSM`, 'wBisquare', `wChen` and `wSELF`.

levels <- c(" good", " margin", " snow&ice", " cloud")
df     <- fread("vignettes/phenofit_MOD13A1_st_20.csv", strip.white = F)
df[, `:=`( t = ymd(t),
    SummaryQA = factor(SummaryQA, levels))]

Ids        <- unique(df$Id)
Id         <- df$Id[1]
d          <- df[Id == Ids[1], ] # get the first site data

print = T
# i = 248; a <- get_phenofit(i, df)
par_sbatch(Ids, get_phenofit, df = df, save = T, outdir = "result_v3")

# library(parallel)
# par_sbatch <- function (X, FUN, ..., return.res = F, save = F, outdir = "."){
#     nparams <- length(X)
#     I_node <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#     nodes <- as.numeric(Sys.getenv("SLURM_JOB_NUM_NODES"))
#     cpus_per_node <- as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))

#     if (is.na(I_node))
#         I_node <- 0
#     if (nparams < cpus_per_node * nodes) {
#         nchunk <- cpus_per_node
#     }
#     else {
#         nchunk <- ceiling(nparams/nodes)
#     }

#     I_beg <- I_node * nchunk + 1
#     I_end <- min((I_node + 1) * nchunk, nparams)
#     if (I_beg > I_end) {
#         fprintf("It is empty in this node!")
#         return(NULL)
#     }
#     print(I_node, nodes, cpus_per_node)

#     print('debug1')

#     res <- tryCatch({
#         print(names(X[I_beg:I_end]))
#         mclapply(X[I_beg:I_end], FUN, mc.cores = cpus_per_node,
#             ...)
#     }, error = function(e) {
#         message(sprintf("[error]: %s", e$message))
#     })
#     if (save) {
#         if (!dir.exists(outdir)) dir.create(outdir)
#         saveRDS(res, file = paste0(outdir, "/results_", I_node,
#             ".RDS"))
#     }
#     if (return.res)
#         return(res)
#     NULL
# }
