
source('test/stable/load_pkgs.R')
source('R/smooth_whit.R')
source('R/smooth_whit_lambda.R')
# source('test/GEE/V-pack.r')


file ="data_test/whit_lambda/MOD13A1_st_1e3_20180725.rda"

if (file.exists(file)){
    load(file)
}else{
    library(sf)

    indir <- "data_test/whit_lambda/raw-csv/"
    files <- dir(indir, "*.csv", full.names = T)

    lst <- llply(files, fread, .progress = "text")
    dt <- do.call(rbind, lst)

    dt[, `:=`(
        t = str_sub(`system:index`, 1, 10) %>% ymd,
        index = str_sub(`system:index`, 12, 31),
        w = qc_summary(SummaryQA, wmin = 0.2),
        SummaryQA = factor(SummaryQA, levels = qc_values, labels = qc_levels),
        `system:index` = NULL,
        .geo = NULL
    )]
    setkeyv(dt, c("index", "t"))

    # site info
    st <- read_sf("data_test/whit_lambda/shp/st_1e3_gee.shp") %>% data.table() %>% .[, 1:3]

    df <- merge(st, dt, by = "index")
    df <- df[order(site), .(site, y = EVI/1e4, t, w, SummaryQA)]

    save(df, st, file = file)
}

sites      <- unique(df$site)
nptperyear <- 23

sitename <- sites[2]
# 1. I need to know whether lambda values are significant different among
# different \code{dt}.
#
# years <- 2000:2017

deltaT = 1
subfix <- sprintf("_grp%d", deltaT)

optim_lambda_FUN <- function(sitename){
    optim_lambda(sitename, df, deltaT = deltaT, extent = T, IsPlot = F, IsSave = F,
                 wFUN = wBisquare, file = "whit_formual_wBisquare.pdf")
}

# deltaT <- 1 # current is 4 at GEE
res.bisquare <- optim_lambda(sitename, df, deltaT = 1, extent = T, IsPlot = F, IsSave = T,
                    wFUN = wBisquare, file = "whit_formual_wBisquare.pdf")
# res.TSM <- optim_lambda(sitename, df, deltaT = 1, extent = T, IsPlot = F, IsSave = T,
#                     wFUN = wTSM, file = "whit_formual_wTSM.pdf")
# res.self <- optim_lambda(sitename, df, deltaT = 1, extent = T, IsPlot = F, IsSave = T,
#                     wFUN = wSELF, file = "whit_formual_wSELF.pdf")

# group = F # three year group
# outdir <- ifelse(group, '_grp', '')
# outdir <- paste0("result", outdir)

# cpus_per_node <- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE'))
# par_sbatch(sites, optim_lambda_FUN, save = T,
#            outdir = paste0("result/whit_lambda/wBisquare", subfix) )
# par_sbatch(sites, optim_lambda_FUN, save = T,
#            outdir = paste0("result/whit_lambda/wTSM", subfix) )
# par_sbatch(sites, optim_lambda_FUN, save = T,
#            outdir = paste0("result/whit_lambda/wSELF", subfix) )
#
# system.time({ res <- pbmclapply(sites, optim_lambda, df = df, mc.cores = cpus_per_node) })
