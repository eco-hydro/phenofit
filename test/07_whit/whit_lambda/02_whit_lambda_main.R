source('test/stable/load_pkgs.R')
source('R/smooth_whit.R')
source('R/smooth_whit_lambda.R')
# source('test/GEE/V-pack.r')
file = 'data_whit_lambda_01.csv'

if (!file.exists(file)){
    dt = fread('data_test/temp/MOD13A1_st_1e3.csv')
    dt[, `:=`(y    = EVI/1e4,
              t    = ymd(date),
              w    = qc_summary(dt$SummaryQA))]
    dt[, per := sum(!is.na(EVI))/.N, site]
    df <- dt[per > 0.3, .(site, y, t, w, IGBPcode)]
    fwrite(df, file)
}

df         <- fread(file)
sites      <- unique(df$site)
nptperyear <- 23

# 1. I need to know whether lambda values are significant different among
# different \code{dt}.
#
# years <- 2000:2017

# deltaT <- 1 # current is 4 at GEE
res <- optim_lambda(sitename, df, deltaT = 1, extent = T, IsPlot = F, IsSave = T)

group = F # three year group
outdir <- ifelse(group, '_grp', '')
outdir <- paste0("result", outdir)

# cpus_per_node <- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE'))
par_sbatch(sites, optim_lambda, df = df, save = T, outdir = "result")
# system.time({ res <- pbmclapply(sites, optim_lambda, df = df, mc.cores = cpus_per_node) })
