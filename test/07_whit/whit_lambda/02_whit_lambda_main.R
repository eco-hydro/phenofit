source('test/stable/load_pkgs.R')
source('test/07_whit/whit_lambda/smooth_whit_lambda.R')
# source('R/smooth_whit_lambda.R')
# source('test/07_whit/whit_lambda/02_whit_lambda_main.R')
# source('R/smooth_whit.R')
# source('test/GEE/V-pack.r')

nptperyear <- 23

# file ="data_test/whit_lambda/MOD13A1_st_1e3_20180725.rda"
file ="data_test/whit_lambda/MOD13A1_st_1e3_20180731.rda"

if (file.exists(file)){
    load(file)
}else{
    library(sf)

    indir <- "data_test/whit_lambda/raw-csv/mask/"
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
    st <- read_sf("data_test/whit_lambda/shp/st_1e3_mask.shp") %>% data.table() %>% .[, 1:3]

    df <- merge(st, dt, by = "index")
    df <- df[order(site), .(site, y = EVI/1e4, t, w, SummaryQA)]

    # remove sites less then 3y valid data
    info <- df[, .N, site][order(N)]
    site_rm <- info[N < nptperyear*3, site]
    I_rm <- which(df$site %in% site_rm)
    df <- df[-I_rm, ]

    save(df, st, file = file)
}

sites      <- unique(df$site) %>% set_names(., .)
sitename <- sites[2]

################################################################################

#' fill_missdate
#' fill missing date
fill_missdate <- function(){
    years <- 2000:2018
    doy   <- seq(1, 366, 16)
    date  <- sprintf("%4d%03d", rep(years, each = 23), doy) %>% parse_date_time("%Y%j") %>% date()
    if (years[1] == 2000) date <- date[-(1:3)]
    date  <- date[1:(length(date)-12)] # for 2018
}

# fill missing template
temp = data.table(t = date, site = rep(sites, rep(length(date), length(sites))))
    # t as here is image date, other than pixel data.
df_org = merge(df, temp, by = c("t", "site"), all = T) # fill missing values
################################################################################
# 1. I need to know whether lambda values are significant different among
# different \code{dt}.
#
# years <- 2000:2017

## parameters
# deltaT    = 1
# is_extent = F

# subfix <- sprintf("_grp%d", ifelse(is_extent, deltaT, 0))

deltaTs    = c(1, 2, 3, 4, 6, 18)
is_extends = c(T, F)

pars <- expand.grid(deltaT = deltaTs, is_extend = is_extends)

# i <- 1:nrow(pars) %>% set_names(pars$deltaT, )

# set extent = false, it will not enclude previous and subsequent year' data.
optim_lambda_FUN <- function(sitename, wFUN = wSELF, deltaT = 1, extend){
    optim_lambda(sitename, df = df_org, deltaT, extend,
                 IsPlot = F, IsSave = F, file = "whit_formual_wBisquare.pdf",
                 wFUN = wFUN)
}

for (i in 1:6){ # length(deltaTs)
    deltaT <- deltaTs[i]

    for (j in 1:length(is_extends)){
        is_extend  <- is_extends[j]
        extend_str <- ifelse(is_extend, "Extend", "nonExtend")
        subfix <- sprintf("grp%02d_%s", deltaTs[i], extend_str)
        print(subfix)

        # optim_lambda_FUN(sitename)
        res <- par_sbatch(sites, optim_lambda_FUN, wFUN = wTSM, deltaT, is_extend,
                          return.res = T, Save = T,
                          outdir = paste0("result/whit_lambda/v013/whit2", subfix))
    }
}

# a <- optim_lambda_FUN(sitename)
# res <- optim_lambda_FUN(102)
# deltaT <- 1 # current is 4 at GEE
# res.bisquare <- optim_lambda(sitename, df, deltaT = 1, extent = T, IsPlot = F, IsSave = F,
#                     wFUN = wBisquare, file = "whit_formual_wBisquare.pdf")
# res.TSM <- optim_lambda(sitename, df, deltaT = 1, extent = T, IsPlot = F, IsSave = T,
#                     wFUN = wTSM, file = "whit_formual_wTSM.pdf")
# res.self <- optim_lambda(sitename, df, deltaT = 1, extent = T, IsPlot = F, IsSave = T,
#                     wFUN = wSELF, file = "whit_formual_wSELF.pdf")

# group = F # three year group
# outdir <- ifelse(group, '_grp', '')
# outdir <- paste0("result", outdir)

# cpus_per_node <- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE'))
#par_sbatch(sites, optim_lambda_FUN, wFUN = wBisquare, Save = T,
 #          outdir = paste0("result/whit_lambda/wBisquare", subfix) )
# par_sbatch(sites, optim_lambda_FUN, wFUN = wTSM, Save = T,
#            outdir = paste0("result/whit_lambda/wTSM", subfix) )

# res <- par_sbatch(sites, optim_lambda_FUN, wFUN = wBisquare,
#                   return.res = F, Save = T,
#             outdir = paste0("result/whit_lambda/wBisquare", subfix))
# res <- par_sbatch(sites, optim_lambda_FUN, wFUN = wTSM,
#                   return.res = F, Save = T,
#             outdir = paste0("result/whit_lambda/wTSM", subfix))
# res <- par_sbatch(sites, optim_lambda_FUN, wFUN = wSELF,
#                   return.res = F, Save = T,
#             outdir = paste0("result/whit_lambda/wSELF", subfix))

# a <- readRDS("result/whit_lambda/wSELF_grp1/results_0.RDS")
# a <- get_sbatch("result/whit_lambda/wSELF_grp1/")
#
# sapply(a, length) %>% {which(. == 1)}
#l
# res <- list()
# for (i in seq_along(sites)[101:1000]){
#     sitename <- sites[i]
#     res[[i]] <- optim_lambda_FUN(sitename, wSELF)
# }
# system.time({ res <- pbmclapply(sites, optim_lambda, df = df, mc.cores = cpus_per_node) })

# x$`system:index` %<>% str_sub( 1, 31)
