# Update 20181024
# ---------------
# This function is used to test the rough curve fitting methods, i.e. (wSG,
# wHANTS, wWHd, wWH2) at 16000 points, sampled at global scale.
# 
# After this, pls run:
# - Fig09_(a)_GET_fitting&gof_data.R
# - Fig09_(b)_tidy_gof_data.R
#
# Update 20180910
# ---------------
## 1. Bug found about Whittaker lambda formula
# Parameter for Whittaker lambda was not updated according to the latest result.
#
# v013 : check_fit first, get statistics of this checked data
# v014 : check_fit first, get statistics of original data
# -----
# v013 is more reasonable. Because, contaminated points are often negative bias.
## 2. HANTs weights updating not work
# wTSM not work for HANTS, unknown reason, need to check in the future
#
## Examples
# a <- rough_fitting(sitename, df, st, get(method))
# a$whit[, .(w = witer1 - witer2, z = ziter1 - ziter2)] %>% unique

source('test/stable/load_pkgs.R')
source('test/07_whit/main_gee_Whittaker.R')

adj.param = TRUE
version   = paste0("v3_", 
    ifelse(adj.param, "adjparam", "noadjparam"))

date <- fill_missdate()

################################################################################
is_flux  <- F

if (is_flux){
    # about 2 minutes
    source("test/07_whit/dat_flux&cam_phenofit.R") # load data at flux sites
    df <- df_org # Get data
}else{
    # about 1 hour
    load("data_test/whit_lambda/MOD13A1_st_1e3_20180731.rda")
}
dir_flux <- ifelse(is_flux, "flux/", "")

## examples
sites    <- unique(df$site) %>% set_names(., .)
sitename <- sites[1]

# 1.1 lambda formula coefs

noise_percs = c(0.1, 0.3, 0.5)
noise_perc  = 0 # default is zero

methods  <- c("wHANTS", "wSG", "wWHIT", "wWHIT")
methods2 <- c("wHANTS", "wSG", "wWH", "wWH2")

lst <- list()
k = 4 # k = 4 is corresponding to `grp01_Extend`

source("R/season_3y.R")
source("R/curvefits.R")

# for (k in 4){ # 1:nrow(coefs)
    # noise_perc <- noise_percs[k]
    # df <- select_valid(df, noise_perc = noise_perc)[, 1:10]
pattern <- rownames(coefs)[k]
param   <- as.list(coefs[k, ])

############################################################################
for (i in 1:4){ # only wWH2 this time
    method <- methods[i] #"sgfitw", "whitsmw2" and "wHANTS".
    FUN    <- get(method, envir = as.environment("package:phenofit"))

    runningId(k, prefix = method)

    if (IsPlot){
        file <- sprintf("%s_%s_%2d%%.pdf", prefix, methods2[i], noise_perc*100)
        CairoPDF(file, 8, 10)
        par(mar = c(1.5, 2, 2, 1), mgp = c(3, 0.6, 0), mfrow = c(5, 1), ann = F)
    }

    lambda <- NULL
    if (i >= 4) lambda <- 2
    # a <- llply(sites[23:100], rough_fitting,
    #              df = df, st = st, FUN = get(method), lambda = lambda,
    #              .progress = "text") #lst[[i]]
    # rough_fitting(sitename, df, st, .FUN = wWHIT, lambda = lambda)

    outdir <- sprintf("%sresult/valid_%s/%s%s_%d%%",
                      dir_flush, version, dir_flux, methods2[i], noise_perc*100)
    # print(outdir)
    # outdir <- sprintf("/flush1/kon055/result/whit_lambda/%s_0", methods2[i])
    # outdir <- sprintf("/flush1/kon055/result/valid_whit_lambda2/%s", pattern)
    temp <- par_sbatch(sites, rough_fitting,
                       df = df, st = st, .FUN = get(method), lambda = lambda,
                       return.res = F, Save = T, outdir = outdir, overwrite = T)

    if (IsPlot) dev.off()
    # brks2 <- season_3y(INPUT, nptperyear, south = sp$lat[1] < 0, FUN =FUN,
    #                    nf = nf, frame = frame,
    #                    IsPlot = IsPlot, print = print, partial = F)
}
# }

## SHOW THE GOF Immediately
get_roughFitting_GOF <- function(){
    files <- dirname(outdir) %>% list.files(recursive = T, full.names = T) %>% {
        names <- basename(dirname(.)) %>% gsub("_0%", "", .)
        set_names(., names)
    }

    lst <- llply(files, get_Fitting, .progress = "text")
    df_rough <- melt_list(lst, "meth")
    df_sim   <- merge(df, df_rough)
    info     <- df_sim[, as.list(GOF_extra2(y, value)[-15]), .(meth)]
    info
}

# info <- get_roughFitting_GOF()
############################################################################

# sitename <- "GF-Guy"
# a <- rough_fitting(sitename, df, st, get(method))
# plot(a$whit$ziter2, type = "b")

# files <- dir(outdir, recursive = T, full.names = T)
# names <- basename(dirname(files)) %>% gsub("_0", "", .)

# lst <- llply(files, get_Fitting, .progress = "text") %>% set_names(names)
# d_rough <- melt_list(lst, "meth")

# save(d_rough, file = "rought_fitting.rda")
# lst %<>% set_names(methods2)
# save(df, lst, file = outfile)

# i = 1
# if (i == 1){
#     infile <- file_flux
#     st     <- fread(file_st_flux)
#     outfile <- sprintf("data_test/phenoflux166_rough_val_%02d%%.rda", noise_perc*100)
#     prefix  <- "phenoflux"
# } else{
#     infile <- file_cam
#     st     <- fread(file_st_cam)
#     outfile <- sprintf("data_test/phenocam133_rough_val_%02d%%.rda", noise_perc*100)
#     prefix  <- "phenocam"
# }

# sitename  <- sites[21]
# lst <- llply(sites, rough_fitting,
#              df = df, st = st, FUN = get(method),
#              .progress = "text")
# rough_fitting(sitename, df, st, FUN = get(method))
# df <- fread(infile) # , strip.white = T
# setkeyv(df, c("site", "t"))
# no noise in this version

# lst %<>% set_names(methods2)
# save(df, lst, file = outfile)
