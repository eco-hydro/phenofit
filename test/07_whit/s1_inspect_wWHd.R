# Update 20181024
# ---------------
# This function is used to test the rough curve fitting methods, i.e. (wSG,
# wHANTS, wWHd, wWH2) at 16000 points, sampled at global scale.
source('test/stable/load_pkgs.R')

adj.param = TRUE
version   = paste0("v3_",
                   ifelse(adj.param, "adjparam", "noadjparam"))

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

################################################################################

source("R/season_3y.R")
source("R/curvefits.R")
source('R/smooth_wWHIT_lambda.R')
source('test/07_whit/main_gee_Whittaker.R')

date <- fill_missdate()


# for (k in 4){ # 1:nrow(coefs)
# noise_perc <- noise_percs[k]
# df <- select_valid(df, noise_perc = noise_perc)[, 1:10]

k = 4 # k = 4 is corresponding to `grp01_Extend`
pattern <- rownames(coefs)[k]
param   <- as.list(coefs[k, ])

info  <- readRDS("wWHd_inspect_sites.RDS")

sites <- info$site
sitename <- sites[1]
sitename <- 9971

sp <- st[site == sitename]
k = 4
iters = 1
param   <- as.list(coefs[k, ]); print(str(param))
par(mfrow = c(3, 1), mar = c(1.5, 2, 2, 1), mgp = c(3, 0.6, 0))
a <- rough_fitting(sitename, df, st, .FUN = wWHIT, lambda = 2, IsPlot = T, iters = iters); title("wWH2")

b1 <- rough_fitting(sitename, df, st, .FUN = wWHIT, lambda = NULL, IsPlot = T, iters = iters); title("wWHd")

b2 <- rough_fitting(sitename, df, st, .FUN = wWHIT, lambda = NULL, T, Ioptim_lambda = T, iters = iters); title("wWHd_optim")


info <- list(wWH2 = a$GOF,
             wWHd_v0 = b1$GOF,
             wWHd_latest = b2$GOF) %>%
    melt_list("meth") %>%
    # .[iter == "iter1"] %>%
    .[order(type, iter), -16]
print(info)
