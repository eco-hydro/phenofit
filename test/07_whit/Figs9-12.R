rm(list = ls())
source('test/stable/load_pkgs.R')
source("test/07_whit/main_phenofit_test.R")
source("test/07_whit/main_Whittaker_Figures.R")
source("test/07_whit/dat_flux&cam_phenofit.R")

load("data_test/whit_lambda/MOD13A1_st_1e3_20180731.rda")
methods  <- c('AG', 'BECK', 'ELMORE', 'ZHANG', 'whit_R', 'wHANTS', 'wSG', 'wWH', "wWH2")
methods2 <- c('AG', 'ZHANG', 'wHANTS', 'wSG', 'wWH', "wWH2")

st$site %<>% as.character()
st$IGBPname %<>% factor(IGBPnames_006)

indice       <- c("R2", "Bias", "RMSE", "Rg")
indice_label <- c("bold((a)*' '*R^2)", "bold((b)*' '*Bias)",
                  "bold((c)*' '*RMSE)", "bold((d)*' '*Roughness)")

###############################################################################
indir <- "V:/result/val_info/info_0_v2/"
files <- dir(indir, full.names = T, recursive = T)

df <- llply(files, tidy_info, .progress = "text") %>%
    # set_names(c("p10%", "p30%", "p50%")) %>%
    # melt_list("perc")
    do.call(rbind, .)
# df$cv %<>% unlist()

################################################################################

source('test/07_whit/Fig09_whit_GOF_spatial.R')

source('test/07_whit/Fig10_wWH vs wWH2 at 16000 sites.R')
source('test/07_whit/Fig11_wWH vs other methods at 16000 sites.R')
source('test/07_whit/Fig12_wWH vs other methods at flux&cam.R')
