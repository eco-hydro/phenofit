source("test/load_pkgs.R")
source("test/07_whit/main_phenofit_test.R")
# library(ggrepel)

################################################################################
# noise_percs = c(0.1, 0.3, 0.5, 0) %>% set_names(paste0("p", .*100, "%"))
# lst %<>% set_names(paste0("p",noise_percs*100, "%"))
# df_org <- melt_list(lst, "perc")

################################################################################
dir_root <- dir_flush

load("data_test/whit_lambda/MOD13A1_st_1e3_20180731.rda")
df$site %<>% as.character()
df_org <- select_valid(df, noise_perc = 0)[, 1:10]


indir <- "result/whit_lambda2" %>%  paste0(dir_root, .)
dirs_raw  <- list.dirs(indir)[-1] %>% set_names(basename(.)) # 1th is indir
files <- dir(indir, full.names = T, recursive = T)
file  <- files[1]

# res <- list()
# for (i in 1:length(dirs_raw)){
#     indir <- dirs_raw[i]
#     lst <- get_sbatch(indir)
#     # lst <- readRDS(files[1])
#     df_fit <- get_Fitting(lst)
#     info   <- get_GOF_fromFitting_I(df_fit, df_org)
#     info   <- info %$% merge(info, rough, by = "site")
#     res[[i]] <- info
# }

getGOF_whit <- function(indir){
    lst <- get_sbatch(indir)
    # lst <- readRDS(files[1])
    df_fit <- get_Fitting(lst)
    info   <- get_GOF_fromFitting_I(df_fit, df_org)
    info   <- info %$% merge(info, rough, by = "site")
    return(info)
}

outdir <- sprintf("/flush1/kon055/result/whit_lambda_GOF")
temp <- par_sbatch(dirs_raw, getGOF_whit,
                   return.res = F, Save = T, outdir = outdir)

##
# 1, 2, 3, 4, 6, 18 group and extend previous and subsequent year or not, 12 
# Whittaker lambda formula models were test here. 
# 
# Result indicate grp01_Extend and grp02_nonExtend is better than others. 
# Hence, only grp01_Extend was further compared with other methods.