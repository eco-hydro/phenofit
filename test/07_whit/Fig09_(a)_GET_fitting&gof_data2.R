source('test/stable/load_pkgs.R')
source("test/07_whit/main_phenofit_test.R")

version     = "v3"
noise_percs = c(0.1, 0.3, 0.5, 0) %>% set_names(paste0("p", .*100, "%"))

################################################################################
dir_root    <- dir_flush
# dirs_raw    <- dir(dir_raw, full.names = T) %>% set_names(basename(.)) # 1th is flux

load("data_test/whit_lambda/MOD13A1_st_1e3_20180731.rda")
df$site %<>% as.character()
setkeyv(df, c("site", "t"))

get_RoughFitting_GOF <- function(file){
    d_fit <- get_Fitting(file)
    d <- merge(df, d_fit)

    by <- .(site, iters, meth) %>% names()
    by <- intersect(by, colnames(d))
    # by <- ""

    all  <- d[      , as.list(GOF_extra2(y, value)), by]
    part <- d[w == 1, as.list(GOF_extra2(y, value)), by]
    info <- listk(all, part)
    info
}

RoughtFitting_GOF.sbatch <- function(indir){
    dir_raw     <-  paste0(dir_root, indir)

    files <- dir(dir_raw, "^w", full.names = T) %>%
        list.files(recursive = T, full.names = T) %>% {
        names <- basename(dirname(.)) %>% gsub("_0%", "", .)
        set_names(., names)
    }
    file <- files[1]
    # print(files)

    temp <- par_sbatch(files, get_RoughFitting_GOF,
                       Save=T,
                       prefix = basename(dir_raw),
                       outdir = dirname(dir_raw))
}

# lst  <- llply(files, get_RoughFitting_GOF, .progress = "text")
# info <- melt_list(lst, "meth")

# indir <- "result/gee_whittaker/valid_v3_noadjparam"
#
# RoughtFitting_GOF.sbatch("result/gee_whittaker/valid_v3_adjparam")
# RoughtFitting_GOF.sbatch("result/gee_whittaker/valid_v3_noadjparam")
#
# RoughtFitting_GOF.sbatch("result/gee_whittaker/valid (20180910)")
# RoughtFitting_GOF.sbatch("result/gee_whittaker/valid (20180912)")
# RoughtFitting_GOF.sbatch("result/gee_whittaker/valid (20181024)")
# RoughtFitting_GOF.sbatch("result/gee_whittaker/valid (20181025)")
RoughtFitting_GOF.sbatch("result/gee_whittaker/valid (20180912) d714")

# indir <- "result/valid"
# RoughtFitting_GOF.sbatch(indir)

# lst <- get_sbatch(dirname(dir_raw),
#            "valid_v3_adjparam.*.RDS")
# df_info  <- melt_list(lst, "meth")
# df_info2 <- df_info[, .(site, iters, meth, NSE, R2, Rg, Rg_norm_by_obs, Rg_norm_by_pred)] %>%
#     melt(.(site, iters, meth) %>% names(), variable.name = "index") %>%
#     spread(meth, value)
