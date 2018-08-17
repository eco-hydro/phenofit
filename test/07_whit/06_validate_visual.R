source('test/stable/load_pkgs.R')
source("test/07_whit/main_phenofit_test.R")
# library(ggrepel)

tidy_rough_fitting <- function(x){
    # if (is.character(file)){
    #     x <- readRDS(file)
    # } else{
    #     x <- file
    # }
    x <- rm_empty(x)
    d_rough <- map(x, "whit") %>% melt_list("site") #%>% melt_list("meth")

    d_melt <- d_rough %>% .[, .(site, t, iter1 = ziter1, iter2 = ziter2)] %>%
        melt(id.vars = c("site", "t"),
             measure.vars = c("iter1", "iter2"), variable.name = "iters") %>%
        .[, .(site, t, iters, value)]

    # list(rough = d_rough, melt = d_melt)
    d_melt
}

get_GOF <- function(d, is_valid = 1, by_site = T){
    grp <- if(by_site) .(meth, site) else .(meth)
    info <- ddply(d[I_valid == is_valid & iters == "iter1"], grp,function(d){
        with(d, GOF(y0, value))
    }) %>% data.table()
}

# newest
get_GOF2 <- function(d, is_valid = 1, by_site = T){
    # grp <- if(by_site) .(meth, site) else .(meth)
    grp <- .(site)
    info <- ddply(d[I_valid == is_valid & iters == "iter1"], grp,function(d){
        with(d, GOF(y0, value))
    }) %>% data.table()
}

# The neckness of speed is GOF
get_GOF3 <- function(df_org, df_fit){
    d <- merge(df_org[, .(site, t, y0, w0, I_valid)], df_fit, by = c("site", "t"))
    setkeyv(d, c("site", "t"))

    byname = c("site", "meth") %>% intersect(colnames(d)) # make sure by name exist

    # user  system elapsed
    # 18.87    0.03   19.09
    get_info <- function(is_valid){
        info <- d[I_valid == is_valid & iters == "iter1",
            .(info = list(GOF(y0, value))), byname]
        info$info %>% do.call(rbind, .) %>% data.table() %>% cbind(info[, ..byname], .)
    }
    list(cal = get_info(is_valid = 0), val = get_info(is_valid = 1)) %>% melt_list("type")
}

# Get GOF info for every
get_GOF_rds <- function(file, df_org){
    lst  <- readRDS(file)
    is_fine <- grepl("phenofit", basename(dirname(file))) # Is fine curve fitting?
    if (is_fine){
        # 1. read phenofit object
        df_fit <- llply(rm_empty(lst), getFittings2, .progress = "text") %>% melt_list("site")
    }else{
        # 2. read rough fitting
        df_fit <- readRDS(file) %>% rm_empty() %>% tidy_rough_fitting()
    }
    # info <- tryCatch(
    #     get_GOF3(df_org, df_fit), # GOF info
    #     error = function(e){ message(sprintf("%s : %s", file, e$message)) }
    # )
    df_fit
    # return(list(fit = df_fit, info = info)) # df_fit also input
}

# get GOF from df_fitting RDS files
get_GOF_fromFitting <- function(file, df_org){
    lst <- readRDS(file)
    info_lst <- llply(lst, function(df_fit){
        varnames <- c("site", "t", "iters", "value", "meth") %>%
            intersect(colnames(df_fit))
        df_fit <- df_fit[, ..varnames] %>% setkeyv(c("site", "t"))
        info <- suppressMessages(get_GOF3(df_org, df_fit))
        return(info)
    }, .progress = "text")
    info_lst
}

################################################################################
noise_percs = c(0.1, 0.3, 0.5) %>% set_names(paste0("p", .*100, "%"))
# lst %<>% set_names(paste0("p",noise_percs*100, "%"))
# df_org <- melt_list(lst, "perc")

################################################################################
# dir_root <- ifelse(.Platform$OS.type == "windows", "Y:/Github/phenofit_cluster/", "")
if (.Platform$OS.type == "windows"){
    dir_root <- "V:/"
} else {
    dir_root <- "/flush1/kon055/"
}

dirs_raw  <- "result/valid" %>%  paste0(dir_root, .) %>%
    {list.dirs(.)[-1]} %>% set_names(basename(.)) # 1th is indir
dirs_fit  <- "result/fitting" %>%  paste0(dir_root, .) %>%
    {list.dirs(.)[-1]} %>% set_names(basename(.)) # 1th is indir

file_info = "val_fittings.rda" #%>% paste0(dir_root, .)

load("data_test/whit_lambda/MOD13A1_st_1e3_20180731.rda")
df$site %<>% as.character()

# To save memory, evaluate GOF for every method
info <- list()
for (i in 1:3){
    runningId(i)
    ############################################################################
    noise_perc <- noise_percs[i]
    pattern    <- sprintf("%2d%%", noise_perc*100)
    # df <- fread(infile) # , strip.white = T

    df_org <- select_valid(df, noise_perc = noise_percs)[, 1:10]
    setkeyv(df_org, c("site", "t"))
    
    ## 1. save df_fit 
    issave_fitting = F
    if (issave_fitting){
        dirs = dirs_raw %>% .[grep(pattern, basename(.))] #%>% melt_list('meth')
        files <- llply(dirs, dir, pattern = "*.RDS", full.names = T) %>% unlist()
        # reoder files to balance speed
        files %<>% .[order(str_extract(basename(.), "\\d{1,}"))]

        temp  <- par_sbatch(files, get_GOF_rds, df_org = df_org, 
            Save=T, outdir=sprintf("%sresult/fitting/fitting_%s", dir_root, pattern))        
    }

    ## 2. get GOF
    iscal_gof = T
    if (iscal_gof){
        dirs  <- dirs_fit %>% .[grep(pattern, basename(.))]
        files <- llply(dirs, dir, pattern = "*.RDS", full.names = T) %>% unlist()
        temp  <- par_sbatch(files, get_GOF_fromFitting, df_org = df_org, 
            Save=T, outdir=sprintf("%sresult/val_info/info_%s", dir_root, pattern))    
    }
    # res <- list()
    # for (j in 1:length(files)){
    #     runningId(j, prefix = "j = ")
    #     res[[j]] <- get_GOF_rds(file = files[j], df_org = df_org)
    # }
    # res <- mclapply(files, get_GOF_rds, df_org = df_org, mc.cores = 16)
    # res <- llply(files, get_GOF_rds, .progress = "text")
    # list(cal = get_GOF2(d, 0), val = get_GOF2(d, 1))
    # lst[[k]] <- list(cal = cal, val = val)
    # saveRDS(res, file = sprintf("val_info_%s.RDS", pattern))
    # info[[i]] <- res
}
# info %<>% set_names(names(noise_percs))
# save(info, file = file_info)

# if (file.exists(file)){
#     load(file)
# } else{
#     ## run the following script in cluster
#     indir <- "result/valid" %>%  paste0(dir_root, .)#Y:/Github/phenofit_cluster/
#     dirs <- list.dirs(indir)[-1] %>% set_names(basename(.))
#
#     # lst <- llply(dirs[1], get_sbatch, .progress = "text")
#     # lst <- llply(dirs, function(indir){
#     #     d <- get_sbatch(indir) %>% tidy_rough_fitting()
#     # }, .progress = "text") #, .progress = "text"
#
#     lst <- mclapply(dirs, function(indir){
#         d <- get_sbatch(indir) %>% tidy_rough_fitting()
#     }, mc.cores = 4) #, .progress = "text"
#     save(lst, file = file)
# }

# load oringal data
# 
# methods <- basename(dirs_k[-1]) %>% str_extract(".*(?=_)")

# df_info <- info %>% map(function(infoI){
#     infoI %>% rm_empty() %>% set_names(methods) %>%
#         map(~melt_list(., "type")) %>% melt_list("meth")
# }) %>% melt_list("perc")
# colnames(df_info)[1] <- "site"

## Whittaker is overfitting in some extent.
# lst_tidy <- lst[[1]] %>% tidy_rough_fitting()
# d <- merge(lst_tidy$melt, df, by = c("site", "t"))

# info    <- list(cal = get_GOF(d, 0), val = get_GOF(d, 1))
# info_df <- melt_list(info, "type")[order(site)]


# dmin = 0.2
# meth <- c("wHANTS", "wSG", "wWH2")[1]# "wWH",
# pdat <- info_df %>% melt(id.vars = c("site", "meth", "type"), variable.name = "index") %>%
#     .[, name := sprintf("%s ~ %s", type, index)] %>%
#     .[index %in% c("R2", "NSE", "RMSE", "AI")] %>%
#     dcast(site+index+type+name~meth, value.var = "value")

# eval(parse(text = sprintf("pdat[, `:=`(diff = wWH - %s, kind = 0)]", meth)))
# pdat[diff >  dmin, kind := 1]
# pdat[diff < -dmin, kind := -1]

# pdat$kind %<>% factor(levels = c(-1, 0, 1), labels = c("bad", "mid", "good"))
# pdat[, table(kind), .(index)]
# # table(pdat$kind)

# p1 <- ggplot(pdat[kind != "mid",  ],
#        aes_string(meth, "wWH", color = "kind")) +
#     geom_abline(slope = 1, col = "red", size = 1) +
#     geom_point(data = pdat, aes(color = NULL), color = "black") +
#     geom_point() +
#     # geom_text_repel(aes(label = site), show.legend = F) +
#     facet_wrap(~name, scale = "free", nrow = 2) +
#     theme(aspect.ratio=1) +
#     ggtitle(basename(outfile))

# p2 <- ggplot(pdat[kind != "mid" & index == "NSE",  ],
#        aes_string(meth, "wWH", color = "kind")) +
#     geom_abline(slope = 1, col = "red", size = 1) +
#     geom_point(data = pdat, aes(color = NULL), color = "black") +
#     geom_point() +
#     geom_text_repel(aes(label = site), show.legend = F) +
#     # facet_wrap(~index, scale = "free") +
#     theme(aspect.ratio=1) +
#     coord_cartesian(xlim = c(-0.2, 1), ylim = c(-0.2, 1))

# print(p2)
# print(p1)

# meth       RMSE       NSE        R2        MAE        AI         Bias    Bias_perc         R pvalue n_sim
# 1 wHANTS 0.05392916 0.8887526 0.8951324 0.03717092 0.9690407 -0.011877982 -0.036886719 0.9461144      0 13650
# 2    wSG 0.08038025 0.7528611 0.7549066 0.05399808 0.9290617 -0.001550526 -0.004815114 0.8688536      0 13650
# 3    wWH 0.06589713 0.8338977 0.8366325 0.04496132 0.9510410 -0.005110335 -0.015869993 0.9146762      0 13650
# 4   wWH2 0.06503765 0.8382024 0.8423709 0.04347108 0.9529469 -0.008935293 -0.027748285 0.9178077      0 13650

# brks   <- season(INPUT, nptperyear,
#                FUN = whitsmw2, wFUN = wFUN, iters = 2,
#                lambda = lambda,
#                IsPlot = IsPlot, plotdat = d,
#                south = d$lat[1] < 0,
#                rymin_less = 0.6, ymax_min = ymax_min,
#                max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5) #, ...
# get growing season breaks in a 3-year moving window

# d_rough_flux <- tidy_rough_fitting(lst)
# d[I_valid == 1, GOF(y0, value), .(meth)]

# d_rough_flux <- tidy_rough_fitting("data_test/phenoflux166_rough_val.RDS")
# d_rough_cam  <- tidy_rough_fitting("data_test/phenocam133_rough_val.RDS")

# save(d_rough_cam, d_rough_flux, file = "data_test/phenofit_rough.rda")
# d2 <- readRDS("data_test/phenocam133_rough.RDS")
