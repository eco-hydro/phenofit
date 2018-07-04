source('test/stable/load_pkgs.R')
library(magrittr)

indir <- "D:/SciData/PhenoCam/PhenoCam_V1_1511/data"
files <- dir(indir, "*1day.csv")

## 1. count csv numbers for every site
indir <- "D:/SciData/PhenoCam/PhenoCam_V1_1511/data"
files <- dir(indir, "*1day.csv", full.names = T)

# sites_grp  <- str_extract(files, ".*(?=_\\d{4})")
sites_type <- str_extract(basename(files), ".*(?=_\\w{2}_)")

files_grp   <- split(files, sites_type)

# make sure land cover is consistent
files_grp_check <- files_grp %>% {.[sapply(., length) > 1]} #%>% names()
stat <- map(files_grp_check, ~str_extract(.x, "_\\w{2}_") %>% gsub("_", "", .) %>% unique()) %>%
    sapply(length) %>% unique()
if (length(stat) == 1 && stat == 1) cat("Only single landcover left for all sites.")


lst <- llply(files_grp, function(files){ llply(files, fread) %>% do.call(rbind, .) }, .progress = "text")
df  <- melt_list(lst, "site")

# define vegetation index
# a. green chromatic coordinate (GCC)
# b. vegetation contrast index (VCI)


# calculating the 90th percentile of GCC or VCI for each day (Sonnentag et al., 2012)
df_sm <- df[, .(site, date, image_count, gcc = gcc_90, vci = gcc_90/(1-gcc_90),
       snow_flag, outlierflag_gcc_90,
       smooth_gcc = smooth_gcc_90, smooth_vci = smooth_gcc_90/(1-smooth_gcc_90),
       int_flag)]

df_sm %<>% add_dn(days = c(8, 16))
df_sm <- reorder_name(df_sm)

fwrite(df_sm, "phenocam133_1day.csv")
