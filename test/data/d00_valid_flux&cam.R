source("test/load_pkgs.R")

## validation data for phenofit
#  Dongdong Kong, 2018-07-05

# 1. phenocam dataset -----------------------------------------------------
## 1. count csv numbers for every site
indir <- "D:/SciData/PhenoCam/PhenoCam_V1_1511/data"
files <- dir(indir, "*1day.csv", full.names = T)

# sites_grp  <- str_extract(files, ".*(?=_\\d{4})")
sites_type <- str_extract(basename(files), ".*(?=_\\w{2}_)")
files_grp  <- split(files, sites_type)

# make sure land cover is consistent
files_grp_check <- files_grp %>% {.[sapply(., length) > 1]} #%>% names()
stat <- map(files_grp_check, ~str_extract(.x, "_\\w{2}_(?=00)") %>% gsub("_", "", .) %>% unique()) %>%
    sapply(length) %>% unique()
if (length(stat) == 1 && stat == 1) cat("Only single landcover left for all sites.")

lst <- llply(files_grp, function(files){ llply(files, fread) %>% do.call(rbind, .) }, .progress = "text")
df  <- melt_list(lst, "site")

## Define vegetation index
# a. green chromatic coordinate (GCC)
# b. vegetation contrast index (VCI)

# calculating the 90th percentile of GCC or VCI for each day (Sonnentag et al., 2012)
df_sm <- df[, .(site, date, image_count, gcc = gcc_90, vci = gcc_90/(1-gcc_90),
       snow_flag, outlierflag_gcc_90,
       smooth_gcc = smooth_gcc_90, smooth_vci = smooth_gcc_90/(1-smooth_gcc_90),
       int_flag)]

df_sm %<>% add_dn(days = c(8, 16))
df_sm <- reorder_name(df_sm)

# If the midday image was not evaluated, a value of NA is assigned
# replace NA flag with 6 (new category),
df_sm[is.na(snow_flag), snow_flag := 6]
df_sm[is.na(int_flag) , int_flag  := 0]

# Only no snow or snow on ground are used (2, 3, 4).
cam_16d <- df_sm[int_flag == 0 & snow_flag %in% c(2:4, 6),
      .(gcc = median(gcc), vci = median(vci)), .(site, year, d16)]
cam_16d[, date := format(parse_date_time(paste0(year, (d16-1)*16+1), "%Y%j"))]
# fwrite(df_sm, "valid_phenocam133_1day.csv")

# setdiff(names(files_grp), df_sm_16d[, .N, .(site)]$site)

# 2. fluxnet2015 tier1 GPP data -------------------------------------------
# See fluxsites_tidy repository for detail:
# https://github.com/kongdd/fluxsites_tidy/blob/master/R/FLUXNET/fluxsite212_daily.R

flux <- fread('Z:/fluxsite212/raw/fluxsites166_daily_v20180312 (60%).csv')
flux$date %<>% ymd

flux <- flux[, .(site, date, GPP_NT, GPP_DT, LE, LE_CORR)]
flux %<>% add_dn(16)

flux_16d <- flux[, lapply(.SD, median, na.rm = T), .(site, year, d16),
                .SDcols = c("GPP_NT", "GPP_DT")]
flux_16d[, date := format(parse_date_time(paste0(year, (d16-1)*16+1), "%Y%j"))]
# flux_16d <- merge(stations, flux_16d)

## rm sites less than 1y
cam_16d[, .N, site][N>=23, ]
flux_16d[, .N, site][N>=23, ]

# cedarcreek, rosemount
## actually, 166 and 131 sites
# fwrite(cam_16d, paste0(dir_data, "valid_phenocam133_16day.csv"))
# fwrite(flux_16d, paste0(dir_data, "valid_phenoflux166_16day.csv"))
