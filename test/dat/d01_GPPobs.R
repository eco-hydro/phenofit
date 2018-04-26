rm(list = ls())

reorder_name <- function (d,
                          headvars = c("site", "date", "year", "doy", "d8", "d16"),
                          tailvars = "")
{
    headvars %<>% intersect(colnames(d))
    tailvars %<>% intersect(colnames(d))
    varnames <- c(headvars, setdiff(colnames(d), union(headvars,
                                                       tailvars)), tailvars)
    if (is.data.table(d)) {
        d[, ..varnames]
    }
    else {
        d[, varnames]
    }
}
# prepare gpp input data, 04 March, 2018
source("test/stable/load_pkgs.R")
stations <- fread("F:/Github/MATLAB/PML/data/flux_166.csv")

# 1. fluxsite observations -----------------------------------------------------
df_raw <- fread("F:/Github/PML_v2/fluxsites_tidy/data/fluxsites166_official_dd.csv")
df_raw <- merge(df_raw, stations[, .(site, lat, IGBP)])
# South Hemisphere
df_raw$date %<>% ymd()

info_raw <- df_raw[GPP_NT >= 0, .N, .(site, year)]
info1 <- info_raw[N >= 365*0.6, ]
info2 <- info_raw[N <  365*0.6, ] #too less, all year data set to NA

d_obs1 <- merge(info1, df_raw, by = c("site", "year"))
# d_obs2 <- merge(info2, df_raw, by = c("site", "year"))
# d_obs2$GPP_NT <- NA
# d_obs <- rbind(d_obs1, d_obs2)

d_obs <- d_obs1
d_obs %<>% reorder_name(c("site", "IGBP", "lat","date", "year", "month", "growing", "N", "GPP_NT", "GPP_vpm"))
setkeyv(d_obs, c("site", "date"))

# d_obs <- d_obs[lat > 0, length(unique(site))] #remove south hemisphere
# x <- df[site == sitename, ]

fprintf("left %s sites...\n", length(unique(d_obs$site)))

# clip selected site
# d_vpm <- clip_selectedSite(d_vpm)

# x[, .N, .(year(date))][N>365*0.6] %>% nrow
# d_obs <- merge(d_obs, d_vpm, by = c("site", "date"), all.x = T)

save("d_obs", file = "Y:/R/phenofit/data/phenofit_INPUT_flux136_obs.rda")

rm(list = ls())

# 2. zhang yao, 2017, scientific data, VPMGPP, 8day -------------------------------
files <- dir("C:/Users/kon055/Desktop/VPMGPP", full.names = T) %>%
    set_names(gsub(".csv","", basename(.)))
d_vpm <- ldply(files, fread, .id = "site")[, c(1, 3, 4)] %>%
    set_names(c("site", "date", "GPP")) %>% as.data.table()
d_vpm[, date := as.Date(date, "%Y%j")]
d_vpm %<>% set_colnames(c("site", "date", "GPP_vpm"))
d_vpm$w <- 1

# main functions ----------------------------------------------------------
clip_selectedSite <- function(INPUT){
    d_temp <- ldply(sites, function(name){
        di_obs <- d_obs[site == name, ]
        date_begin <- min(di_obs$date)
        date_end   <- max(di_obs$date)

        INPUT[site == name & (date >= date_begin & date <= date_end), ]
    }, .id = NULL)
    merge(siteInfo, d_temp)#RETURN directly
}

load("Y:/R/phenofit/data/phenofit_INPUT_flux136_obs.rda")
load("Y:/R/phenofit/data/phenofit_MultipleINPUT_flux212.rda")

sites    <- unique(d_obs$site) %>% set_names(., .)
siteInfo <- unique(d_obs[, .(site, lat, IGBP)])

d_vpm %<>% clip_selectedSite

# Multiple Products -------------------------------------------------------
INPUT <- lst$MOD13A1
INPUT_lst <- llply(lst[1:4], clip_selectedSite, .progress = "text")
INPUT_lst %<>% c(list(GPP_vpm = d_vpm))

save(INPUT_lst, file = "Y:/R/phenofit/data/phenofit_MultipleINPUT_flux136.rda")
