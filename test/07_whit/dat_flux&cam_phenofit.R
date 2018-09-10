# load station data of MOD13A1 at flux and phenocam sites
source('test/stable/load_pkgs.R')
# source("test/07_whit/dat_flux&cam_phenofit.R")

st_flux <- fread(file_st_flux)[, .(ID, site, lat, lon, IGBPname)]
st_cam  <- fread(file_st_cam)[, .(ID, site, lat, lon, IGBPname)]

types <- c("(a) FLUXNET", "(b) PhenoCam")
st <- list(st_flux, st_cam) %>% set_names(types) %>% melt_list("type")

df_org.flux <- fread(file_flux)
df_org.cam  <- fread(file_cam)
df_org <- list(df_org.flux, df_org.cam) %>% set_names(types) %>% melt_list("type")
df_org$date %<>% ymd()
df_org$t    %<>% ymd()

## merge GEE whit result
df_org[is.na(SummaryQA), SummaryQA := "cloud"]
df_org$SummaryQA %<>% factor(qc_levels)

## visualization
sites <- unique(df_org$site)
