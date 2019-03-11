# load station data of MOD13A1 at flux and phenocam sites
source("test/load_pkgs.R")
# source("test/07_whit/dat_flux&cam_phenofit.R")

################################################################################
# if (.Platform$OS.type == "windows"){
#     st_flux <- fread("F:/Github/MATLAB/PML/data/flux_166.csv")
#     st_cam  <- fread("D:/SciData/PhenoCam/PhenoCam_V1_1511/phenocam133_site.csv")

#     st_flux <- st_flux[, .(ID = 1:.N, site, lat, lon = long, IGBPname = IGBP)]
#     st_cam  <- st_cam[,  .(ID = 1:.N, site = sitename, lat, lon, primary_veg_type, secondary_veg_type,
#                IGBPname = IGBPnames_006[landcover_igbp])]
#     fwrite(st_cam , file_st_cam)
#     fwrite(st_flux, file_st_flux)
# } else if (.Platform$OS.type == "unix"){
#     st_cam  <- fread(file_st_cam)
#     st_flux <- fread(file_st_flux)
# }

# # format raw MODIS VI data exported frorm GEE
# if (!file.exists(file_flux) || !file.exists(file_cam)){
#     # phenoflux
#     infile_flux  <- "file:///D:/Document/GoogleDrive/phenoflux212_MOD13A1_006_0m_buffer.csv"
#     tidy_MOD13.gee(infile_flux) %>% merge(st_flux[, .(site)]) %>% fwrite(file_flux)

#     ## phenocam
#     infile_cam  <- "file:///D:/Document/GoogleDrive/phenocam133_MOD13A1_006_0m_buffer.csv"
#     tidy_MOD13.gee(infile_cam) %>% merge(st_cam[, .(site)]) %>% fwrite(file_cam)
# }

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
