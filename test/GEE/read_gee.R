# read GEE json data, clipped by buffered points
source("test/load_pkgs.R")

indir <- paste0(dir_flush, "ET&GPP/fluxnet212")
df <- read_jsons.gee(indir, pattern = ".*MOD16A2.*.geojson", IsSave = T)  # ET
df <- read_jsons.gee(indir, pattern = ".*MOD17A2H.*.geojson", IsSave = T) # GPP
df <- read_jsons.gee(indir, pattern = ".*MOD13A1.*.geojson", IsSave = T)  # 500m EVI & NDVI
df <- read_jsons.gee(indir, pattern = ".*MOD13Q1.*.geojson", IsSave = T)  # 250m EVI & NDVI
df <- read_jsons.gee(indir, pattern = ".*MCD15A3H.*.geojson", IsSave = T) # 500m, 4-day, LAI
