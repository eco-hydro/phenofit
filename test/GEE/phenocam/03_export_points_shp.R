library(maptools)

source('test/stable/load_pkgs.R')

df <- fread(file_st_cam)
sp <- df2sp(df)

writePointsShape(sp, "phenocam_st133.shp")
