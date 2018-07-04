library(maptools)

st133 <- fread("file:///F:/Github/phenology/phenofit/phenocam133_site.csv")
prj <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
coordinates(st133) <- ~lon + lat
proj4string(st133) <- prj

writePointsShape(st133, "phenocam_st133.shp")
