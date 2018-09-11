library(Ipaper)
library(magrittr)
library(maptools)
library(rgdal)

## point files
file <- "F:/Github/lc_005/shp/st_1000.shp"
sp <- readOGR(file) #%>% st_sfc()
sp@data <- data.frame(site = 1:nrow(sp), IGBPcode_005deg = sp$IGBPcode)

writePointsShape(sp, "st_1e3.shp")
## raster
file_left = "F:/Github/lc_005/grid/MCD12Q1_2010_006_20_left.tif"
r <- read_stars(file_left) %>% st_as_sfc


res <- over(sp[, "IGBPcode"], r)

r <- rgdal::readGDAL(file_left)
