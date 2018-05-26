indir <- "C:/Users/kon055/Google Drive/01.Papers/06_Phenofit/Official/ArcGIS_prj/shp/"

x <- read.dbf("land_poly.dbf")
x <- data.table(x)
x <- x[!is.na(area)]

setkeyv(x, c("GRIDCODE", "area"))
I <- x[area > 100, ]

Ids <- dlply(x, .(GRIDCODE), function(d){
    n <- nrow(d)
    d[(n - 20 + 1):n, ]
})

poly <- rgdal::readOGR("land_poly.shp")

data <- data.table(poly@data)
I <- Ids[[1]]
a <- poly[I$ID, ]

Is <- map(Ids, "ID") %>% do.call(c, .)
a <- poly[Is, ]
a@data <- a@data[, c("GRIDCODE", "area")]
writePolyShape(a, "lc_top20.shp")
