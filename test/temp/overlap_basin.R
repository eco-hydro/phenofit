library(rgdal)
library(maptools)

r <- readGDAL("Z:/GanRong/AUNDVI/1982-01-01.tif")
# proj4string(r) <- prj84
# r@data <- data.frame(id = 1:prod(840, 680))

basin <- readShapePoly("Z:/AU_basin/AU_basin.shp")
proj4string(basin) <- proj4string(r)

d <- over(r, basin[, 'StationID'])

I <- which(!is.na(d$StationID))
Ids <- d$StationID[I]

# 780 basin ID
Id_lst <- lapply(levels(Ids), function(x) which(Ids == x)) %>%
    set_names(levels(Ids)) %>% rm_empty()

## 2. read tif files
files <- dir("Z:/GanRong/AUNDVI/", "*.tif", full.names = T)
# month tif files,
# map, reduce
lst <- llply(files, function(file){
    r <- readGDAL(file)
    r@data$band1[I] #return
}, .progress = "text") %>%
    set_names(gsub(".tif", "", basename(files)))

## 3. aggregate info basin
df_lst <- llply(lst, function(x){
    d <- data.table(x, Ids)[x >= 0, .(NDVI = mean(x, na.rm = T)), by = Ids]
    merge(data.table(Ids = unique(Ids)), d, all.x = T)#return
})

BasinId <- df_lst[[1]]$Ids
df <- map(df_lst, ~.x$NDVI) %>% do.call(cbind, .)

# reduce(df, merge, .init = df[[1]], by = )
par(mfrow=c(4,4))
for (i in 1:16){
    plot(df[i, ]/1e4, type = "l"); grid()
}
# spplot(r[which(I == '003303'), ], sp.layout = list("sp.polygons", basin, fill = "transparent", col = "red"))
