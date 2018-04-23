#' prj84
#' Default Spatial object projection: WGS84.
#' @export
prj84 <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
# CRS("+init=epsg:4326")

#' extractId
#' @export
extractId <- function(station, shpFname){
    coordinates(station) <- ~long + lat
    proj4string(station) <- prj

    shp   <- readShapePoly(shpFname, proj4string = prj)
    bound <- SpatialPolygons(shp@polygons, proj4string = prj)
    ## clipped station
    clipId <- which(!is.na(over(station, bound))) %>% as.numeric

    plot(shp, axes = T)
    plot(station[clipId, ], add = T)
    clipId#return clipId
}

#' @title  construct SpatialGridDataFrame
#'
#' @param range A numeric vector, [lat_min, lat_max, lon_min, lon_max]
#' @param cellsize A numeric, grid cell size.
#' @param midgrid A vector, [midgrid_long, midgrid_lat]. If midgrid = false, then
#' begin point and end point locate on grid lines; If true, then begin point and
#' end point in the middle of grid.
#'
#' @examples
#' range <- c(25, 40, 73, 105) # Tibetan Plateau
#' Grid  <- get_grid(range, cellsize = 1/12, midgrid = TRUE)
#' @export
get_grid <- function(range, cellsize, midgrid = c(TRUE, TRUE)) {
    lat_range  <- range[1:2]
    long_range <- range[3:4]
    # lat_range <- c(25, 40); long_range <- c(73, 105)
    
    if (length(midgrid) == 1){
        midgrid <- rep(midgrid, 2)
    }else if(length(midgrid) != 2){
        midgrid <- midgrid[1:2]
        message("error: midgrid length should be 1 or 2!")
    }
    
    offset  <- c(long_range[1], lat_range[1]) + cellsize/2 * (midgrid)
    dims    <- c(diff(long_range), diff(lat_range)) / cellsize + !midgrid

    grid <- GridTopology(
        cellcentre.offset = offset,
        cellsize = c(1, 1) * cellsize, cells.dim = dims)
    #SpatialPixelsDataFrame, other than GirdDataframe. They have a big difference!
    grid <- SpatialPixelsDataFrame(grid, data = data.frame(id = seq.int(1, prod(dims))),
                                   proj4string = prj)
    return(grid)
}

#' df2sp
#' Convert dataframe data into SpatialPointsDataframe
#' @export
df2sp <- function(x, prj){  
  if (missing(prj)) prj <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
  coordinates(x) <- ~long + lat
  proj4string(x) <- prj
  return(x)
}

#' rdist.earth
#' 
#' @param x1 A position matrix [lat, lon], degree unit.
#' @param x2 A position matrix [lat, lon], degree unit.
#' @export
rdist.earth <- function (x1, x2 = NULL){
  R <- 6378.388
  # coslat1 <- cos(x1[, 2])
  # sinlat1 <- sin(x1[, 2])
  # coslon1 <- cos(x1[, 1])
  # sinlon1 <- sin(x1[, 1])
  coslat1 <- cos((x1[, 2] * pi)/180)
  sinlat1 <- sin((x1[, 2] * pi)/180)
  coslon1 <- cos((x1[, 1] * pi)/180)
  sinlon1 <- sin((x1[, 1] * pi)/180)

  coslat2 <- cos((x2[, 2] * pi)/180)
  sinlat2 <- sin((x2[, 2] * pi)/180)
  coslon2 <- cos((x2[, 1] * pi)/180)
  sinlon2 <- sin((x2[, 1] * pi)/180)
  pp <- cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1) %*%
    t(cbind(coslat2 * coslon2, coslat2 * sinlon2, sinlat2))
  return(R * acos(ifelse(abs(pp) > 1, 1 * sign(pp), pp)))
}

#' get_nearGrids
#' 
#' find nearest 3*3 grids of points
#'
#' with the help of RANN package to quickly select nearest 3*3 grids
#' use line dist to replace sphere distance
#'
#' @return nearId nearest order
#' @export
get_nearGrids <- function(grid, station, ngrid = 3*3){
  loc_grid    <- coordinates(grid)
  loc_station <- coordinates(station)
  x <- RANN::nn2(loc_grid, loc_station, ngrid)
  res <- data.frame(stationId = rep(station@data$stationId, ngrid),
                    nearId = rep(1:ngrid, rep(nrow(station), ngrid)),
                    idx = as.numeric(x$nn.idx),
                    dists = as.numeric(x$nn.dists) *6371 * pi/180)
  # rdist.earth(loc_grid[res$idx, ], loc_station)
  return(res)
}

#' gridSaveToPoly
#'
#' Write spatialPixelDataframe into polygon shpfile for the convenience of
#' checking in google earth.
#'
#' @param grid The spatialPixelDataframe or spatialGridDataframe object.
#' @param file The string of output shapefile name including path.
#' @export
gridSaveToPoly <- function(grid, file){
  poly <-  as(grid, "SpatialPolygonsDataFrame")
  # writePolyShape(poly, file) #, can't write prj info
  rgdal::writeOGR(poly, dsn = file, layer = basename(file), driver = "ESRI Shapefile")
}
