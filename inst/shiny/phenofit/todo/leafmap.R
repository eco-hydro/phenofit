library(leaflet)
library(Ipaper)
library(raster)
library(shiny)

# SpatialObjects
range <- c(25, 40, 73, 105) # Tibetan Plateau
grid  <- get_grid(range, cellsize = 2, midgrid = TRUE) %>% raster()
poly_grid  <- as(grid, "SpatialPolygonsDataFrame")
coord <- coordinates(poly)

# leaflet map for TP pheonlogy research
lc_colors_005 = c("#aec3d6", "#162103", "#235123", "#399b38", "#38eb38", "#39723b", 
    "#6a2424", "#c3a55f", "#b76124", "#d99125", "#92af1f", "#10104c", 
    "#cdb400", "#cc0202", "#332808", "#d7cdcc", "#f7e174", "#743411")
lc_names_005  = c('WATER', 'ENF', 'EBF', 'DNF', 'DBF', 'MF', 
    'CSH', 'OSH', 'WSA', 'SAV', 'GRA', 'WET', 
    'CRO', 'URB', 'CNV', 'SNOW', 'BSV', 'UNC')

lc_colors_006 = c("#743411", "#162103", "#235123", "#399b38", "#38eb38", "#39723b", 
    "#6a2424", "#c3a55f", "#b76124", "#d99125", "#92af1f", "#10104c", 
    "#cdb400", "#cc0202", "#332808", "#d7cdcc", "#f7e174", "#aec3d6")
lc_names_006  = c('UNC', 'ENF', 'EBF', 'DNF', 'DBF', 'MF', 
    'CSH', 'OSH', 'WSA', 'SAV', 'GRA', 'WET', 
    'CRO', 'URB', 'CNV', 'SNOW', 'BSV', 'WATER')


basemap <- function(grid, poly_grid, poly){
    leaflet() %>% 
    addTiles(
        "http://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}.jpg",
        attribution = 'Tiles &copy; Esri &mdash; Source: Esri, i-cubed, USDA, USGS, AEX, GeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP, and the GIS User Community',
        group = "World Imagery"
    ) %>%
    addWMSTiles(
        "http://webmap.ornl.gov/ogcbroker/wms?",
        layers = "10004_31", # Land cover 2007
        options = WMSTileOptions(format = "image/png", transparent = TRUE),
        attribution = "MODIS Land Cover (MCD12Q1) &copy NASA",
        group = "MODIS Land Cover"
    ) %>% 
    addLegend(colors = lc_colors_005, labels = lc_names_005, 
        group = "MODIS Land Cover", 
        position = "bottomleft", 
        title = div("MCD12Q1 2007", br(), "IGBP Land Cover", align = "center")  %>% as.character(), 
        opacity = 1) %>% 
    addProviderTiles(
        "OpenTopoMap",
        group = "Open Topo Map"
    ) %>% 
    addLayersControl(
        baseGroups = c("World Imagery", "Open Topo Map"),
        overlayGroups = c("Markers", "Polygons", "MODIS Land Cover")
    ) %>% 
    addScaleBar(position = "bottomright", options = scaleBarOptions(imperial = FALSE)) %>% 
    addPolygons(data = poly_grid, group = "Polygons",
                color = "#444444",
                fillColor = "transparent",
                weight = 1, smoothFactor = 0.5,
                # stroke = FALSE,
                opacity = 1, fillOpacity = 0.5,
                # fillColor = ~pal(log10(pop)),
                labelOptions = labelOptions(textsize = "15px"),
                highlightOptions = highlightOptions(color = "white", weight = 2,
                                                    bringToFront = TRUE),
                label = ~paste0(": ", formatC(id, big.mark = ","))) %>% 
    addPolygons(data = poly, group = "Polygons",
                color = "#444444",
                fillColor = "transparent",
                weight = 2, smoothFactor = 0.5,
                # stroke = FALSE,
                opacity = 1, fillOpacity = 0.5)
}
