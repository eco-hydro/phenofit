library(xml2)
library(magrittr)


x <- read_xml('<RasterSymbolizer> <ColorMap type="intervals" extended="false" >
              </ColorMap></RasterSymbolizer>')


colormap <- xml_child(x)

brks <- c(-Inf, 1, 2, 5, 8, 10, Inf)
ninterval <- length(brks) - 1
levels <- cut(1, brks) %>% levels()

colors <- RColorBrewer::brewer.pal(ninterval, "Spectral")

for (i in 1:ninterval){
    node <- read_xml(sprintf('<ColorMapEntry color="%s" quantity="%f" label="%s"/>', colors[i], brks[i+1], levels[i]))
    xml_add_child(colormap, node)
}

write_xml(x, "lgd.xml", options = c("no_declaration", "format"))
