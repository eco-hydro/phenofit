library(jsonlite)
library(tidyverse)
library(data.table)
library(magrittr)
library(phenofit)
library(plyr)
library(lubridate)

indir <- "D:/Documents/GoogleDrive" # /phenofit/data
files <- dir(indir, "*.geojson", full.names = T)

read_data <- function(file){
    lst <- read_json(file)$features
    data <- map(lst, function(l){
        l$properties$array %>% purrr::transpose() %>%
            map(unlist) %>%
            do.call(cbind.data.frame, .) %>%
            set_colnames(c("raw", "iter1", "iter2")) %>%
            data.table()
    })

    sites <- map_chr(lst, "id") %>% str_extract(".*(?=_)")

    df <- set_names(data, sites) %>% melt_list(var.name = "site")
    df
}

df <- fread("file:///D:/Documents/GoogleDrive/phenofit/data/fluxsites212_MOD13A1_006_0m_buffer.csv")
df <- df[, .(site = substr(`system:index`, 12, 17),
             date = ymd(date),
             EVI  = EVI/1e4,
             qc = SummaryQA)]

lst <- llply(files, read_data)
df_sm <- do.call(rbind, lst) %>% {.[order(site), ]}
df_sm$date <- df[site == "FR-LBr", date][1:409]

df_full = merge(df_sm, df, by = c("site", "date"))
df_full[is.na(qc), qc := 3]
df_full$qc %<>% as.factor()
## visualization
sites <- unique(df_full$site)

file <- "whit_GEE.pdf"
Cairo::CairoPDF(file, 10, 4)
# par(mfrow = c(4, 1), mar = c(1, 2, 3, 1), mgp = c(1.5, 0.6, 0))

## sometimes sample and reduceRegions result is different
for (i in seq_along(sites)){
    runningId(i)
    d <- df_full[site == sites[i], ]
    titlestr <- sprintf("[%03d] %s", i, sites[i])
    p <- ggplot(d, aes(date, raw, color = iters)) +
        geom_point(aes(color = qc, shape = qc)) +
        geom_line (aes(date, iter1), color = "blue", size = 0.8) +
        geom_line (aes(date, iter2), color = "red", size = 0.8) +
        geom_vline(xintercept = ymd(0101 + (2001:2017)*1e4), aes(color = NULL, shape = NULL),
                   size = 0.4, linetype=2, color = "grey60") +
        # scale_x_date(breaks = ymd(0101 + (2001:2017)*1e4))
        scale_color_manual(values = c("0" = "grey60", "1" = "#00BFC4",
                                      "2" = "#F8766D", "3" = "#C77CFF",
                                      "iter1" = "blue", "iter2" = "red"), drop = F) +
        scale_shape_manual(values = c(19, 15, 4, 17), drop = F) +
        theme_light()+
        theme(legend.position="none",
              panel.grid.minor = element_blank()) +
        ggtitle(titlestr) +
        labs(y = "EVI")
    p2 <- gridExtra::arrangeGrob(p, lgd, nrow = 2, heights = c(15, 1), padding = unit(1, "line"))

    if (i != 1) grid.newpage()
    grid.draw(p2)
}

dev.off()
file.show(file)

