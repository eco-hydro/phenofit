library(jsonlite)
library(grid)
library(gridExtra)

source("test/load_pkgs.R")

dir_gdrive   <- "D:/Document/GoogleDrive/whit" #data/gee_phenofit/v2/

files    <- dir(dir_gdrive, "*.geojson", full.names = T)
patterns <- str_extract(basename(files), ".*(?=_\\d{4}_)") %>% unique()

df <- llply(patterns, function(pattern) readwhitMAT(dir_gdrive, pattern),
            .progress = "text") %>% set_names(patterns) %>% melt_list("meth")
df$meth %<>% as.factor()

df_flux <- df[grep("phenoflux", meth), ]
df_cam  <- df[grep("phenocam", meth), ]

## rename meth
sub_levels <- . %>% {
    level <- levels(.); mapvalues(., level, gsub("phenoflux166_|phenocam133_", "", level))
}
df_flux$meth %<>% sub_levels
df_cam$meth  %<>% sub_levels

fwrite(df_flux, "data_test/gee_whit_phenoflux166.csv")
fwrite(df_cam , "data_test/gee_whit_phenocam133.csv")
################################################################################

k = 1
if (k == 1){
    mat_whit <- df_flux
    st       <- fread(file_st_flux)
    full     <-  fread(file_flux)
} else{
    mat_whit <- df_cam
    st       <- fread(file_st_cam)
    full     <- fread(file_cam)
}
full$date %<>% ymd()

df_full = merge(full, mat_whit, by = c("site", "date"))
df_full[is.na(SummaryQA), SummaryQA := "cloud"]
df_full$SummaryQA %<>% factor(qc_levels)
## visualization
sites <- unique(df_full$site)
lgd   <- phenofit:::make_legend(linename = c("iter1", "iter2"),
                                linecolor = c("blue", "red"))

file <- "whit_GEE.pdf"
# Cairo::CairoPDF(file, 10, 4)
# par(mfrow = c(4, 1), mar = c(1, 2, 3, 1), mgp = c(1.5, 0.6, 0))
## sometimes sample and reduceRegions result is different
ps  <- list()
lwd <- 0.6

fmt_label <- function(x) sprintf("%.2f", x)
for (i in seq_along(sites)){
    runningId(i)
    d <- df_full[site == sites[i], ]
    titlestr <- sprintf("[%03d] %s", i, sites[i])

    p <- ggplot(d, aes(date, y)) +
        geom_vline(xintercept = ymd(0101 + (2001:2017)*1e4), aes(color = NULL, shape = NULL),
                   size = 0.3, linetype=2, color = "white") +
        geom_point(aes(color = SummaryQA, shape = SummaryQA), size = 1.2) +
        geom_line (aes(date, iter1), color = "blue", size = lwd) +
        geom_line (aes(date, iter2), color = "red", size = lwd) +
        geom_text(data = d[1, ], label = titlestr, x = -Inf, y =Inf, vjust = 1.5, hjust = -0.08) +
        # scale_x_date(breaks = ymd(0101 + (2001:2017)*1e4))
        scale_color_manual(values = qc_colors, drop = F) +
        scale_shape_manual(values = qc_shapes, drop = F) +
        # theme_light()+
        theme(legend.position = "none",
              panel.grid.minor = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_blank(),
              plot.margin = margin(2, 3, 2, 0))+
        scale_y_continuous(labels = fmt_label)
        # ggtitle(titlestr) +
        # labs(y = "")
    # p2 <- gridExtra::arrangeGrob(p, lgd, nrow = 2, heights = c(15, 1), padding = unit(1, "line"))
    ps[[i]] <- p
    # if (i != 1) grid.newpage()
    # grid.draw(p2)
}

# dev.off()
# file.show(file)

file = "Fig3_whit_point_example.pdf"
