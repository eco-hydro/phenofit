source("test/load_pkgs.R")
library(magrittr)
library(magick)

merge_image <- function(file_data, file_mask, outdir = "./"){
    outfile <- gsub(".jpg$", "_merge.jpg", file_data) %>% {paste0(outdir, basename(.))}

    if (file.exists(outfile)) return()
        
    p_data <- image_read(file_data)
    p_mask <- image_read(file_mask)
    
    # mask <- image_convert(image, format = NULL, type = NULL, colorspace = NULL,
    # depth = NULL, antialias = NULL)
    
    edge <- image_edge(p_mask)
    mask <- image_convert(edge, format = "rgba", colorspace = "sRGB")
    
    img <- image_fill(mask, rgb(0,0,0,0.5),fuzz = 10) %>% 
        image_transparent("black", fuzz = 10) %>%
        {image_mosaic(c(p_data, .))}
    
    image_write(img, outfile)
}

merge_roiInfo <- function(site){
    pattern <- sites_check[i]
    files   <- dir(indir, paste0(pattern, "_.*_roi.csv"), full.names = T)

    d <- llply(files, read.table, sep = ",", header = T, stringsAsFactors = F) %>% 
        do.call(rbind, .) %>% data.table()
    d$type <- str_extract(d$maskfile, "_\\w{2}_") %>% gsub("_", "", .)
    d$grp  <- str_extract(d$maskfile, ".*(?=_\\d{2})")

    x <- d[, .(
        start_date = ymd(min(start_date)),
        start_time = min(start_time),
        end_date = ymd(max(end_date)),
        end_time = max(end_time),
        maskfile = first(maskfile),
        sample_image = first(sample_image)
    ), .(grp, type)]
    # x$interval <- llply(1:nrow(x), function(i) interval(x$start_date[i],x$end_date[i]))
    x$ngrp <- x[, length(unique(type))]
    # x$interval <- interval(x$start_date, x$end_date)
    x
}

################################################################################
## 1. count csv numbers for every site
indir <- "D:/SciData/PhenoCam/PhenoCam_V1_1511/data"
files <- dir(indir, "*1day.csv")

# sites_grp  <- str_extract(files, ".*(?=_\\d{4})")
sites_type <- str_extract(files, ".*(?=_\\w{2}_)")

files_grp   <- split(files, sites_type)
sites_check <- files_grp %>% {.[sapply(., length) > 1]} %>% names()

# 2. Get mask information, and merge mask information if landcover type and mask 
# region are approximately same.
# 
# Then multiple landcover or multiple mask sites are further visually inspected. 
# Only the dominant landcover is left.
n   <- length(sites_check)

lst <- llply(sites_check, merge_roiInfo, .progress = "text") %>% set_names(sites_check)
df <- melt_list(lst, "site")
# It has been comfirmed that among 133 sites, only `silverton` with single landcover, 
# but have two masks sometime.
df <- df[ngrp != 1 | site %in% "silverton", ]

write.xlsx(df, "Phenocam_multiple_landcover_sites.xlsx")
# df$interval <- interval(df$start_date, df$end_date)
sites_multiple <- unique(df$site)


# json files --------------------------------------------------------------
## GET SITE INFO
files  <- dir(indir, "*.json", full.names = T)
info <- llply(files, function(file){
    res <-read_json(file) %>% {c(.[1], .$phenocam_site)}
    res[sapply(res, is.null)] <- NA
    res
})
info <- purrr::transpose(info) %>% lapply(unlist) %>%
    do.call(cbind.data.frame, .) %>% data.table()
info$multi <- 0

Id <- match(sites_multiple, info$sitename)
info$multi[Id] <- 1
fwrite(info, "phenocam133_site.csv")

## remove files to check
outdir <- "D:/SciData/PhenoCam/PhenoCam_V1_1511/check_multi/"
for (i in 1:length(sites_multiple)){
    runningId(i)
    j <- i

    pattern <- sites_multiple[i]

    files_cvs <- dir(indir, paste0(pattern, ".*1day.csv$"), full.names = T)
    files_data <- dir(indir, paste0(pattern, ".*jpg$"), full.names = T)
    files_mask <- gsub(".jpg$", ".tif", files_data)

    # files_cvs %>% { file.rename(., paste(dirname(dirname(.)), basename(.), sep = "/") ) }
    # retry({
        temp <- mapply(merge_image, files_data, files_mask, outdir = outdir)
    # }, maxTimes = 2)
}
# x <- sites_check[[1]]

## after manually select
info <- read.xlsx("Phenocam_multiple_landcover_sites.xlsx", 2) %>% data.table()
info_del <- info[is.na(selected), ]

files_del <- sprintf("%s/%s_1day.csv", indir, info_del$grp)
files_del %>% { file.rename(., paste(dirname(dirname(.)), basename(.), sep = "/") ) }


## write station info into shapefile and upload to GEE
# library(maptools)
df <- fread(file_st_cam)
sp <- df2sp(df)

writePointsShape(sp, "phenocam_st133.shp")


# b_overlap <- function(ints){
#     n <- length(ints)
#     flag <- FALSE
#     for (i in 1:(n-1)){
#         for (j in (i+1):n){
#             flag <- flag || int_overlaps(ints[[i]], ints[[j]])
#         }
#         # print(flag)
#         if (flag) break
#     }
#     flag
# }
