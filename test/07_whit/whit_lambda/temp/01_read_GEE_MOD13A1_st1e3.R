source('test/stable/load_pkgs.R')

indir <- "test/GEE/data"
files <- dir(indir, '*.csv', full.names = T)
files_json <- dir(indir, '*.json$', full.names = T)

siteorder <- function(sites){
    factor(sites) %>% as.numeric()
}

read_json_gee <- function(file){
    temp <- read_json(file)$features
    # get system:id
    sites_raw <- map(temp, "id") %>% substr(12, 100) %>% siteorder()
    res <- ldply(temp, function(x){
        map(x$properties, first, default =NA) %>% unlist()
    }, .progress = "text")
    res %<>% data.table()
    res$site <- sites_raw
    res[, `:=`(IGBPcode= as.integer(CID) + 1,
               CID = NULL)]
    res %<>% reorder_name(c("site", "IGBPcode"))
    setkeyv(res, c("site", "date"))
    res#quickly return
}

res <- read_json_gee(file)
fwrite(res, "MOD13A1_st_1e3_IGBPcode17.txt")

lst <- llply(files, fread, .progress = "text")
(I_del <- which(sapply(lst, ncol) != 16))
basename(files)[I_del]

if (is_empty(I_del)){
    dt <- do.call(rbind, lst)
    siteorder <- substr(dt$`system:index`, 12, 100) %>% factor() %>% as.numeric()

    dt <- dt[, 2:(ncol(dt)-1)]
    dt$site <- siteorder
    setkeyv(dt, c("site", "date"))
    dt[, `:=`(IGBPcode = CID + 1)]

    dt$CID <- NULL
}
head(dt)

files <- dir("D:/Documents/Github/phenology/phenofit", "*.txt", full.names = T)

lst <- llply(files, fread)

lst[[2]]$site %<>% `+`(15000)
lst[[3]]$site %<>% `+`(16000)

dt <- do.call(rbind, lst)
dt <- reorder_name(dt, c("site", "IGBPcode"))
setkeyv(dt, c("site", "date"))
site_del <- dt[, all(is.na(EVI)), .(site)][V1 == TRUE, ]$site
fwrite(dt[!(site %in% site_del), ], "MOD13A1_st_1e3.csv")
# 13:Urban and Built-up Lands, 17:Water Bodies, no data


