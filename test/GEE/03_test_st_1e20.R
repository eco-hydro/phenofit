siteorder <- function(sites){ factor(sites) %>% as.numeric() }

tidy_gee_MOD13A1 <- function(file){
    dt <- fread(file)
    sites_raw <- substr(dt$`system:index`, 12, 100) %>% siteorder()

    dt$site <- sites_raw
    dt[, `:=`(IGBPcode= as.integer(CID),
              CID  = NULL,
              w    = qc_summary(SummaryQA),
              date = ymd(date),
              year = year(date),
              doy  = as.integer(yday(date)),
              y    = EVI/1e4)]

    ## renew data according to DayOfYear
    dt[is.na(DayOfYear), DayOfYear := doy] #DayOfYear has missing value
    # in case of last scene of year
    dt[abs(DayOfYear - doy) >= 300, date := as.Date(sprintf("%d-%03d", year+1, DayOfYear), "%Y-%j")]
    dt[abs(DayOfYear - doy) <  300, date := as.Date(sprintf("%d-%03d", year  , DayOfYear), "%Y-%j")]

    dt %<>% reorder_name(c("site", "IGBPcode"))
    setkeyv(dt, c("site", "date"))

    return(dt)
}

dt <- tidy_gee_MOD13A1("file:///C:/Users/kon055/Google Drive/Github/MOD13A1_st_20_0m_buffer.csv")
dt[, per := sum(!is.na(EVI))/.N, site]
dt <- dt[per > 0.4, ]

colnames(dt)[14] <- "t"
dt <- dt[order(IGBPcode, site, t), ]

sites      <- unique(dt$site)
nptperyear <- 23
# dt[all(is.na(EVI)), unique(site), .(site)]
wFUN <- wTSM
ymax_min <- 0.1
