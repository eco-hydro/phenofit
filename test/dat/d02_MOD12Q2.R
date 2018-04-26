library(jsonlite)
library(magrittr)
library(tidyverse)
library(plyr)
library(data.table)
library(lubridate)

files <- dir("C:/Users/kon055/Google Drive/Github/data/phenology/json/", ".*MCD12Q2.*json", full.names = T)

lst <- read_json(files[1])$features %>% map("properties")

len <- length(lst[[1]])
df <- map_df(lst, function(x){
    res <- x[1:(len-2)] %>% map(~first(.x)) %>% c(., x[(len-1):len])#quickly return
    res[sapply(res, is_empty)] <- NA
    as_tibble(res)
}) %>% do.call(cbind.data.frame, .) %>% as.data.table()

df[, date := ymd(date)]
df[, diff := as.integer(difftime(date, ymd("2000-01-01"), units = "days"))]

scols <- contain(df, "^Onset") %>% .[order(str_extract(., "\\d"))]

names <- c("site", "date", "flag", "Increase", "Maximum", "Decrease", "Minimum")
pheno <- df[, lapply(.SD, function(x) {x - df$diff}), ,.SDcols = scols] %>%
{set_names(., gsub("Onset_Greenness_|\\d", "", colnames(.)))} %>%
{rbind(cbind(df[, .(site, date)], .[, 1:4]) %>% .[, flag := paste0(year(date), "_1")] ,
       cbind(df[, .(site, date)], .[, 5:8]) %>% .[, flag := paste0(year(date), "_2")])} %>%
    .[order(site, flag), ..names]

# remove sites all is na
I_del <- which(rowSums(is.na(pheno[, 4:7])) == 4)
pheno <- pheno[-I_del, ]

# fwrite(pheno, "data/MCD12Q2_flux212_2001-2014.csv")
# openxlsx::write.xlsx(pheno, "data/MCD12Q2_flux212_2001-2014.xlsx")

# plot data
p_dat <-  gather(pheno, phenophase, doy, -site, -date, -flag)
p_dat$phenophase %<>% factor(c("Increase", "Maximum", "Decrease", "Minimum"))
p_dat$season <- str_extract(p_dat$flag, "(?<=_).*")
# %>% transpose()

ggplot(p_dat, aes(phenophase, doy)) +
    geom_violin(aes(fill = phenophase)) +
    geom_boxplot(alpha = 0.4, width = 0.1) +
    scale_fill_brewer() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
    facet_wrap(~season)

#     geom_dotplot(binwidth = 5, binaxis = "y")
# + geom_jitter(width = 0.2)
