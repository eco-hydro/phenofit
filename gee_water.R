files <- dir("//clw-03-cdc.it.csiro.au/OSM_CBR_CoRE_working/timeseries/Climate/WATCH/", "*.nc$")
vars <- str_extract(files, "\\w{1,}(?=_\\d)")

lst <- split(files, vars)

## gz files
files_gz <- dir("//clw-03-cdc.it.csiro.au/OSM_CBR_CoRE_working/timeseries/Climate/WATCH/", "*.nc.gz$") %>%
    gsub(".gz","", .)
vars <- str_extract(files, "\\w{1,}(?=_\\d)")

lst_gz <- split(files, vars)


setdiff(files_gz, files)

a <- lst %>% .[sapply(., length) != 456]

str_extract(a$SWdown_daily_WFDEI, "\\d{4}") %>% table
