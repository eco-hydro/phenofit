library(phenopix)
library(tidyverse)
library(data.table)
library(plyr)
library(lubridate)

#' transform data into TIMESAT input format.
write_TSM.MOD13 <- function(d, nptperyear = 23){
    sitename <- d$site[1]
    file_y <- sprintf("TSM_%s_y.txt", sitename)
    file_w <- sprintf("TSM_%s_w.txt", sitename)
    
    write_TSM(d$EVI/1e4  , file_y, nptperyear)
    write_TSM(d$SummaryQA, file_w, nptperyear)    
}

#' Curve fitting by phenopix
phenopix_year <- function(yeari, d){
    di <- d[year == yeari, ]
    x  <- di$y
    tout <- t <- difftime(di$date, date.origin, units="days") %>% as.numeric()

    r <- FitDoubleLogBeck(x, t, tout, return.par = T)
    as.numeric(r$predicted)
}

## 1. TIMESAT ##################################################################
## Preparing INPUT data for TIMESAT
data("MOD13A1")
df <- MOD13A1$dt[date >= ymd("20010101") & date <= ymd("20171231")]
sitename <- "CA-NS6"
d  <- df[site == sitename]
# d <- tidy_MOD13.gee(d)
nptperyear <- 23
write_TSM.MOD13(d)

d <- d[, .(site, date, year = year(date), SummaryQA, y = EVI/1e4)]
# returned by TSM
d$TIMESAT <- read.table("../TSM_ex/gauss.txt") %>% as.numeric()

## 2. PHENOPIX #################################################################
years <- 2001:2017
date.origin <- ymd('2001-01-01')

l <- llply(years, phenopix_year, d = d, .progress = "text")

d$phenopix <- do.call(c, l)

## 3. show result ##############################################################

plot(y~date, d[year>=2011], type = "b")
lines(TIMESAT~date, d[year>=2011], type = "l", col = "red")
lines(phenopix~date, d[year>=2011], type = "l", col = "blue")


fwrite(d, "../phenofit_comp.txt")
