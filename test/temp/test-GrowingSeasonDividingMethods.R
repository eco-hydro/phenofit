file <- sprintf("optim_lambda, nptperyear=%d_season_bad_4.pdf", nptperyear)
CairoPDF(file, width = 10, height = 2*6)
par(mfrow = c(6, 1), mar = c(2, 3, 2, 1), mgp = c(1.5, 0.6, 0))
# dev.off()

source("R/PhenoBrks.R")
source("test/GEE/test-func.R")
source("test/GEE/pkg_seasonality.R")

library(phenofit)
library(data.table)
library(lubridate)
library(magrittr)
library(Cairo)
library(tidyverse)
library(plyr)

dt <- dt_MOD13A1[, .(site, IGBPname = IGBP, lat, long, t = date, y = EVI, SummaryQA, Tn)]
dt[, ":="(w = qc_summary(SummaryQA))]
sites      <- unique(dt$site)
sitename   <- dt$site[1]
nptperyear <- 23

################################################################################

methods <- c("sgfitw", "whitsmw2", "smooth_wHANTS")
method  <- methods[3] #"sgfitw", "whitsmw2" and "smooth_wHANTS".
FUN  <- get(method)
file <- sprintf("st10_%s.pdf", method)

CairoPDF(file, width = 10, height = 12)
par(mfrow = c(5, 1), mar = c(1, 2, 3, 1), mgp = c(1.5, 0.6, 0))

for (i in 1:length(sites)){
    runningId(i)
    sitename <- sites[i]
    # sitename <- "AU-How"
    d <- dt[site == sitename, ]
    res <- whit_brks(d, nptperyear, FUN, frame = 16, partial = F)
}
dev.off()
file.show(file)

################################################################################

# add one more year in head and tail
res   <- list()
stats <- list()
for (i in seq_along(sites)){
    runningId(i)
    sitename <- sites[i]#; grp = 1
    year <- year(d$t)
    year_begin <- year[1]
    year_end   <- last(year)

    d    <- dt[site == sitename]
    d    <- d[1:(23*3), ]
    lat  <- d$coords_x2[1] #d$lat[1]
    IGBP_code <- d$IGBPcode[1]#d$IGBP[1]
    IGBP_name <- IGBPnames[IGBP_code]

    if ( sum(d$w >= 0.5, na.rm = T) > nrow(d) * 0.3){
        INPUT_SI <- INPUT
        INPUT_SI$t %<>% {as.numeric(difftime(., ymd(20000101)))}
        # vc <- v_curve(INPUT$y, w = INPUT$w, llas = seq(-2, 2, by = 0.01), d = 2, show = F); lambda <- vc$lambda
        str_title <- sprintf("[%02d] %s, IGBP = %s, log10(lambda) = %.3f, lat = %.2f", i, sitename, IGBP_name, log10(lambda), lat)

        if (is.null(brks)) next()
        stat       <- GOF(d$y, brks$whit$iter3, INPUT$w) #, INPUT$w
        stats[[i]] <- c(site = sitename, stat)
        if (stat[['NSE']] < 0.4){
            params$IsPlot <- T
            brks <- do.call(season, params)
            title(str_title)
        }
        res[[i]]  <- listk(site = sitename, lat, IGBP_name, lambda = lambda)#, vc$
    } else{
        message(str_title)
    }
}

dev.off()
# res %<>% set_names(sites)
file.show(file)
info <- do.call(rbind, stats) %>% data.table()
