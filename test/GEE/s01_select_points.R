source('test/stable/load_pkgs.R')
source("test/GEE/pkg_seasonality.R")

file <- "file:///C:/Users/kon055/Google Drive/Github/data/MOD13A1_st_20_0m_buffer.csv"
wmin <- 0.1

dt <- tidy_gee_MOD13A1(file)
dt[, per := sum(!is.na(y) & w > wmin)/.N, site]

df <- dt[per > 0.3]
df <- df[order(IGBPcode, site, t), ]
# fwrite(df, "phenofit_seasonality_st20.csv")

## global parameters for phenofit
sites      <- unique(df$site)
nptperyear <- 23
wFUN       <- wTSM
ymax_min   <- 0.1

info <- rm_nonSIPoints(df, IsPlot = T, file = 'SI_2.pdf')

