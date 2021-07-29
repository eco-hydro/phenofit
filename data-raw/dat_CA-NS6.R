library(phenofit)
data("MOD13A1")

df <- tidy_MOD13(MOD13A1$dt)
st <- MOD13A1$st

date_start <- as.Date("2010-01-01")
date_end <- as.Date("2016-12-31")

sitename <- "CA-NS6" # df$site[1]
d <- df[site == sitename & (date >= date_start & date <= date_end), ]

CA_NS6 = d
use_data(CA_NS6, overwrite = TRUE)
