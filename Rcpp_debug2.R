#! /usr/bin/Rscript
# library(data.table)
# library(magrittr)
# library(foreach)
# library(phenofit)
devtools::load_all()
library(data.table)
library(plyr)
data("MOD13A1")

dt <- tidy_MOD13(MOD13A1$dt)
st <- MOD13A1$st

sitename <- dt$site[1]
d <- dt[site == sitename, ] # get the first site data
sp <- st[site == sitename, ] # station point
# global parameter
IsPlot <- TRUE
print <- FALSE
nptperyear <- 23
ypeak_min <- 0.05
south <- sp$lat < 0

dnew <- add_HeadTail(d, south, nptperyear) # add one year in head and tail
INPUT <- check_input(dnew$t, dnew$y, dnew$w, dnew$QC_flag,
    nptperyear,
    south = south,
    maxgap = nptperyear / 4, alpha = 0.02, wmin = 0.2
)

# source('helper_MOD13A1.R')
wFUN <- wTSM # wBisquare #

# The `maxExtendMonth` in season_mov and curvefits is different
# lambda   <- init_lambda(INPUT$y) # lambda for whittaker
brks2 <- season_mov(INPUT,
    rFUN = smooth_wWHIT, wFUN = wFUN,
    plotdat = d, IsPlot = IsPlot, IsPlot.OnlyBad = F
)

brks2 <- season_mov(INPUT,
    rFUN = smooth_wWHIT, wFUN = wFUN,
    plotdat = d, IsPlot = IsPlot, IsPlot.OnlyBad = F
)
print(brks2)

# param <- list(
#     INPUT, brks2,
#     methods = c("AG", "Zhang", "Beck", "Elmore", "Gu"), # ,"klos",
#     wFUN = wFUN,
#     nextend = 2, maxExtendMonth = 2, minExtendMonth = 1,
#     minPercValid = 0.2
# )
