library(data.table)
library(plyr)

data("MOD13A1")

dt <- tidy_MOD13.gee(MOD13A1$dt)
st <- MOD13A1$st

sitename <- dt$site[1]
d     <- dt[site == sitename, ] # get the first site data
sp    <- st[site == sitename, ] # station point
# global parameter
IsPlot = TRUE
print  = FALSE
nptperyear = 23
ypeak_min  = 0.05

dnew  <- add_HeadTail(d, nptperyear = nptperyear) # add one year in head and tail
INPUT <- check_input(dnew$t, dnew$y, dnew$w, nptperyear, south = sp$lat < 0, 
    maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)

#' data("MOD13A1")
#'
#' dt <- tidy_MOD13.gee(MOD13A1$dt)
#' st <- MOD13A1$st
#'
#' sitename <- dt$site[1]
#' d     <- dt[site == sitename, ] # get the first site data
#' sp    <- st[site == sitename, ] # station point
#' # global parameter
#' IsPlot = T
#' print  = F
#' nptperyear = 23
#' ypeak_min  = 0.05
#'
#' dnew     <- add_HeadTail(d, nptperyear = 23) # add one year in head and tail
#' INPUT    <- check_input(dnew$t, dnew$y, dnew$w, nptperyear,
#'                         maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
#' INPUT$y0 <- dnew$y   # raw time-series, for visualization
#'
