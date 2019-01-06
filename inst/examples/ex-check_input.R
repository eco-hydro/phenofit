library(phenofit)
data("MOD13A1")

dt <- tidy_MOD13.gee(MOD13A1$dt)
st <- MOD13A1$st

sitename <- dt$site[1]
d     <- dt[site == sitename, ] # get the first site data
sp    <- st[site == sitename, ] # station point
south <- sp$lat < 0
# global parameter
IsPlot = TRUE
print  = FALSE
nptperyear = 23
ypeak_min  = 0.05
wFUN = wTSM

# add one year in head and tail
dnew     <- add_HeadTail(d, south = south, nptperyear = nptperyear) 
INPUT    <- check_input(dnew$t, dnew$y, dnew$w, QC_flag = dnew$QC_flag,
     nptperyear = nptperyear, south = south, 
     maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
