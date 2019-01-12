library(phenofit)
data("MOD13A1")

df <- tidy_MOD13.gee(MOD13A1$dt)
st <- MOD13A1$st

sitename <- 'CA-NS6' # df$site[1]
d     <- df[site == sitename, ] # get the first site data
sp    <- st[site == sitename, ] # station point
south <- sp$lat < 0
# global parameter
IsPlot = TRUE
print  = FALSE
nptperyear = 23
ypeak_min  = 0.05
wFUN = wTSM

## 1. check_input
# add one year in head and tail
dnew     <- add_HeadTail(d, south = south, nptperyear = nptperyear)
INPUT    <- check_input(dnew$t, dnew$y, dnew$w, QC_flag = dnew$QC_flag,
     nptperyear = nptperyear, south = south,
     maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)

## 2. Rough fitting and growing season dividing
# Rough fitting and growing season dividing
brks2 <- season_3y(INPUT,
    rFUN = wWHIT, wFUN = wFUN,
    plotdat = d, IsPlot = IsPlot, print = FALSE, IsPlot.OnlyBad = FALSE)

## 3. Fine fitting and growing season dividing
fits <- curvefits(
    INPUT, brks2,
    methods = c("AG", "Beck", "Elmore", "Zhang"), #,"klos", "Gu"
    wFUN = wFUN,
    nextent = 2, maxExtendMonth = 2, minExtendMonth = 1, minPercValid = 0.2,
    print = TRUE, verbose = FALSE)

## 4. Phenological metric extraction
l_param   <- get_param(fits)
d_GOF     <- get_GOF(fits)
d_fitting <- get_fitting(fits)
l_pheno   <- get_pheno(fits, "AG", IsPlot=TRUE)

cairo_pdf('phenofit_MOD13A1.pdf', 12, 6)
g <- plot_phenofit(d_fitting, brks2);
grid::grid.draw(g)
dev.off()

# for shiny app: phenofit
# save(df, st, file = "flux10_MOD13A1.rda")
# shiny::runApp("inst/shiny/phenofit/", 80, TRUE)
