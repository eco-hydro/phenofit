context("Phenology Extraction")

source('helper_MOD13A1.R')
wFUN =  wTSM # wBisquare #

# The `maxExtendMonth` in season_3y and curvefits is different
# lambda   <- init_lambda(INPUT$y) # lambda for whittaker
brks2 <- season_3y(INPUT, south = sp$lat[1] < 0, 
    rFUN = wWHIT, wFUN = wFUN,
    plotdat = d, IsPlot = IsPlot, print = F, IsOnlyPlotbad = F)

param <- list(
    INPUT, brks2,
    methods = c("AG", "zhang", "beck", "elmore", 'Gu'), #,"klos",
    debug = F, 
    wFUN = wFUN,
    nextent = 2, maxExtendMonth = 2, minExtendMonth = 1,
    qc = as.numeric(dnew$SummaryQA), minPercValid = 0.2,
    print = FALSE
)

fit  <- do.call(curvefits, param)

fit$INPUT   <- INPUT
fit$seasons <- brks2

## check the curve fitting parameters
params <- getparam(fit)
print(str(params, 1))
print(params$AG)

## Get GOF information

expect_silent({
    suppressWarnings({
        stat  <- ldply(fit$fits, function(fits_meth){
            ldply(fits_meth, statistic.phenofit, .id = "flag")
        }, .id = "meth")
        fit$stat <- stat

        g <- plot_phenofit(fit, d)
        grid::grid.newpage(); grid::grid.draw(g)    
    })
})
