context("Phenology Extraction")

source('helper_MOD13A1.R')
wFUN = wTSM

# lambda   <- init_lambda(INPUT$y) # lambda for whittaker
# # param = listk(
# #     INPUT, nptperyear,
# #     FUN = whitsmw2, wFUN = wBisquare, iters = 2,
# #     lambda,
# #     IsPlot = IsPlot, plotdat = d,
# #     south = sp$lat[1] < 0,
# #     rymin_less = 0.6, ypeak_min = ypeak_min,
# #     max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5
# # ) 
brks2 <- season_3y(INPUT, south = sp$lat[1] < 0, 
    FUN = wWHIT, wFUN = wFUN,
    plotdat = d, IsPlot = IsPlot, print = F, IsOnlyPlotbad = F)

param <- list(
    INPUT, brks2,
    methods = c("AG", "zhang", "beck", "elmore", 'Gu'), #,"klos",
    debug = F, 
    wFUN = wFUN,
    nextent = 2, maxExtendMonth = 3, minExtendMonth = 1,
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
stat  <- ldply(fit$fits, function(fits_meth){
    ldply(fits_meth, statistic.phenofit, .id = "flag")
}, .id = "meth")
fit$stat <- stat

expect_silent({
    suppressWarnings({
        g <- plot_phenofit(fit, d)
        grid::grid.newpage(); grid::grid.draw(g)    
    })
})
