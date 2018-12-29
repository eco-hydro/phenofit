context("Phenology Extraction")

source('helper_MOD13A1.R')
wFUN = wTSM # wBisquare #

# The `maxExtendMonth` in season_3y and curvefits is different
# lambda   <- init_lambda(INPUT$y) # lambda for whittaker
brks2 <- season_3y(INPUT,
    rFUN = wWHIT, wFUN = wFUN,
    plotdat = d, IsPlot = IsPlot, print = F, IsOnlyPlotbad = F)

param <- list(
    INPUT, brks2,
    methods = c("AG", "Zhang", "Beck", "Elmore", 'Gu'), #,"klos",
    debug = F, 
    wFUN = wFUN,
    nextent = 2, maxExtendMonth = 2, minExtendMonth = 1,
    qc = as.numeric(dnew$SummaryQA), minPercValid = 0.2,
    print = FALSE
)

fit  <- do.call(curvefits, param)

# test plot.phenofit
expect_silent({
    fit$fits %>% first() %>% # 1th method
        first() %>% plot()   # 1th year
})

fit$INPUT   <- INPUT
fit$seasons <- brks2

## check the curve fitting parameters
params <- getparam(fit)
print(str(params, 1))
print(params$AG)

# Test GOF 
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


# analytical
expect_silent({
    p <- lapply(fit$fits, PhenoExtract, 
        analytical = TRUE, smoothed.spline = FALSE,
        IsPlot = T)
    pheno  <- map(p, tidyFitPheno, origin = INPUT$t[1]) %>% purrr::transpose()
})

# numDeriv
expect_silent({
    p <- lapply(fit$fits, PhenoExtract, 
        analytical = FALSE, smoothed.spline = FALSE,
        IsPlot = F)
    pheno  <- map(p, tidyFitPheno, origin = INPUT$t[1]) %>% purrr::transpose()
})

# spline
expect_silent({
    p <- lapply(fit$fits, PhenoExtract, 
        analytical = FALSE, smoothed.spline = TRUE,
        IsPlot = F)
    pheno  <- map(p, tidyFitPheno, origin = INPUT$t[1]) %>% purrr::transpose()
})
# print(pheno)
