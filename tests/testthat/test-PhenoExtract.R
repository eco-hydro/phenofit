context("Phenology Extraction")

# source('helper_MOD13A1.R')
wFUN = wTSM # wBisquare #

# The `maxExtendMonth` in season_3y and curvefits is different
# lambda   <- init_lambda(INPUT$y) # lambda for whittaker
brks2 <- season_3y(INPUT,
    rFUN = wWHIT, wFUN = wFUN,
    plotdat = d, IsPlot = IsPlot, IsPlot.OnlyBad = F, print = F)

param <- list(
    INPUT, brks2,
    methods = c("AG", "Zhang", "Beck", "Elmore", 'Gu'), #,"klos",
    verbose = F,
    wFUN = wFUN,
    nextent = 2, maxExtendMonth = 2, minExtendMonth = 1,
    QC_flag = dnew$QC_flag, minPercValid = 0.2,
    print = FALSE
)

fit  <- do.call(curvefits, param)

# test plot.fFIT
expect_silent({
    fit[[1]] %>% plot()   # 1th year
})

## check the curve fitting parameters
params <- get_param(fit)
print(str(params, 1))
print(params$AG)

# Test GOF
expect_silent({
    suppressWarnings({
        stat <- ldply(fit, GOF_fFITs, .id = "flag") %>% data.table()

        df_fit <- get_fitting(fit)

        g <- plot_phenofit(df_fit, brks2)
        grid::grid.newpage(); grid::grid.draw(g)
    })
})


# analytical
expect_silent({
    p <- PhenoExtract(fit,
        analytical = TRUE, smoothed.spline = FALSE,
        IsPlot = T)
})

# numDeriv
expect_silent({
    p <- PhenoExtract(fit,
        analytical = FALSE, smoothed.spline = FALSE,
        IsPlot = F)
})

# spline
expect_silent({
    p <- PhenoExtract(fit,
        analytical = FALSE, smoothed.spline = TRUE,
        IsPlot = F)
})
