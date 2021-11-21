context("Phenology Extraction")

# source('helper_MOD13A1.R')
wFUN = wTSM # wBisquare #

# The `maxExtendMonth` in season_mov and curvefits is different
# lambda   <- init_lambda(INPUT$y) # lambda for whittaker
brks <- season_mov(INPUT,
    options = list(rFUN = "smooth_wWHIT", wFUN = wFUN)
)
param <- list(
    INPUT, brks,
    options = list(
        methods = c("AG", "Beck", "Elmore", "Gu", "Zhang"), # ,"klos",
        wFUN = wFUN, nextend = 2, maxExtendMonth = 3, minExtendMonth = 1,
        minPercValid = 0.2, use.rough = TRUE
    )
)

fit  <- do.call(curvefits, param)[1:6]

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
        stat  <- get_GOF(fit)
        d_fit <- get_fitting(fit)

        g <- plot_curvefits(d_fit, brks)
        grid::grid.newpage(); grid::grid.draw(g)
    })
})

test_that("get_pheno.fFIT works", {
    # analytical
    expect_silent({
        p <- get_pheno(fit,
            analytical = TRUE, smoothed.spline = FALSE, IsPlot = TRUE)
    })
    # numDeriv
    expect_silent({
        p <- get_pheno(fit,
            analytical = FALSE, smoothed.spline = FALSE, IsPlot = FALSE)
    })
    # spline
    expect_silent({
        p <- get_pheno(fit,
            analytical = FALSE, smoothed.spline = TRUE, IsPlot = FALSE)
    })
})

test_that("get_pheno.rfit works", {
    rfit = brks2rfit(brks)
    r = get_pheno(rfit)
    expect_equal(names(r), c("doy", "date"))
    expect_equal(nrow(r$doy), nrow(brks$dt))
    expect_equal(colnames(r$doy)[1:2], c("flag", "origin"))    
})
