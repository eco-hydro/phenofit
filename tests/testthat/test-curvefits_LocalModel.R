context("curvefits")

source('helper_MOD13A1.R')

wFUN = "wTSM"
brks2 <- season_mov(INPUT,
    options = list(
        rFUN = "smooth_wWHIT", wFUN = wFUN,
        r_min = 0.05, ypeak_min = 0.05,
        lambda = 10,
        verbose = FALSE
    )
)

param <- list(
    INPUT, brks2,
    options = list(
        methods = c("AG", "Beck", "Elmore", "Gu", "Zhang"), #,"klos",
        wFUN = wFUN, nextend = 2, maxExtendMonth = 3, minExtendMonth = 1,
        minPercValid = 0.2, use.rough = TRUE, verbose = FALSE
    )
)

test_that("curvefits_LocalModel works", {
    expect_silent({
        suppressWarnings({
            # Fine fitting
            fits <- curvefits_LocalModel(
                INPUT, brks2,
                options = list(
                    methods = c("AG", "Beck", "Elmore", "Zhang", "Gu"), # ,"klos", "Gu"
                    wFUN = wFUN,
                    nextend = 2, maxExtendMonth = 2, minExtendMonth = 1, minPercValid = 0.2
                ),
                constrain = TRUE
            )
            # merge local model function into global model function
            fits_merged = merge_GlobalModels(fits)
        })
    })
})
