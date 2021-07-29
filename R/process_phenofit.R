#' Extract Vegetation Phenology at site scale
#'
#' @param d_obs data.table with the columns of y, t, w and QC_flag (optional).
#' @inheritParams check_input
#' @inheritParams curvefits
#' @inheritParams season
#' @param ... other parameters to [curvefits()]
#'
#' @keywords internal
#' @export
process_phenofit <- function(
    d_obs,
    nptperyear = 36,
    south = FALSE,
    .v_curve = FALSE,
    options_season = list(
        # rFUN = "smooth_wWHIT",
        wFUN = "wTSM",
        # wmin = 0.1,
        # iters = 2,
        # .lambda_vcurve = TRUE, lambda = NULL,
        maxExtendMonth = 12, # maxExtendMonth,
        MaxPeaksPerYear = 3,
        MaxTroughsPerYear = 4
    ),
    options_fitting = list(
        methods = c("AG", "Zhang", "Beck", "Elmore", "Gu"),
        # methods = methods, # ,"klos",, 'Gu'
        wFUN = "wTSM",
        # iters = 2,
        # wmin = 0.1,
        maxExtendMonth = 12, minExtendMonth = 0.5,
        # minPercValid = 0,
        use.y0 = FALSE
    ),
    brks = NULL,
    TRS = c(0.1, 0.2, 0.5, 0.6, 0.8, 0.9),
    # ymin = 0.1, used for check_input
    # wsnow = 0.8,
    # use.y0 = FALSE,
    # overwrite = FALSE,
    run.curvefit = TRUE,
    ...)
{
    options_season %<>% modifyList(list(...))
    options_fitting %<>% modifyList(list(...))

    ## 2.1 load site data
    # d_obs <- listk(t, y, w, QC_flag) %>% as.data.table()
    if (!("QC_flag" %in% colnames(d_obs))) {
        d_obs %<>% mutate(QC_flag = ifelse(w >= 0.5, "good", "cloud"))
    }

    INPUT <- check_input(d_obs$t, d_obs$y, d_obs$w,
        QC_flag = d_obs$QC_flag, nptperyear,
        maxgap = ceiling(nptperyear / 12 * 1.5),
        south = south,
        date_start = d_obs$t[1],
        date_end = last(d_obs$t)
    )
    # frame = floor(nptperyear/8) * 2 + 1 # wSG
    if (.v_curve) {
        lg_lambdas <- seq(3.3, 5, 0.1) # 2000-
        r <- v_curve(INPUT, lg_lambdas, d = 2, IsPlot = FALSE)
        lambda <- r$lambda
    }
    # wFUN <- "wBisquare", "wTSM", threshold_max = 0.1, IGBP = CSH

    brks2 <- season_mov(INPUT, options_season, ...)
    # brks2 <- season_mov(INPUT,
    #     maxExtendMonth = 2,
    #     # minExtendMonth = minExtendMonth,
    #     wmin = wmin,
    #     r_min = 0.1,
    # )
    if (!is.null(brks)) brks2$dt <- brks$dt
    # plot_season(INPUT, brks2)

    if (!run.curvefit) {
        listk(brks = brks2, data = d_obs) # , INPUT
    } else {
        dfit <- pheno <- NULL
        ## 2.4 Curve fitting
        fit <- curvefits(INPUT, brks2, options_fitting, ...)

        ## check the curve fitting parameters
        # l_param <- get_param(fit)
        # d_gof <- get_GOF(fit)
        fitted.values <- get_fitting(fit)

        ## 2.5 Extract phenology
        l_pheno <- get_pheno(fit, TRS = TRS, IsPlot = F) # %>% map(~melt_list(., "meth"))
        pheno <- l_pheno$doy %>% melt_list("meth")
        listk(brks = brks2, data = d_obs, pheno, fit, fitted.values)
    }
}

## visualization
# if (write.fig) {
#     check_dir(dirname(outfile))
#     # dfit2 = merge(d_obs, dfit[, -(3:4)], all.x = TRUE) # fill gaps in growing seasons
#     g <- plot_curvefits(dfit, brks2, d_obs = d_obs, title = title, cex = 1.5, ylab = ylab, yticks = yticks)
#     write_fig(g, outfile, 9, length(methods)*1.4, show = show)
# }
