#' Extract Vegetation Phenology at site scale
#'
#' @inheritParams check_input
#' @inheritParams curvefits
#' @inheritParams season
#' @param ... other parameters to [curvefits()]
#'
#' @export
process_phenofit <- function(
    y, t, w, QC_flag = NULL, nptperyear = 36,
    brks = NULL,
    rFUN = "smooth_wWHIT",
    wFUN = wTSM,
    lambda = NULL,

    TRS = c(0.1, 0.2, 0.5, 0.6, 0.8, 0.9),
    iters = 2,
    maxExtendMonth = 12, minExtendMonth = 0.5,
    minPercValid = 0,
    south = FALSE,
    verbose = TRUE,
    # lg_lambdas = seq(1, 4, 0.1),
    methods = c("AG", "Zhang", "Beck", "Elmore", "Gu"),
    outfile = NULL, # outfile for pdf figure
    title = NULL,
    ylab = "EVI", yticks = NULL,
    IsPlot.brks = FALSE,
    write.fig = TRUE, show = FALSE,
    ymin = 0.1, wmin = 0.1,
    wsnow = 0.8,
    use.y0 = FALSE,
    alpha = 0.02,
    cex = 1.5,
    overwrite = FALSE,
    .v_curve = FALSE,
    run.curvefit = TRUE,
    ...)
{
    if (!is_empty(outfile)) {
        if (file.exists(outfile) && !overwrite) return()
    } else {
        write.fig = FALSE
    }

    ## 2.1 load site data
    d_obs <- listk(t, y, w, QC_flag) %>% as.data.table()
    if (!("QC_flag" %in% colnames(d_obs))) {
        d_obs %<>% mutate(QC_flag = ifelse(w >= 0.5, "good", "cloud"))
    }

    INPUT <- check_input(d_obs$t, d_obs$y, d_obs$w, QC_flag = NULL, nptperyear,
        maxgap = ceiling(nptperyear/12*1.5),
        south = south,
        date_start = d_obs$t[1],
        date_end = last(d_obs$t))
    # frame = floor(nptperyear/8) * 2 + 1 # wSG

    if (.v_curve) {
        lg_lambdas <- seq(3.3, 5, 0.1) # 2000-
        r <- v_curve(INPUT, lg_lambdas, d = 2, IsPlot = FALSE)
        lambda = r$lambda
    }
    # print(lambda)
    # wFUN <- "wBisquare", "wTSM", threshold_max = 0.1, IGBP = CSH
    brks2  <- season_mov(INPUT,
        rFUN = get(rFUN),
        wFUN = wFUN,
        iters = iters, wmin = 0.1,
        .lambda_vcurve = TRUE,
        lambda = lambda, 
        # nf = nf, frame = frame,
        maxExtendMonth = 12, #maxExtendMonth,
        # r_max = r_max, r_min = r_min,
        # r_minPeakHeight = r_minPeakHeight,
        # calendarYear = calendarYear,
        # ...,
        # IsPlot.vc = FALSE,
        # plotdat = INPUT, print = TRUE,
        # titlestr = "")
        IsPlot = TRUE,
        IsPlot.OnlyBad = FALSE,
        # minpeakdistance = nptperyear/36*2, # 20 days
        MaxPeaksPerYear = 3,
        MaxTroughsPerYear = 4,
        # ypeak_min = 0.08,
        ...
    )
    # save(params, file = "params.rda")
    # print(list(...))

    # brks2 <- season_mov(INPUT,
    #     FUN = smooth_wWHIT, wFUN = wFUN,
    #     maxExtendMonth = 2,
    #     # minExtendMonth = minExtendMonth,
    #     wmin = wmin,
    #     .lambda_vcurve = TRUE,
    #     lambda = lambda,
    #     MaxPeaksPerYear = 3, MaxTroughsPerYear = 4,
    #     r_min = 0.1,
    #     IsPlot = TRUE, IsPlot.OnlyBad = FALSE, print = FALSE #, ...
    # )
    # IsPlot.brks

    if (!is.null(brks)) brks2$dt <- brks$dt
    # plot_season(INPUT, brks2)

    if (!run.curvefit) {
        listk(brks = brks2, data = d_obs) # , INPUT
    } else {
        dfit <- pheno <- NULL
        ## 2.4 Curve fitting
        fit <- curvefits(INPUT, brks2,
            methods = methods, # ,"klos",, 'Gu'
            wFUN = wFUN,
            iters = 2,
            wmin = wmin,
            maxExtendMonth = maxExtendMonth, minExtendMonth = minExtendMonth,
            minPercValid = minPercValid,
            print = verbose, verbose = FALSE,
            use.y0 = use.y0,
            ...
        )

        ## check the curve fitting parameters
        l_param <- get_param(fit)
        dfit    <- get_fitting(fit)
        # d_gof <- get_GOF(fit)

        ## visualization
        if (write.fig) {
            check_dir(dirname(outfile))
            # dfit2 = merge(d_obs, dfit[, -(3:4)], all.x = TRUE) # fill gaps in growing seasons
            g <- plot_curvefits(dfit, brks2, d_obs = d_obs, title = title, cex = 1.5, ylab = ylab, yticks = yticks)
            write_fig(g, outfile, 9, length(methods)*1.4, show = show)
        }
        ## 2.5 Extract phenology
        l_pheno <- get_pheno(fit, TRS = TRS, IsPlot = F) # %>% map(~melt_list(., "meth"))
        pheno <- l_pheno$doy %>% melt_list("meth")
        list(brks = brks2, data = d_obs, pheno = pheno, fit = fit, dfit = dfit)
    }
}
