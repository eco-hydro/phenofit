#' Extract Vegetation Phenology at site scale
#'
#' @inheritParams check_input
#' @inheritParams curvefits
#' @inheritParams season
#' @param ... other parameters to [curvefits()]
#'
#' @export
phenofit_process <- function(y, t, w, QC_flag, nptperyear = 36,
    brks = NULL,
    wFUN = wTSM,
    fineFit = TRUE,
    TRS = c(0.1, 0.2, 0.5, 0.6, 0.8, 0.9),
    maxExtendMonth = 1.5, minExtendMonth = 0.5,
    minPercValid = 0,
    south   = FALSE,
    verbose = TRUE,
    lambda  = NULL, lg_lambdas = seq(1, 4, 0.1),
    methods = c("AG", "Zhang", "Beck", "Elmore"),
    prefix  = "",
    titlestr = NULL, ylab = NULL,
    IsPlot.brks = FALSE,
    write.fig = TRUE, show = FALSE,
    ymin = 0.1, wmin = 0.1,
    wsnow = 0.8,
    use.y0 = FALSE,
    alpha = 0.02,
    cex = 1.5,
    overwrite = FALSE,
    ...)
{
    file_pdf = sprintf("Figure/%s%s.pdf", prefix, titlestr)
    if (file.exists(file_pdf) && !overwrite) return()
    ## Parameters
    # lambda   <- 5    # non-parameter Whittaker, only suit for 16-day. Other time-scale should assign a lambda.
    # nptperyear <- 36   # How many points for a single year
    ymax_min   <- 0.1  # the maximum ymax shoud be greater than `ymax_min`
    rymin_less <- 0.8  # trough < ymin + A*rymin_less
     #wTSM #wBisquare # Weights updating function, could be one of `wTSM`, 'wBisquare', `wChen` and `wSELF`.

    ## 2.1 load site data
    # south    = FALSE
    # print      = FALSE # whether print progress

    ## 2.2 Check input data
    dnew <- data.table(y, t, date = t, w, QC_flag)
    # dnew  <- add_HeadTail(d, south, nptperyear = nptperyear) # add additional one year in head and tail
    INPUT <- check_input(dnew$t, dnew$y, dnew$w, dnew$QC_flag,
                         nptperyear, south,
                         maxgap = nptperyear/4, alpha = alpha,
                         ymin = ymin, wmin = wmin, wsnow = wsnow)

    ## 2.3 Divide growing seasons
    # if (is.null(lambda)) lambda <- v_curve(INPUT, lg_lambdas)$lambda
    brks2 <- season_mov(INPUT,
                    FUN = smooth_wWHIT, wFUN = wFUN,
                    maxExtendMonth = 3,
                    # minExtendMonth = minExtendMonth,
                    wmin = wmin,
                    IsOptim_lambda = TRUE,
                    lambda = lambda,
                    r_min = 0.1,
                    IsPlot = IsPlot.brks, IsPlot.OnlyBad = FALSE, print = FALSE, ...)
      if (!is.null(brks)) brks2$dt = brks$dt
    # }
    # plot_season(INPUT, brks2)
    dfit <- pheno <- NULL
    if (fineFit) {
        ## 2.4 Curve fitting
        fit  <- curvefits(INPUT, brks2,
                      methods = methods, #,"klos",, 'Gu'
                      wFUN = wFUN,
                      iters = 2,
                      nextend = 2,
                      wmin = wmin,
                      maxExtendMonth = maxExtendMonth, minExtendMonth = minExtendMonth,
                      minPercValid = minPercValid,
                      print = verbose, verbose = FALSE,
                      use.y0 = use.y0,
                      ...)

        ## check the curve fitting parameters
        l_param <- get_param(fit)
        dfit   <- get_fitting(fit)
        # d_gof   <- get_GOF(fit)

        ## visualization
        if (write.fig) {
            check_dir(dirname(file_pdf))

            d_obs = INPUT[c("t", "y", "QC_flag")] %>% as.data.table()
            # dfit2 = merge(d_obs, dfit[, -(3:4)], all.x = TRUE) # fill gaps in growing seasons
            g <- plot_curvefits(dfit, brks2, d_obs = d_obs, title = titlestr, cex = 1.5, title.ylab = ylab)
            write_fig(g, file_pdf, 9, 6, show = show)
        }
        ## 2.5 Extract phenology
        l_pheno <- get_pheno(fit, TRS = TRS, IsPlot = F) 
        pheno <- l_pheno$doy %>% melt_list("meth")
    }
    list(brks = brks2, fit = fit, dfit = dfit, pheno = pheno)
}

