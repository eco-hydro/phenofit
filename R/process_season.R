#' divide_seasons
#'
#' @inheritParams check_input
#' @param d data.frame, with the columns of `t`, `y` and `w`.
#' @param options_season options of [season_mov()]
#' @param .v_curve If true, it will use V-curve theory to optimize Whittaker
#' parameter, lambda.
#'
#' @note site-year may be not continuous.
#'
#' @rdname season_mov
#' @export
process_season <- function(
    d_obs,
    options = list(
        # rFUN = "smooth_wWHIT",
        wFUN = "wTSM",
        # wmin = 0.1,
        # iters = 2,
        # .lambda_vcurve = TRUE, lambda = NULL,
        maxExtendMonth = 12, # maxExtendMonth,
        MaxPeaksPerYear = 3,
        MaxTroughsPerYear = 4
    ),
    nptperyear = 36, south = FALSE,
    .v_curve = FALSE,
    ...)
{
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
    # print(lambda)
    # wFUN <- "wBisquare", "wTSM", threshold_max = 0.1, IGBP = CSH
    brks2 <- season_mov(INPUT, options, ...)
    # if (!is.null(brks)) brks2$dt <- brks$dt
    # plot_season(INPUT, brks2)
    listk(INPUT, brks = brks2, data = d_obs) # , INPUT
    # listk(INPUT, brks = brks2, lambda = lambda)
}
