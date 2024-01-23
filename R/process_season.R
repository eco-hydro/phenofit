# ' @param .v_curve If true, it will use V-curve theory to optimize Whittaker
# ' parameter, lambda.

#' divide_seasons
#'
#' @inheritParams check_input
#' @param d_obs data.frame, with the columns of `t`, `y` and `w`.
#' @param options options of [season_mov()]
#'
#' @note site-year may be not continuous.
#' 
#' @keywords internal
#' @export
process_season <- function(
    d_obs,
    options = list(
        # rFUN = "smooth_wWHIT",
        wFUN = "wTSM",
        # wmin = 0.1,
        # iters = 2,
        # lambda = NULL,
        maxExtendMonth = 12, # maxExtendMonth,
        MaxPeaksPerYear = 3,
        MaxTroughsPerYear = 4
    ),
    nptperyear = 36, south = FALSE,
    ...)
{
    set_options(season = options, ...)
    opt = .options$season
    
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
    if (is.null(lambda)) {
        lg_lambdas <- seq(1, 5, 0.1) # 2000-
        r <- v_curve(INPUT, lg_lambdas, plot = FALSE)
        # lambda <- r$lambda
        options %<>% modifyList(r["lambda"])
    }
    # wFUN <- "wBisquare", "wTSM", threshold_max = 0.1, IGBP = CSH
    brks2 <- season_mov(INPUT, options, ...)
    # if (!is.null(brks)) brks2$dt <- brks$dt
    # plot_season(INPUT, brks2)
    listk(INPUT, brks = brks2, data = d_obs, lambda = options$lambda) # , INPUT
    # listk(INPUT, brks = brks2, lambda = lambda)
}
