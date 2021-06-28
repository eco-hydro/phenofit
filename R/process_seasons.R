#' divide_seasons
#'
#' @param d data.frame, with the columns of `t`, `y` and `w`.
#'
#' @param sp one line data.frame
#' @param sitename character
#' @param .movmean moving mean by `wSG`
#' @param .v_curve If true, it will use V-curve theory to optimize Whittaker
#' parameter, lambda.
#'
#' @note site-year may be not continuous.
#'
#' @export
divide_seasons <- function(d_obs, nptperyear = 23,
    iters = 2,
    south = FALSE, 
    # igbp = "GRA",
    lambda = 100,
    # nf = 3,
    rFUN = "smooth_wWHIT",
    wFUN = wTSM,
    # r_max = 0.2, r_min = 0.0,
    # r_minPeakHeight = 0.05,
    calendarYear = FALSE,
    maxExtendMonth = 12,
    # .movmean = TRUE,
    .v_curve = FALSE,
    is.plot = FALSE,
    ...)
{
    # sp = sp[1, , drop = FALSE]
    # south <- sp$lat < 0
    # maxExtendMonth <- ifelse(igbp == "EBF", 2, 2)
    ## moving average first
    # check_season(sitename, df_raw, stations)
    # dnew  <- add_HeadTail(d, south = south, nptperyear)
    # dnew = d_obs
    INPUT <- check_input(d_obs$t, d_obs$y, d_obs$w, QC_flag = NULL, nptperyear,
        maxgap = ceiling(nptperyear/12*1.5),
        south = south,
        date_start = d_obs$t[1],
        date_end = last(d_obs$t))

    # frame = floor(INPUT$nptperyear/8) * 2 + 1 # wSG
    if (.v_curve) {
        lg_lambdas <- seq(3.3, 5, 0.1) # 2000-
        r <- v_curve(INPUT, lg_lambdas, d = 2, IsPlot = FALSE)
        lambda = r$lambda
    }
    print(lambda)
    # if (.movmean) {
    #     r_wSG <- smooth_wSG(INPUT$y, INPUT$w, nptperyear = nptperyear, INPUT$ylu, iters = 2, frame = frame)
    #     INPUT$y = r_wSG$zs %>% last()
    # }
    # plot_input(INPUT)
    # browser()
    # parameters for season_mov
    # threshold_max = 0.1
    # threshold_max  <- ifelse(cv_coef(d$y)[3] >= 1, 0.1, 0.2) # empirical param
    # FUN_fit <- ifelse(sp$IGBP %in% IGBP_forest, "wHANTS", "wWHIT")
    # "wBisquare"
    
    # wFUN <- "wBisquare", "wTSM", threshold_max = 0.1, IGBP = CSH
    brks2  <- season_mov(INPUT,
        rFUN = get(rFUN),
        wFUN = wFUN,
        iters = iters, wmin = 0.1,
        .lambda_vcurve = TRUE,
        lambda = lambda, 
        # nf = nf, frame = frame,
        maxExtendMonth = maxExtendMonth,
        # r_max = r_max, r_min = r_min,
        # r_minPeakHeight = r_minPeakHeight,
        # calendarYear = calendarYear,
        # ...,
        # IsPlot.vc = FALSE,
        # plotdat = INPUT, print = TRUE,
        # titlestr = "")
        IsPlot = is.plot,
        IsPlot.OnlyBad = FALSE,
        # minpeakdistance = nptperyear/36*2, # 20 days
        MaxPeaksPerYear = 3,
        MaxTroughsPerYear = 4,
        # ypeak_min = 0.08,
        ...
    )
    print(list(...))
    listk(INPUT, brks = brks2, lambda = lambda)
}

plot_rough <- function(rFUN = "smooth_wSG", show.arrow = FALSE, ...) {

    l <- divide_seasons(dat_gpp, 365, is.plot = TRUE,
                        rFUN = rFUN,
                        .v_curve = TRUE, iters = 3, ...)
    d_rough = l$brks$fit %>% dplyr::select(t, y, starts_with("ziter")) %>%
        melt(c("t", "y"), variable.name = "ziter")

    p <- ggplot(d, aes(date, GPP)) +
        # geom_point() +
        geom_rect(data = d_ribbon, aes(x = NULL, y = NULL, xmin = xmin, xmax = xmax,
                                       group = I, fill = crop),
                  ymin = -Inf, ymax = Inf, alpha = 0.2, show.legend = F) +
        geom_line(color = "grey60") + # aes(color = cut(GPP_QC, brks_qc), group = 1)
        # geom_point(data = d[GPP_QC < 0.7], aes(shape = qc_lev, color = qc_lev), size = 1) +
        # facet_wrap(~type, nrow = 2) +
        scale_shape_manual(values = c(4, 8, 16)) +
        geom_vline(xintercept = brks_year, color = "yellow3") +
        # geom_vline(xintercept = brks_mid, color = "red", linetype = 2, size = 0.2) +
        anno_cropInfo(show.arrow = show.arrow) +
        new_scale_color() +
        geom_line(data = d_rough,aes(t, value, color = ziter)) +
        scale_color_manual(values = c("blue", "green", "red")) +
        geom_point(data = l$brks$dt, aes(beg, y_beg), color = "blue") +
        geom_point(data = l$brks$dt, aes(end, y_end), color = "blue") +
        geom_point(data = l$brks$dt, aes(peak, y_peak), color = "red") +
        labs(y = expression("GPP (gC"~m^-2~d^-1*")"))
    p
}
