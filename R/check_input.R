# .normalize and .backnormalize function referenced by phenopix package
.normalize     <- function(x, sf) (x-sf[1])/(sf[2]-sf[1])
.backnormalize <- function(x, sf) (x+sf[1]/(sf[2]-sf[1]))*(sf[2]-sf[1])

# ' @param perc_wc critical percentage of good- and marginal- quality points for
# ' `wc`.

#' check_input
#'
#' Check input data, interpolate NA values in y, remove spike values, and set
#' weights for NA in y and w.
#'
#' @param t Numeric vector, `Date` variable
#' @param y Numeric vector, vegetation index time-series
#' @param w (optional) Numeric vector, weights of `y`. If not specified,
#' weights of all `NA` values will be `wmin`, the others will be 1.0.
#' @param QC_flag Factor (optional) returned by `qcFUN`, levels should be
#' in the range of `c("snow", "cloud", "shadow", "aerosol", "marginal",
#' "good")`, others will be categoried into `others`. `QC_flag` is
#' used for visualization in [get_pheno()] and [plot_phenofit()].
#' @param nptperyear Integer, number of images per year.
#' @param south Boolean. In south hemisphere, growing year is 1 July to the
#' following year 31 June; In north hemisphere, growing year is 1 Jan to 31 Dec.
#' @param Tn Numeric vector, night temperature, default is null. If provided,
#' Tn is used to help divide ungrowing period, and then get background value in
#' ungrowing season (see details in [phenofit::backval()]).
#' @param wmin Double, minimum weight of bad points, which could be smaller
#' the weight of snow, ice and cloud.
#' @param wsnow Doulbe. Reset the weight of snow points, after get `ylu`.
#' Snow flag is an important flag of ending of growing
#' season. Snow points is more valuable than marginal points. Hence, the weight
#' of snow should be great than that of `marginal`.
#' @param missval Double, which is used to replace NA values in y. If missing,
#' the default vlaue is `ylu[1]`.
#' @param ymin If specified, `ylu[1]` is constrained greater than ymin. This
#' value is critical for bare, snow/ice land, where vegetation amplitude is quite
#' small. Generally, you can set ymin=0.08 for NDVI, ymin=0.05 for EVI,
#' ymin=0.5 gC m-2 s-1 for GPP.
#' @param maxgap Integer, nptperyear/4 will be a suitable value. If continuous
#' missing value numbers less than maxgap, then interpolate those NA values by
#' zoo::na.approx; If false, then replace those NA values with a constant value
#' `ylu[1]`. \cr
#' Replacing NA values with a constant missing value (e.g. background value ymin)
#' is inappropriate for middle growing season points. Interpolating all values
#' by na.approx, it is unsuitable for large number continous missing segments,
#' e.g. in the start or end of growing season.
#' @param alpha Double value in `[0,1]`, quantile prob of ylu_min.
#' @param alpha_high Double value in `[0,1]`, quantile prob of `ylu_max`. If not 
#' specified, `alpha_high=alpha`.
#' @param date_start,date_end starting and ending date of the original vegetation
#' time-sereis (before `add_HeadTail`)
#' @param ... Others will be ignored.
#' @param mask_spike Boolean. Whether to remove spike values?
#'
#' @return A list object returned:
#' * `t` : Numeric vector
#' * `y0`: Numeric vector, original vegetation time-series.
#' * `y` : Numeric vector, checked vegetation time-series, `NA` values are interpolated.
#' * `w` : Numeric vector
#' * `Tn`: Numeric vector
#' * `ylu`: = `[ymin, ymax]`. `w_critical` is used to filter not too bad values.
#'
#'      If the percentage good values (w=1) is greater than 30\%, then `w_critical`=1.
#'
#'      The else, if the percentage of w >= 0.5 points is greater than 10\%, then
#'      `w_critical`=0.5. In boreal regions, even if the percentage of w >= 0.5
#'      points is only 10\%, we still can't set `w_critical=wmin`.
#'
#'      We can't rely on points with the wmin weights. Then,  \cr
#'      `y_good = y[w >= w_critical]`,  \cr
#'      `ymin = pmax( quantile(y_good, alpha/2), 0)`  \cr `ymax = max(y_good)`.
#'
#' @seealso [phenofit::backval()]
#' @example inst/examples/ex-check_input.R
#' @export
check_input <- function(t, y, w, QC_flag,
    nptperyear, south = FALSE, Tn = NULL,
    wmin = 0.2,
    wsnow = 0.8,
    ymin, missval,
    maxgap, alpha = 0.02, alpha_high = NULL, 
    date_start = NULL, date_end = NULL,
    mask_spike = TRUE,
    ...)
{
    if (missing(QC_flag)) QC_flag <- NULL
    if (missing(nptperyear)){
        nptperyear <- ceiling(365/as.numeric(difftime(t[2], t[1], units = "days")))
    }
    if (missing(maxgap)) maxgap = ceiling(nptperyear/12*1.5)
    if (is.null(alpha_high)) alpha_high = alpha

    y0  <- y
    n   <- length(y)
    if (missing(w) || is.null(w)) w <- rep(1, n)

    # ylu   <- quantile(y[w == 1], c(alpha/2, 1 - alpha), na.rm = TRUE) #only consider good value
    # ylu   <- range(y, na.rm = TRUE)
    # only remove low values
    w_critical <- wmin
    if (sum(w == 1, na.rm = TRUE) >= n*0.4){
        w_critical <- 1
    }else if (sum(w >= 0.5, na.rm = TRUE) > n*0.4){
        # Just set a small portion for boreal regions. In this way, it will give
        # more weights to marginal data.
        w_critical <- 0.5
    }
    y_good <- y[w >= w_critical] %>% rm_empty()
    # alpha/2, alpha_high set to 0.05 for remote sensing (20200322)
    ylu    <- c(pmax( quantile(y_good, alpha_high/2), 0), 
               quantile(y_good, 1 - alpha/2))

    if (!missing(ymin) && !is.na(ymin)){
        # constrain back ground value
        ylu[1] <- pmax(ylu[1], ymin)
    }
    A = diff(ylu)
    # When check_ylu, ylu_max is not used. ylu_max is only used for dividing
    # growing seasons.

    # adjust weights according to ylu
    # if (trim){
    #     I_trim    <- y < ylu[1] #| y > ylu[2]
    #     w[I_trim] <- wmin
    # }
    # NAN values check
    if (missing(missval))
        missval <- ylu[1] #- diff(ylu)/10

    # generally, w == 0 mainly occur in winter. So it's seasonable to be assigned as minval
    ## 20180717 error fixed: y[w <= wmin]  <- missval # na is much appropriate, na.approx will replace it.
    # values out of range are setted to wmin weight.
    # #based on out test marginal extreme value also often occur in winter
    # #This step is really dangerous! (checked at US-Me2)
    y[y < ylu[1]] <- missval
    y[y > ylu[2] & w < pmin(w_critical + 0.01, 1)] <- missval

    w[y0 < ylu[1] | y0 > max(y_good)] <- wmin # | y > ylu[2],

    # 整治snow
    if (!is.null(QC_flag)) {
        I_snow = QC_flag == "snow" & y >= (ylu[1] + 0.2*A)
        y[I_snow] = missval
        w[QC_flag == "snow"] = wsnow
    }

    ## 2. rm spike values
    if (mask_spike) {
        # 强化除钉值模块, 20191127
        std   <- sd(y, na.rm = TRUE)
        ymov <- cbind(y[c(1, 1:(n - 2), n-1)], y[c(2, 3:n, n)]) %>% rowMeans(na.rm = TRUE)
        # ymov2 <- movmean(y, 1)
        halfwin <- ceiling(nptperyear/36) # about 10-days
        ymov2   <- movmean(y, halfwin = halfwin)
        # which(abs(y - ymean) > std) & w <= w_critical
        I_spike <- which(abs(y - ymov) > 2*std | abs(y - ymov2) > 2*std) # 95.44% interval, `(1- 2*pnorm(-2))*100`

        y[I_spike]  <- NA # missval
        y0[I_spike] <- missval # for debug
    }
    ## 3. gap-fill NA values
    w[is.na(w) | is.na(y)] <- wmin
    w[w <= wmin] <- wmin
    # left missing values were interpolated by `na.approx`
    # browser()
    y <- na.approx(y, maxgap = maxgap, na.rm = FALSE)
    # If still have na values after na.approx, just replace it with `missval`.
    y[is.na(y)] <- missval

    if (!is_empty(Tn)){
        Tn <- na.approx(Tn, maxgap = maxgap, na.rm = FALSE)
    }
    list(t = t, y0 = y0, y = y, w = w, QC_flag = QC_flag, Tn = Tn, ylu = ylu,
        nptperyear = nptperyear, south = south,
        date_start = date_start, date_end = date_end)
}
# write_fig(expression({
#     Ind = 1:length(y)
#     # Ind = t <= "2004-01-01"
#     lwd = 0.8
#     print(I_spike)
#     plot(t[Ind], y[Ind], type = "b", lwd = lwd)
#     grid()
#     lines(t, ymov, type = "b", col = "blue", lwd = lwd)
#     points(t[I_spike], y[I_spike], col = "red")
#     # points(t[I_spike2], y[I_spike2], col = "red")
#     y[I_spike]  <- NA # missval
#     y0[I_spike] <- missval # for debug
#     plot(t[Ind], y[Ind], type = "l", lwd = lwd, col = "red")
# }), "check_input.pdf", 10, 5)
# write_fig(expression({
    # plot(y0, type = "l")
    # points(I_spike, y0[I_spike])
    # lines(ymov, col = "blue")
    # lines(ymov2, col = "green")
# }), "b.pdf")

#' check_ylu
#'
#' Curve fitting values are constrained in the range of `ylu`.
#' Only constrain trough value for a stable background value. But not for peak
#' value.
#'
#' @param yfit Numeric vector, curve fitting result
#' @param ylu limits of y value, `[ymin, ymax]`
#'
#' @return yfit, the numeric vector in the range of `ylu`.
#'
#' @export
#' @examples
#' check_ylu(1:10, c(2, 8))
check_ylu <- function(yfit, ylu){
    I_max <- yfit > ylu[2]
    I_min <- yfit < ylu[1]
    # yfit[I_max] <- ylu[2]
    yfit[I_min] <- ylu[1]
    return(yfit)
}

# #' check_ylu2
# #'
# #' values out of ylu, set to be na and interpolate it.
# #' @export
# check_ylu2 <- function(y, ylu){
#     I <- which(y < ylu[1] | y > ylu[2])
#     if (length(I) > 0){
#         n    <-length(y)
#         y[I] <- NA
#         y <- na.approx(y, na.rm = F)
#         # if still have na values in y
#         I_nona <- which(!is.na(y)) # not NA id
#         if (length(I_nona) != n){
#             # na values must are in tail or head now
#             iBegin <- first(I_nona)
#             iEnd   <- last(I_nona)
#             if (iBegin > 2) y[1:iBegin] <- y[iBegin]
#             if (iEnd < n)   y[iEnd:n]   <- y[iEnd]
#         }
#     }
#     return(y)
# }
