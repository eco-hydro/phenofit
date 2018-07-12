# .normalize and .backnormalize function referenced by phenopix package
.normalize     <- function(x, sf) (x-sf[1])/(sf[2]-sf[1])
.backnormalize <- function(x, sf) (x+sf[1]/(sf[2]-sf[1]))*(sf[2]-sf[1])

#' check_input
#'
#' Check input data, interpolate NA values in y, remove spike values, and set
#' weights for NA in y and w.
#'
#' @template INPUT
#' @param Tn Numeric vector, night temperature, default is null. If provided,
#' Tn is used to help divide ungrowing period, and then get background value in
#' ungrowing season (see details in \code{\link[phenofit]{backval}}).
#' @param wmin Double, minimum weigth (i.e. weight of snow, ice and cloud).
#' @param missval Double, which is used to replace NA values in y. If missing,
#' the default vlaue is \code{ylu[1]}.
#' @param maxgap Integer, nptperyear/4 will be a suitable value. If continuous
#' missing value numbers less than maxgap, then interpolate those NA values by
#' zoo::na.approx; If false, then replace those NA values with a constant value
#' \code{ylu[1]}. \cr
#' Replacing NA values with a constant missing value (e.g. background value ymin)
#' is inappropriate for middle growing season points. Interpolating all values
#' by na.approx, it is unsuitable for large number continous missing segments,
#' e.g. in the start or end of growing season.
#'
#' @return A list object returned
#' \itemize{
#' \item{y} Numeric vector
#' \item{t} Numeric vector
#' \item{w} Numeric vector
#' \item{Tn} Numeric vector
#' \item{ylu} =[ymin, ymax]. \code{w_critical} is used to filter not too
#'      bad values. If the percentage good values (w=1) is greater than 30\%, then
#'      \code{w_critical}=1. The else, if the percentage of w >= 0.5 points is greater
#'      than 10\%, then \code{w_critical}=0.5. In boreal regions, even if the percentage
#'      of w >= 0.5 points is only 10\%, we still can't set \code{w_critical=wmin}.
#'      We can't rely on points with the wmin weights. Then,
#'      \code{y_good = y[w >= w_critical ]},
#'      \code{ymin = pmax( quantile(y_good, alpha/2), 0)}, \code{ymax = max(y_good)}.
#' }
#'
#' @seealso \code{\link[phenofit]{backval}}
#'
#' @export
check_input <- function(t, y, w, Tn = NULL, wmin = 0.2, missval, maxgap = 10, alpha = 0.01){
    n   <- length(y)
    if (missing(w)) w <- rep(1, n)

    # ylu   <- quantile(y[w == 1], c(alpha/2, 1 - alpha), na.rm = T) #only consider good value
    # ylu   <- range(y, na.rm = T)
    # only remove low values
    w_critical <- wmin + 0.1
    if (sum(w == 1, na.rm = T) >= n * 0.3){
        w_critical <- 1
    }else if (sum(w >= 0.5, na.rm = T) > n * 0.1){
        # Just set a small portion for boreal regions. In this way, it will give
        # more weights to marginal data.
        w_critical <- 0.5
    }
    y_good <- y[w >= w_critical & !is.na(y)]
    ylu    <- c(pmax( quantile(y_good, alpha/2), 0),
               max ( y_good, na.rm = T))
    # adjust weights according to ylu
    # if (trim){
    #     I_trim    <- y < ylu[1] #| y > ylu[2]
    #     w[I_trim] <- wmin
    # }
    # NAN values check
    if (missing(missval))
        missval <- ylu[1] #- diff(ylu)/10

    # generally, w == 0 mainly occur in winter. So it's seasonable to be assigned as minval
    y[w <= wmin]  <- missval # na is much appropriate, na.approx will replace it.
    # values out of range are setted to wmin weight.
    w[y < ylu[1] | y > ylu[2]] <- wmin
    # based on out test marginal extreme value also often occur in winter
    y[y < ylu[1]]                  <- missval
    y[y > ylu[2] & w < w_critical] <- missval

    ## 2. rm spike values
    std   <- sd(y, na.rm = T)
    ymean <- cbind(y[c(1, 1:(n - 2), n-1)], y[c(2, 3:n, n)]) %>% rowMeans(na.rm = T)
    # y[abs(y - ymean) > std]
    # which(abs(y - ymean) > std)
    I_spike <- which(abs(y - ymean) > 2*std & w < w_critical) # 95.44% interval
    y[I_spike] <- NA#missval

    ## 3. gap-fill NA values
    w[is.na(w) | is.na(y)] <- wmin
    w[w <= wmin]  <- wmin
    # left missing values were interpolated by `na.approx`
    y <- na.approx(y, maxgap = maxgap, na.rm = FALSE)
    # If still have na values after na.approx, just replace it with `missval`.
    y[is.na(y)] <- missval

    if (!is_empty(Tn)){
        Tn <- na.approx(Tn, maxgap = maxgap, na.rm = FALSE)
    }
    list(t = t, y = y, w = w, Tn = Tn, ylu = ylu)#quickly return
}

#' check_fit
#' 
#' Curve fitting values are constrained in the range of \code{ylu}.
#' 
#' @param yfit Numeric vector, curve fitting result
#' @param ylu limits of y value, [ymin, ymax]
#' @export
check_fit <- function(yfit, ylu){
    I_max <- yfit > ylu[2]
    I_min <- yfit < ylu[1]
    # yfit[I_max] <- ylu[2]
    yfit[I_min] <- ylu[1]
    return(yfit)
}

# values out of ylu, set to be na and interpolate it.
# Not export
check_fit2 <- function(y, ylu){
    I <- which(y < ylu[1] | y > ylu[2])
    if (length(I) > 0){
        n    <-length(y)
        y[I] <- NA
        y <- na.approx(y, na.rm = F)
        # if still have na values in y
        I_nona <- which(!is.na(y)) # not NA id
        if (length(I_nona) != n){
            # na values must are in tail or head now
            iBegin <- first(I_nona)
            iEnd   <- last(I_nona)
            if (iBegin > 2) y[1:iBegin] <- y[iBegin]
            if (iEnd < n)   y[iEnd:n]   <- y[iEnd]
        }
    }
    return(y)
}