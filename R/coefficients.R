#' skewness and kurtosis
#'
#' Inherit from package `e1071`
#' @param x a numeric vector containing the values whose skewness is to be
#' computed.
#' @param na.rm a logical value indicating whether NA values should be stripped
#' before the computation proceeds.
#' @param type an integer between 1 and 3 selecting one of the algorithms for
#' computing skewness.
#'
#' @examples
#' x = rnorm(100)
#' coef_kurtosis <- kurtosis(x)
#' coef_skewness <- skewness(x)
#' 
#' @export
kurtosis <- function (x, na.rm = FALSE, type = 3) {
    if (any(ina <- is.na(x))) {
        if (na.rm) {
            x <- x[!ina]
        } else {
            return(NA_real_)
        }
    }
    if (!(type %in% (1:3))) stop("Invalid 'type' argument.")
    n <- length(x)
    x <- x - mean(x)
    r <- n * sum(x^4)/(sum(x^2)^2)

    y <- if (type == 1) {
        r - 3
    } else if (type == 2) {
        if (n < 4) {
            warning("Need at least 4 complete observations.")
            return(NA_real_)
        }
        ((n + 1) * (r - 3) + 6) * (n - 1)/((n - 2) * (n - 3))
    } else{
        r * (1 - 1/n)^2 - 3
    }
    y
}

#' @rdname kurtosis
#' @export
skewness <- function (x, na.rm = FALSE, type = 3) {
    if (any(ina <- is.na(x))) {
        if (na.rm) {
            x <- x[!ina]
        } else {
            return(NA_real_)
        }
    }
    if (!(type %in% (1:3))) stop("Invalid 'type' argument.")
    n <- length(x)
    x <- x - mean(x)
    y <- sqrt(n) * sum(x^3)/(sum(x^2)^(3/2))

    if (type == 2) {
        if (n < 3){
            warning("Need at least 3 complete observations.")
            return(NA_real_)
        }
        y <- y * sqrt(n * (n - 1))/(n - 2)
    } else if (type == 3){
        y <- y * ((1 - 1/n))^(3/2)
    }
    y
}


#' weighted CV
#' @param x Numeric vector
#' @param w weights of different point
#' @export
#'
#' @keywords internal
#' @return Named numeric vector, (mean, sd, cv).
#' @examples
#' library(phenofit)
#' x = rnorm(100)
#' coefs <- cv_coef(x)
cv_coef <- function(x, w){
    if (missing(w)) w <- rep(1, length(x))
    if (length(x) == 0){
        return( c(mean = NA_real_, sd = NA_real_, cv = NA_real_) )
    }
    # rm NA_real_
    I <- is.finite(x)
    x <- x[I]
    w <- w[I]

    mean <- sum(x * w) / sum(w)
    sd   <- sqrt(sum((x  - mean)^2 * w) /sum(w))
    cv   <- sd / mean
    c(mean = mean, sd = sd, cv = cv) # quickly return
}


#' Critical value of determined correlation
#'
#' @param n length of observation.
#' @param NumberOfPredictor Number of predictor, including constant.
#' @param alpha significant level.
#'
#' @return `F` statistic and `R2` at significant level.
#'
#' @keywords internal
#' @references
#' Chen Yanguang (2012), Geographical Data analysis with MATLAB.
#' @examples
#' R2_critical <- R2_sign(30, NumberOfPredictor = 2, alpha = 0.05)
#' @export
R2_sign <- function(n, NumberOfPredictor = 2, alpha = 0.05){
    freedom_r = NumberOfPredictor - 1 # regression
    freedom_e = n - NumberOfPredictor # error

    F  = qf(1 - alpha, freedom_r, freedom_e)
    R2 = 1 - 1/(1 + F*freedom_r/freedom_e)

    # F = 485.1
    # F = R2/freedom_r/((1-R2)/freedom_e)
    # Rc = sqrt(/(qf(1 - alpha, 1, freedom) + freedom)) %TRUE>% print  # 0.11215
    return(list(F = F, R2 = R2))
}

#' GOF
#'
#' Good of fitting
#'
#' @param Y_obs Numeric vector, observations
#' @param Y_sim Numeric vector, corresponding simulated values
#' @param w Numeric vector, weights of every points. If w included, when
#' calculating mean, Bias, MAE, RMSE and NSE, w will be taken into considered.
#' @param include.cv If true, cv will be included.
#' @param include.r If true, r and R2 will be included.
#' 
#' @return
#' * `RMSE` root mean square error
#' * `NSE` NASH coefficient
#' * `MAE` mean absolute error
#' * `AI` Agreement index (only good points (w == 1)) participate to
#' calculate. See details in Zhang et al., (2015).
#' * `Bias` bias
#' * `Bias_perc` bias percentage
#' * `n_sim` number of valid obs
#' * `cv` Coefficient of variation
#' * `R2` correlation of determination
#' * `R` pearson correlation
#' * `pvalue` pvalue of `R`
#'
#' @references
#' Zhang Xiaoyang (2015), http://dx.doi.org/10.1016/j.rse.2014.10.012
#'
#' @examples
#' Y_obs = rnorm(100)
#' Y_sim = Y_obs + rnorm(100)/4
#' GOF(Y_obs, Y_sim)
#'
#' @export
GOF <- function(Y_obs, Y_sim, w, include.cv = FALSE, include.r = FALSE){
    if (missing(w)) w <- rep(1, length(Y_obs))

    # remove NA_real_ and Inf values in Y_sim, Y_obs and w
    valid <- function(x) !is.na(x) & is.finite(x)

    I <- which(valid(Y_sim) & valid(Y_obs) & valid(w))
    # n_obs <- length(Y_obs)
    n_sim <- length(I)

    Y_sim <- Y_sim[I]
    Y_obs <- Y_obs[I]
    w     <- w[I]

    if (include.cv) {
        CV_obs <- cv_coef(Y_obs, w)
        CV_sim <- cv_coef(Y_sim, w)
    }
    if (is_empty(Y_obs)){
        out <- c(RMSE = NA_real_, NSE = NA_real_, MAE = NA_real_, AI = NA_real_,
            Bias = NA_real_, Bias_perc = NA_real_, n_sim = NA_real_)

        if (include.r) out <- c(out, R2 = NA_real_, R = NA_real_, pvalue = NA_real_)
        if (include.cv) out <- c(out, obs = CV_obs, sim = CV_sim)
        return(out)
    }

    # R2: the portion of regression explained variance, also known as
    # coefficient of determination
    #
    # https://en.wikipedia.org/wiki/Coefficient_of_determination
    # https://en.wikipedia.org/wiki/Explained_sum_of_squares
    y_mean <- sum(Y_obs * w) / sum(w)

    SSR    <- sum( (Y_sim - y_mean)^2 * w)
    SST    <- sum( (Y_obs - y_mean)^2 * w)
    # R2     <- SSR / SST

    RE     <- Y_sim - Y_obs
    Bias   <- sum ( w*RE)     /sum(w)                     # bias
    Bias_perc <- Bias/y_mean                              # bias percentage
    MAE    <- sum ( w*abs(RE))/sum(w)                     # mean absolute error
    RMSE   <- sqrt( sum(w*(RE)^2)/sum(w) )                # root mean sqrt error

    # https://en.wikipedia.org/wiki/Nash%E2%80%93Sutcliffe_model_efficiency_coefficient
    NSE  <- 1  - sum( (RE)^2 * w) / SST # NSE coefficient

    # Observations number are not same, so comparing correlation coefficient
    # was meaningless.
    # In the current, I have no idea how to add weights `R`.
    
    
    if (include.r){
        R      <- NA_real_
        pvalue <- NA_real_
        
        tryCatch({
            cor.obj <- cor.test(Y_obs, Y_sim, use = "complete.obs")
            R       <- cor.obj$estimate[[1]]
            pvalue  <- cor.obj$p.value
        }, error = function(e){
            message(e$message)
        })

        R2 = R^2
    }    
    # In Linear regression, R2 = R^2 (R is pearson cor)
    # R2     <- summary(lm(Y_sim ~ Y_obs))$r.squared # low efficient

    # AI: Agreement Index (only good values(w==1) calculate AI)
    AI <- NA_real_
    I2 <- which(w == 1)
    if (length(I2) >= 2) {
        Y_obs = Y_obs[I2]
        Y_sim = Y_sim[I2]
        y_mean = mean(Y_obs)
        AI = 1 - sum( (Y_sim - Y_obs)^2 ) / sum( (abs(Y_sim - y_mean) + abs(Y_obs - y_mean))^2 )
    }

    out <- c(RMSE = RMSE, NSE = NSE, MAE = MAE, AI = AI,
             Bias = Bias, Bias_perc = Bias_perc, n_sim = n_sim)

    if (include.r) out <- c(out, R2 = R2, R = R, pvalue = pvalue)
    if (include.cv) out <- c(out, obs = CV_obs, sim = CV_sim)
    return(out)
}
