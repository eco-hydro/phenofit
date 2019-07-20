#' Initial lambda value of Whittaker smoother
#' 
#' This function is only suitable for 16-day EVI time-series.
#' 
#' @param y Numeric vector
#' 
#' @keywords internal
#' @examples
#' library(phenofit)
#' data("MOD13A1")
#' 
#' dt <- tidy_MOD13.gee(MOD13A1$dt)
#' st <- MOD13A1$st
#' 
#' sitename <- dt$site[1]
#' d      <- dt[site == sitename, ] # get the first site data
#' lambda <- init_lambda(d$y)
#' @export
init_lambda  <- function(y){
    y        <- y[!is.na(y)] #rm NA values
    mean     <- mean(y)
    sd       <- sd(y)
    cv       <- sd/mean
    skewness <- skewness(y, type = 2)
    kurtosis <- kurtosis(y, type = 2)
    # lambda was transformed by log10
    # lambda   <- 0.555484 + 1.671514*mean - 3.434064*sd - 0.052609*skewness + 0.009057*kurtosis
    # lambda   <- 0.555465 + 1.501239*mean - 3.204295*sd - 0.031902*skewness # Just three year result
    # lambda <- 0.831120 + 1.599970*mean - 4.094027*sd - 0.035160*cv - 0.063533*skewness # all year
    lambda <- 0.9809 + 0.7247*mean - 2.6752*sd - 0.3854*skewness - 0.0604*kurtosis

    # lambda <- 0.817783 + 1.803588*mean - 4.263469*sd - 0.038240*cv - 0.066914*skewness - 0.011289*kurtosis  #3y
    return(10^lambda)
}
