#' skewness and kurtosis
#'
#' Inherit from package `e1071`
#' @param x a numeric vector containing the values whose skewness is to be
#' computed.
#' @param na.rm a logical value indicating whether NA values should be stripped
#' before the computation proceeds.
#' @param type an integer between 1 and 3 selecting one of the algorithms for
#' computing \code{\link[e1071]{skewness}}.
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
