#' weighted kurtosis
#' 
#' @references
#' [1]. Rimoldini, L., 2014. Weighted skewness and kurtosis unbiased by sample size 
#' and Gaussian uncertainties. Astron. Comput. 5, 1–8. 
#' https://doi.org/10.1016/j.ascom.2014.02.001
kurtosis <- function(x, w){
    # rm NA
    I <- !is.na(w) & !is.na(x)
    x <- x[I]
    w <- w[I]

    if (sum(w) == 0) 
        return( list(mean = NA, sd = NA, skewness = NA, kurtosis = NA) )

    V1 <- sum(w)
    V2 <- sum(w^2)
    V3 <- sum(w^3)
    V4 <- sum(w^4)

    u  <- sum( w * x) / V1
    m2 <- sum( w * (x - u)^2) / V1
    m3 <- sum( w * (x - u)^3) / V1
    m4 <- sum( w * (x - u)^4) / V1

    K2 <- M2 <- V1^2 / (V1^2 - V2) * m2
    K3 <- M3 <- V1^3 / (V1^3 - 3*V1*V2 + 2*V3) * m3
    M4_up   <- V1^2 * (V1^4 - 3*V1^2*V2 + 2*V1*V3 + 3*V2^2 - 3*V4) * m4 -3*V1^2*(2*V1^2*V2 - 2*V1*V3 - 3*V2^2 + 3*V4) * m2^2
    M4_down <- (V1^2 - V2) * (V1^4 - 6*V1^2*V2 + 8*V1*V3 + 3*V2^2 - 6*V4)
    M4 <- M4_up / M4_down
    K4 <- (V1^2 * (V1^4 - 4*V1*V3 + 3*V2^2) * m4 - 3*V1^2*(V1^4 - 2*V1^2*V2 + 4*V1*V3 - 3*V2^2) * m2^2 ) / M4_down

    skewness <- K3/K2^(3/2)
    kurtosis <- K4/K2^2
    list(mean = u, sd = m2, skewness = skewness, kurtosis = kurtosis)
}
