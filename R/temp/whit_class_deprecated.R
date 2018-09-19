#' whittaker smoother
wWHIT_d <- function(y, w, ylu, nptperyear, wFUN = wTSM, iters=1, lambdas=1000, df = NULL, ...,
    validation = FALSE){
    if (all(is.na(y))) return(y)
    n <- sum(w)

    OUT <- list()
    yiter <- y
    for (j in seq_along(lambdas)){
        lambda <- lambdas[j]

        fits <- list()
        for (i in 1:iters){
            if (i > 1) w <- wFUN(y, z, w, iters, nptperyear, ...)
            z_temp <- smooth_HdH(yiter, w, lambda=lambda)$z
            # If curve has been smoothed enough, it will not care about the 
            # second smooth. If no, the second one is just prepared for this
            # situation.
            z <- smooth_HdH(z_temp, w, lambda=lambda)$z #genius move
            
            if (!validation){
                z <- check_fit(z, ylu)
            }
            fits[[i]] <- z
            # yiter <- z# update y with smooth values
        }
        fits %<>% set_names(paste0('iter', 1:iters))

        # CROSS validation
        if (validation){
            h   <- fit$dhat

            df  <- sum(h)
            r   <- (y - z)/(1 - h)
            cv  <- sqrt( sum( r^2*w ) /n )
            gcv <- sqrt( sum( (r/(n-df))^2*w ))
            LV  <- whit_V(y, z, w) #L curve, D is missing now
            OUT[[j]] <- c(list(data = as_tibble(c(list(w = w), fits)),
                df = df, cv = cv, gcv = gcv), LV)
        }else{
            OUT[[j]] <- list(data = as_tibble(c(list(w = w), fits)))
        }
    }
    if (length(lambdas) == 1) OUT <- OUT[[1]]
    return(OUT)
}
#'
#' Faster whittake smoother and yhat
#' 
#' @references
#' [1]. L- and V-curves for optimal smoothing, Statistical Modelling 2015; 15(1): 91-111
smooth_HdH = function(y, w = 0 * y + 1, lambda = 1e4) {
    # Whittaker smoothing with second order differences
    # Computation of the hat diagonal (Hutchinson and de Hoog, 1986)
    # In: data vector (y), weigths (w), smoothing parameter (lambda)
    # Out: list with smooth vector (z), hat diagonal (dhat)
    # Paul Eilers, 2013
    # Prepare vectors to store system
    n     = length(y)
    g0    = rep(6, n)
    g0[1] = g0[n]      = 1
    g0[2] = g0[n - 1]  = 5
    g1    = rep(-4, n)
    g1[1] = g1[n-1]    = -2
    g1[n] = 0
    g2    = rep(1, n)
    g2[n  -1] = g2[n]  = 0
    # Store matrix G = W + lambda * D’ * D in vectors
    g0 = g0 * lambda + w
    g1 = g1 * lambda
    g2 = g2 * lambda
    # Compute U’VU decomposition (upper triangular U, diagonal V)
    u1 = u2 = v = rep(0, n)
    for (i in 1:n) {
        vi = g0[i]
        if (i > 1)
            vi = vi - v[i - 1] * u1[i - 1] ^ 2
        if (i > 2)
            vi = vi - v[i - 2] * u2[i - 2] ^ 2
        v[i] = vi
        if (i < n) {
            u = g1[i]
            if (i > 1)
                u = u - v[i - 1] * u1[i - 1] * u2[i - 1]
            u1[i] = u / vi
        }
        if (i < n - 1)
            u2[i] = g2[i] / vi
    }
    # Solve for smooth vector
    z = 0 * y
    for (i in 1:n) {
        zi = y[i] * w[i]
        if (i > 1)
            zi = zi - u1[i - 1] * z[i - 1]
        if (i > 2)
            zi = zi - u2[i - 2] * z[i - 2]
        z[i] = zi
    }
    z = z / v
    for (i in n:1) {
        zi = z[i]
        if (i < n)
            zi = zi - u1[i] * z[i + 1]
        if (i < n - 1)
            zi = zi - u2[i] * z[i + 2]
        z[i] = zi
    }
    s0 = s1 = s2 = rep(0, n)
    # Compute diagonal of inverse
    for (i in n:1) {
        i1 = i + 1
        i2 = i + 2
        s0[i] = 1 / v[i]
        if (i < n) {
            s1[i] =  - u1[i] * s0[i1]
            s0[i] = 1 / v[i] - u1[i] * s1[i]
        }
        if (i < n - 1) {
            s1[i] =  - u1[i] * s0[i1] - u2[i] * s1[i1]
            s2[i] =  - u1[i] * s1[i1] - u2[i] * s0[i2]
            s0[i] = 1 / v[i] - u1[i] * s1[i] - u2[i] * s2[i]
        }
    }
    return(list(z = z, dhat = s0))
}
smooth_HdH = compiler::cmpfun(smooth_HdH)


#' bisection
#' 
#' @references
#' [1]. https://rpubs.com/aaronsc32/bisection-method-r
#' @export
#' @examples
#' bisection(f_whitdf, 10, 1e6, toly =1e-1)
bisection <- function(f, a, b, tolx = 1e-7, toly = 1e-3, maxit = 1000, trace = T,
                      ln = TRUE, ...)
{
    # f = FUN(par, ...) #for other parameters
    cal_mid <- function(a, b, y0, y1){
        # x <- c(a, b)
        # y <- c(y0, y1)
        # if (ln){
        #     mid    <- exp( predict(lm(log(x)~y), data.frame(y = 0)) )
        # }else{
        #     lm_fit <- lm(x ~ y)
        #     mid    <- predict(lm_fit, data.frame(y = 0))
        # }
        mid <- ifelse(ln, exp( (log(a) + log(b))/2), (a + b)/2)
        return(mid)
        # ifelse(ln, exp( (log(a) + log(b))/2), (a + b)/2)
    }
    # mid  <- (a+b)/2
    # mid  <- ifelse(ln, exp( (log(a) + log(b))/2), (a + b)/2)
    y0   <- f(a)
    y1   <- f(b)
    mid  <- cal_mid(a, b, y0, y1)
    ymid <- f(mid)

    i <- 1
    if(sign(y0) == sign(y1)){
        stop("Starting vaules are not unsuitable!")
    } else {
        while(abs(ymid) > toly){
            if (i > maxit) stop('Max iteration reached!')
            # Print the value obtained in each iteration next line
            if (trace) cat(sprintf('i = %3d: lims = [%5.3f, %5.3f], mid = %.3f, fmid = %.4f\n', i, a, b, mid, ymid)) 

            if(sign(y1) == sign(ymid)){b <- mid}else{a <- mid}
            y0   <- f(a)
            y1   <- f(b)
            mid  <- cal_mid(a, b, y0, y1)
            ymid <- f(mid)

            i<-i+1
        }
        # Print the value obtained in each iteration next line
        if (trace) cat(sprintf('i = %3d: lims = [%5.3f, %5.3f], mid = %.3f, fmid = %.4f\n', i, a, b, mid, ymid)) 
    }
    return(mid)
}
