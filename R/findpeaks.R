#' findpeaks
#'
#' Find peaks (maxima) in a time series. This function is modified from
#' \code{pracma::findpeaks}.
#' 
#' @param IsDiff If want to find extreme values, `IsDiff` should be true; If 
#' just want to find the continue negative or positive values, just set 
#' `IsDiff` as false.
#' @param nups minimum number of increasing steps before a peak is reached
#' @param ndowns minimum number of decreasing steps after the peak
#' @param zero can be `+', `-', or `0'; how to interprete succeeding steps
#' of the same value: increasing, decreasing, or special
#' @param peakpat define a peak as a regular pattern, such as the default
#' pattern ``[+]{1,}[-]{1,}''; if a pattern is provided, the parameters
#' \code{nups} and \code{ndowns} are not taken into account
#' @param minpeakheight The minimum (absolute) height a peak has to have
#' to be recognized as such
#' @param minpeakdistance The minimum distance (in indices) peaks have to have
#' to be counted. If the distance of two maximum extreme value less than
#' `minpeakdistance`, only the real maximum value will be left.
#' @param threshold_min Threshold is defined as the difference of peak value with
#' trough value. There are two threshold (left and right). The minimum threshold
#' should be greater than threshold_min.
#' @param threshold_max Similar as `threshold_min`, The maximum threshold should
#' be greater than `threshold_max`.
#' @param npeaks  the number of peaks to return. If \code{sortstr} = true, the 
#' largest npeaks maximum values will be returned; If \code{sortstr} = false, 
#' just the first npeaks are returned in the order of index.
#' @param sortstr Boolean, Should the peaks be returned sorted in decreasing oreder of 
#' their maximum value?
#'
#' @examples
#' x <- seq(0, 1, len = 1024)
#' pos <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81)
#' hgt <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
#' wdt <- c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005)
#' pSignal <- numeric(length(x))
#' for (i in seq(along=pos)) {
#'     pSignal <- pSignal + hgt[i]/(1 + abs((x - pos[i])/wdt[i]))^4
#' }
#'
#' plot(pSignal, type="l", col="navy"); grid()
#' x <- findpeaks(pSignal, npeaks=3, threshold_min=4, sortstr=TRUE)
#' points(val~pos, x$X, pch=20, col="maroon")
#'
#' @export
findpeaks <- function (x, IsDiff = TRUE, nups = 1, ndowns = nups, zero = "0", peakpat = NULL,
                       minpeakheight = -Inf, minpeakdistance = 1,
                       threshold_min = 0, threshold_max = 0,
                       npeaks = 0, sortstr = FALSE, IsPlot = F)
{
    stopifnot(is.vector(x, mode = "numeric") ||
                  is.vector(x, mode = "logical") || length(is.na(x)) == 0)

    if (minpeakdistance < 1)
        warning("Handling 'minpeakdistance < 1' is logically not possible.")
    if (!zero %in% c("0", "+", "-"))
        stop("Argument 'zero' can only be '0', '+', or '-'.")

    # extend the use of findpeaks:
    # 
    # 
    if (IsDiff){
        xc <- sign(diff(x))
    }else{
        xc <- x
    }
    xc <- paste(as.character(sign(xc)), collapse = "")
    xc <- gsub("1", "+", gsub("-1", "-", xc))
    if (zero != "0")      xc      <- gsub("0", zero, xc)
    if (is.null(peakpat)) peakpat <- sprintf("[+]{%d,}[-]{%d,}", nups, ndowns)

    rc <- gregexpr(peakpat, xc)[[1]]

    if (rc[1] < 0) return(NULL)
    x1 <- rc
    x2 <- rc + attr(rc, "match.length") - 1
    attributes(x1) <- NULL
    attributes(x2) <- NULL
    n <- length(x1)
    xv <- xp <- numeric(n)
    for (i in 1:n) {
        # for duplicated extreme values, get the median
        vals <- x[x1[i]:x2[i]]
        maxI <- which(vals == max(vals, na.rm = T))
        xp[i] <- floor(median(maxI)) + x1[i] - 1
        xv[i] <- x[xp[i]]
    }
    inds <- which(xv >= minpeakheight &
                      xv - pmin(x[x1], x[x2]) >= threshold_max &
                      xv - pmax(x[x1], x[x2]) >= threshold_min)
    X <- cbind(xv[inds], xp[inds], x1[inds], x2[inds])

    if (length(X) == 0) return(NULL)
    # remove near point where dist < minpeakdistance
    rm_near <- function(x){
        I <- which(diff(x[, 2]) < minpeakdistance)[1]
        if (is.na(I) | length(x) == 0){
            return(x)
        }else{
            I_del <- I + as.integer(x[I, 1] > x[I+1, 1]) #remove the min
            x <- x[-I_del, , drop = F]
            if (nrow(x) <= 1) return(x);
            rm_near(x)
        }
    }

    X <- X[order(X[, 2]), ,drop = F] # update 20180122; order according to index
    X <- rm_near(X) # sort index is necessary before `rm_near`

    if (sortstr) { # || minpeakdistance > 1
        sl <- sort.list(X[, 1], na.last = NA, decreasing = TRUE)
        X <- X[sl, , drop = FALSE]
    }

    if (npeaks > 0 && npeaks < nrow(X)) {
        X <- X[1:npeaks, , drop = FALSE]
    }

    X <- setNames(data.table(X), c("val", "pos", "left", "right"))
    if (IsPlot){
        plot(x, type = "b"); grid()
        points(val~pos, X, col = "blue")
    }
    return(list(gregexpr = rc, X = X))
}
