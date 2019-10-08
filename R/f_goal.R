#' Goal function of fine curve fitting methods
#'
#' @inheritParams optim_pheno
#' @inheritParams Logistic
#'
#' @param fun A curve fitting function, can be one of `doubleAG`, 
#' `doubleLog.Beck`, `doubleLog.Elmore`, `doubleLog.Gu`, 
#' `doubleLog.Klos`, `doubleLog.Zhang`, see [Logistic()] 
#' for details.
#' @param ... others will be ignored.
#'
#' @return RMSE Root Mean Square Error of curve fitting values.
#' 
#' @examples
#' library(phenofit)
#' # simulate vegetation time-series
#' fFUN = doubleLog.Beck
#' par  = c(
#'     mn  = 0.1,
#'     mx  = 0.7,
#'     sos = 50,
#'     rsp = 0.1,
#'     eos = 250,
#'     rau = 0.1)
#' t    <- seq(1, 365, 8)
#' tout <- seq(1, 365, 1)
#' y <- fFUN(par, t)
#' 
#' par0 <- c(
#'     mn  = 0.15,
#'     mx  = 0.65,
#'     sos = 100,
#'     rsp = 0.12,
#'     eos = 200,
#'     rau = 0.12)
#' f_goal(par0, y, t, fFUN)
#' @export
f_goal <- function(
    par, y, t,
    fun,
    w, ylu, ...)
{
    # FUN <- match.fun(fun)
    if (!all(is.finite(par))) return(9999)

    pred <- fun(par, t = t)
    # If have no finite values, return 9999
    if (!all(is.finite(pred))) return(9999) # for Klos fitting

    if (!missing(w)) {
        if (!missing(ylu)){
            # points out of ylu should be punished!
            w[pred < ylu[1] | pred > ylu[2]] <- 0
            # pred   <- check_ylu(pred, ylu)
        }
        SSE  <- sum((y - pred)^2 * w)
    } else {
        SSE  <- sum((y - pred)^2)
    }

    # if (missing(w)) w <- rep(1, length(y))
    RMSE <- sqrt(SSE/length(y))
    # NSE  <- SSE/sum((y - mean(pred))^2)

    # 1. better handle low and high values simulation
    # xpred_2 <- sqrt(xpred_2)
    # x_2     <- sqrt(x_2)
    # xpred_2 <- log(xpred_2+1)
    # x_2     <- log(x_2+1)
    # xpred_2 <- 1/pred          # inverse NSE
    # x_2     <- 1/y

    # xpred_2 <- pred - mean(y)
    # x_2     <- y - mean(y)
    # NSE2 <- sum((x_2 - xpred_2)^2 * w)/sum((x_2 - mean(x_2))^2) #NSE

    # const <- ylu[2]
    # xpred_2 <- pred - ylu[1]; xpred_2[xpred_2 < 0] <- const
    # x_2     <- y     - ylu[1]; x_2[x_2 < 0] <- const
    return(RMSE)
}

f_goal2 <- function(
    par, y, t,
    fun,
    w = NULL, ylu = NULL, ...){
    f_goal_cpp(par, y, t, fun, w, ylu)
}
