#' curvefits by local model functions of TIMESAT
#' 
#' @description Local model functions `f_L(t)`, `f_C(t)` and `f_R(t)`
#' describe the VI variation in intervals around the left minima, the central
#' maxima and the right minima.
#' 
#' Local model function are merged into global model function via [merge_LocalModels()]
#' and Per J\"onsson et al. (2004; their Eq. 12),
#' where cut-off function sharply drop from 1 to 0 in small intervals around 
#' `(t_L + t_C)/2` and `(t_C + t_R)/2`.
#' 
#' \deqn{
#' F(t)= \begin{cases}
#' \alpha(t) f_{L}(t)+[1-\alpha(t)] f_{C}(t) & t_{L}<t<t_{C} \\ 
#' \beta(t) f_{C}(t)+[1-\beta(t)] f_{R}(t) & t_{C}<t<t_{R}\end{cases}
#' }
#' 
#' @inheritParams curvefits
#' @param fits List objects returned by [curvefits_LocalModel()] (not [curvefits()]).
#' 
#' @inheritSection curvefits options for fitting
#' 
#' @references 
#' 1. Per J\"onsson, P., Eklundh, L., 2004. TIMESAT - A program for analyzing
#'     time-series of satellite sensor data. Comput. Geosci. 30, 833-845.
#'     \doi{10.1016/j.cageo.2004.05.006}.
#' 
#' @seealso [curvefits()]
#' @example R/examples/ex-curvefits_LocalModel.R
#' @export
curvefits_LocalModel <- function(
    INPUT, brks,
    # methods, wFUN,
    options = list(),
    #   iters = 2, wmin = 0.1,
    #   nextend = 2, maxExtendMonth = 2, minExtendMonth = 1,
    #   minT = 0,
    #   minPercValid = 0,
    #   use.rough = FALSE,
    #   use.y0 = TRUE,
    ...)
{
    if (all(is.na(INPUT$y))) return(NULL)
    set_options(fitting = options)
    opt = .options$fitting

    QC_flag    <- INPUT$QC_flag
    nptperyear <- INPUT$nptperyear
    t          <- INPUT$t
    years <- year(t)
    n     <- length(t)

    # doys  <- as.numeric(t) # days from origin
    doys <- as.numeric(difftime(t, date.origin, units = "day")) # + 1

    # Tn for background module
    w  <- w0 <- INPUT$w
    y0 <- if (opt$use.y0) INPUT$y0 else INPUT$y
    # Tn <- INPUT$Tn # if has no Tn, NULL will be return
    # has_Tn <- ifelse(is_empty(Tn), FALSE, TRUE)

    # possible snow or cloud, replaced with Whittaker smoothing.
    I_all <- match(brks$fit$t, t) %>% rm_empty()

    if (opt$use.rough) {
        # if the range of t is smaller than `fit$t`
        INPUT$y[I_all] <- last(brks$fit)
    } else {
        I_fix <- which(w[I_all] == opt$wmin)
        I_y   <- I_all[I_fix]
        INPUT$y[I_y] <- last(brks$fit)[I_fix]
    }
    # use the weights of last iteration of rough fitting
    # w[I_all] <- brks$fit %>% {.[, contain(., "witer"), with = F]} %>% last()
    # w[I_fix] <- wmin + 0.1 # exert the function of whitaker smoother

    # growing season division
    di <- data.table( beg  = getDateId_before(brks$dt$beg, t),
                      peak = getDateId_before(brks$dt$peak, t),
                      end  = getDateId_after(brks$dt$end, t)) #%>% na.omit()

    # TODO: fix the bug when seasons are not continuous
    s = c(di$beg, di$peak, di$end %>% last()) %>% sort()
    nlocal = length(s) - 2

    width_ylu = nptperyear*2

    y    <- INPUT$y
    fits <- vector(nrow(di), mode = "list")

    curvefit_gs <- function(I_beg, I_end, type = 1L) {
        I     <- I_beg:I_end
        # I_extend <- get_extentI(w0, I_beg, I_end, nptperyear)
        # In TIMESAT global model function, `nextent` is not used.
        I_extend <- I

        ## 2. input data
        ti   <- doys[I_extend]
        yi   <- INPUT$y[I_extend]
        wi   <- w[I_extend]
        # yi_good <- yi[w0[I_extend] > wmin]

        ## update ylu in a three year moving window
        ylu <- get_ylu(y, years, w0, width_ylu, I_beg:I_end, Imedian = TRUE, opt$wmin)
        ylu <- merge_ylu(INPUT$ylu, ylu)
        # yi[yi < ylu[1]] <- ylu[1] # update y value

        # beginI <- ifelse(i == 1, 1, 2) # make sure no overlap
        beginI = 1
        tout   <- doys[I] %>% {.[beginI]:last(.)} # make sure return the same length result.

        fit  <- curvefit(yi, ti, tout, nptperyear = nptperyear,
                         w = wi, ylu = ylu,
                         iters = opt$iters, methods = opt$methods, wFUN = opt$wFUN,
                         type = type,
                         ...)
        # if too much missing values
        if (sum(wi > pmax(opt$wmin+0.1, 0.2))/length(wi) < opt$minPercValid){
            fit$model %<>% map(function(x){
                x$zs %<>% map(~.x*NA) # list()
                return(x)
            })
        }
        return(fit)
    }

    fits = list()
    for (i in 1:nlocal) {
        # odd: t2t (1, 3, 5, ...)
        # even: p2k (2, 4, 6)
        # typeI = ifelse(i %% 2 == 1, "t2t", "p2p")
        typeI = ifelse(i %% 2 == 1, 1L, -1L)

        if (opt$verbose) fprintf("  [curvefits] running %d ... \n", i)
        I = s[i]:s[i+2]
        data = list(t = doys[I], y = y0[I], QC_flag = QC_flag[I]) %>% as.data.table()
        fit = curvefit_gs(s[i], s[i+2], type = typeI)
        fit$data = data
        fit$type = typeI
        fits[[i]] = fit
    }
    names(fits)[seq(1, nlocal, 2)] = brks$dt$flag
    return(fits)
    # L1:curve fitting method, L2:yearly flag
    # return(list(tout = t[first(di$beg):last(di$end)],  # dates for OUTPUT curve fitting VI
    #             fits = fits))
}

#' cutoff
#'
#' @param n1,n2 peak and trough (or trough and peak)
#'
#' @return A numeric vector, in small intervals around `(tL + tC)/2` and `(tC+tR)/2`,
#' respectively, smoothly drop from 1 to 0.
#' @references
#' 1. Per J\"onsson, P., Eklundh, L., 2004. TIMESAT - A program for analyzing
#'     time-series of satellite sensor data. Comput. Geosci. 30, 833-845.
#'     https://doi.org/10.1016/j.cageo.2004.05.006.
#' @keywords internal
cutoff <- function(n1, n2, k = n1:n2) {
    ndiff = n2 - n1;
    nmid  = (n1 + n2) / 2;
    # k = n1:n2;
    1 - (atan(10 * (k - nmid) / ndiff) - atan(10 * (n1 - nmid) / ndiff)) /
        (atan(10 * (n2 - nmid) / ndiff) - atan(10 * (n1 - nmid) / ndiff))
}

#' @rdname curvefits_LocalModel
#' @export
merge_LocalModels <- function(fits) {
    # begin from `t2t`; end at `t2t`
    nlocal = length(fits)
    N_gs = (nlocal + 1) / 2

    fits2 = list()
    for (i in 1:N_gs) {
        # i_c = c(di$beg[i], di$end[i])
        # i_l = c(di$peak[i-1], di$peak[i])
        # i_r = c(di$peak[i], di$peak[i+1])
        k = (i - 1) * 2 + 1
        fit_c = fits[[k]]

        fit_l = if(k <= 1)    NULL else fits[[k-1]]
        fit_r = if(k >= N_gs) NULL else fits[[k+1]]
        fit = .merge_LocalModel(fit_c, fit_l, fit_r)
        fits2[[i]] <- fit
    }
    names(fits2) = names(fits)[seq(1, nlocal, 2)]
    fits2
}

# ' Merge Local Model Functions
# ' @references
# ' TIMESAT Mannual v3.3 (2018), Eq. 10
# ' @rdname curvefits_LocalModel
# ' @export 
.merge_LocalModel <- function(fit_c, fit_l = NULL, fit_r = NULL) {
    if (is.null(fit_l) && is.null(fit_r)) return(fit_c)

    t_c = fit_c$tout
    t_l = fit_l$tout
    t_r = fit_r$tout

    info_l = match2(t_l, t_c)
    info_r = match2(t_r, t_c)

    a = if (!is.null(fit_l)) cutoff(first(info_l$I_x), last(info_l$I_x)) else 0
    b = if (!is.null(fit_r)) cutoff(first(info_r$I_x), last(info_r$I_x)) else 1

    # for z_c
    ind_l = info_l$I_y
    ind_r = info_r$I_y
    if (is.null(fit_l) && !is.null(fit_r)) ind_l = -(ind_r)
    if (is.null(fit_r) && !is.null(fit_l)) ind_r = -(ind_l)

    nmeth = fit_c$model %>% length()
    niter = fit_c$model[[1]]$zs %>% length()
    iterator = set_names(1:niter, paste0("iter", 1:niter))

    for(i in 1:nmeth) {
        zs = map(iterator, function(j){
        # for (j in 1:niter) {
            z_c = fit_c$model[[i]]$zs[[j]]
            z_l = if (is.null(fit_l)) 0 else fit_l$model[[i]]$zs[[j]][info_l$I_x]
            z_r = if (is.null(fit_r)) 0 else fit_r$model[[i]]$zs[[j]][info_r$I_x]

            z_cl = z_c[ind_l]
            z_cr = z_c[ind_r]
            # TODO: 目前算法可能会导致zc变短，L84 `beginI = 1`可消除此错误
            ans = c(a * z_l + (1 - a) * z_cl,
              b * z_cr + (1 - b) * z_r)
            # plot(z_c, col = "grey", type = "l", lwd = 1, ylim = c(0, 0.6))
            # lines(ind_l, z_l)
            # lines(ind_r, z_r)
            # lines(ans, col = "red")
            ans
        })
        fit_c$model[[i]]$zs = zs
    }
    fit_c
}
