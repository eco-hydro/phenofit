#' curvefits_cutoff
#' 
#' @keywords internal
#' @example inst/examples/ex-curvefits.R
#' @export
curvefits_cutoff <- function(
    INPUT, brks,
    # methods, wFUN,
    options = list(),
    globalModel = FALSE,
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

    # growing season dividing
    di <- data.table( beg  = getDateId_before(brks$dt$beg, t),
                      peak = getDateId_before(brks$dt$peak, t),
                      end  = getDateId_after(brks$dt$end, t)) #%>% na.omit()

    browser()
    width_ylu = nptperyear*2

    y    <- INPUT$y
    fits <- vector(nrow(di), mode = "list")
    browser()
    
    for (i in 1:nrow(di)){
        if (opt$verbose) fprintf("  [curvefits] running %d ... \n", i)

        I     <- di$beg[i]:di$end[i]
        I_beg <- di$beg[i]
        I_end <- di$end[i]

        I_extend <- get_extentI(w0, I_beg, I_end, nptperyear)

        ## 2. input data
        ti   <- doys[I_extend]
        yi   <- INPUT$y[I_extend]
        wi   <- w[I_extend]
        # yi_good <- yi[w0[I_extend] > wmin]

        ## update ylu in a three year moving window
        ylu <- get_ylu(y, years, w0, width_ylu, I_beg:I_end, Imedian = TRUE, opt$wmin)
        ylu <- merge_ylu(INPUT$ylu, ylu)
        # yi[yi < ylu[1]] <- ylu[1] # update y value

        beginI <- ifelse(i == 1, 1, 2) # make sure no overlap
        tout   <- doys[I] %>% {.[beginI]:last(.)} # make sure return the same length result.

        fFITs  <- curvefit(yi, ti, tout, nptperyear = nptperyear,
                         w = wi, ylu = ylu,
                         iters = opt$iters, methods = opt$methods, wFUN = opt$wFUN,
                         ...)
        # add original input data here, global calculation can comment this line
        # `y` is original time-series without checked, This for plot
        data <- list(t = doys[I], y = y0[I], QC_flag = QC_flag[I]) %>% as.data.table()
        fFITs$data <- data

        # if too much missing values
        if (sum(wi > pmax(opt$wmin+0.1, 0.2))/length(wi) < opt$minPercValid){
            fFITs$model %<>% map(function(x){
                x$zs %<>% map(~.x*NA) # list()
                return(x)
            })
        }
        fits[[i]] <- fFITs
    }
    # L1:curve fitting method, L2:yearly flag
    fits %<>% set_names(brks$dt$flag)
    fits
    # return(list(tout = t[first(di$beg):last(di$end)],  # dates for OUTPUT curve fitting VI
    #             fits = fits))
}
