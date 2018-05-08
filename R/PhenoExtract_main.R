#'
#' curve fit vegetation index (VI) time-series
#'
#' curve fit VI using 'spline', 'beck', 'elmore', 'klos', 'AG', 'Gu', 'zhang' methods
#'
#' @param x Vegetation time-series index, numeric vector
#' @param t The corresponding doy of x
#' @param tout The output interpolated time.
#'
#' @export
curvefit <- function(x, t = index(x), tout = t, meth = 'BFGS',
    methods = c('spline', 'beck', 'elmore', 'klos', 'AG', 'zhang'), ...)
{
    if (all(is.na(x))) return(NULL)
    if (length(methods) == 1 && methods == 'all')
        methods <- c('spline', 'beck', 'elmore', 'klos', 'AG', 'Gu', 'zhang')

    params <- list(x, t, tout, optimFUN = I_optim, method = meth, ...)

    if ('spline' %in% methods) fit.spline <- splinefit(x, t, tout)

    if ('beck'   %in% methods) fit.beck   <- do.call(FitDL.Beck,  c(params, pfun = p_nlminb))  #nlminb
    if ('elmore' %in% methods) fit.elmore <- do.call(FitDL.Elmore,c(params, pfun = p_nlminb))  #nlminb

    # best: BFGS, but its speed lower than other function, i.e. nlm
    if ('Gu'     %in% methods) fit.Gu     <- do.call(FitDL.Gu,    c(params, pfun = p_nlminb))  #nlm, ucminf
    if ('zhang'  %in% methods) fit.zhang  <- do.call(FitDL.Zhang, c(params, pfun = p_nlminb))  #nlm
    if ('AG'     %in% methods) fit.AG     <- do.call(FitAG,       c(params, pfun = p_nlminb))     #nlm
    if ('klos'   %in% methods) fit.klos   <- do.call(FitDL.Klos,  c(params, pfun = p_optim))   #BFGS, Nelder-Mead, L-BFGS-B

    names <- ls(pattern = "fit\\.") %>% set_names(., .)
    fits  <- lapply(names, get, envir = environment()) %>%
        set_names(toupper(gsub("fit\\.", "", names))) #remove `fit.` and update names
    return(fits)
}

#' plot_phenofit
#'
#' @param fit data from phenofit_site
#' @importFrom dplyr left_join
#' @export
plot_phenofit <- function(fit, d, title = NULL, plotly = F){
    #global variables
    origin.date <- fit$INPUT$t[1] #
    I   <- match(fit$tout, fit$INPUT$t)
    org <- fit$INPUT[I, ]

    getFitVI <- function(fit_years){
        map(fit_years, ~window(.x$pred, .x$data$t)) %>%
            do.call(c, .) %>% set_names(NULL)
    }
    # fit_years: fits for one method, multiple years
    getFitVI_iters <- function(fit_years){
        if (is.null(fit_years[[1]]$fun)){
            out <- getFitVI(fit_years) %>% as.numeric() %>% {tibble(iter1 = .)}
        }else{
            out <- map(fit_years, function(fit){
                I       <- match(fit$data$t, fit$tout)
                d_iters <- map_df(fit$fits, ~.x[I])
                d_iters$t <- fit$data$t
                return(d_iters)
            }) %>% bind_rows()
        }
        # out
        # # A tibble: 2,515 x 3
        # iter1 iter2 t
        # <dbl> <dbl> <date>
        # 1 0.594 0.780 2005-02-06
        # 2 0.595 0.780 2005-02-07
        # 3 0.595 0.780 2005-02-08
        out$t %<>% add(origin.date - 1)
        left_join(org, out, by = "t")
    }

    fits_years <- map(fit$fits, getFitVI_iters)
    pdat1 <- melt(fits_years, id.vars = c("t", "y")) %>%
        set_names(c("t", "y","iters", "val", "meth"))
    # fits_years <- map(fit$fits, getFitVI)

    # t_fit <- index(fits_years[[1]]) + t[1] - 1
    # new   <- as_tibble(c(list(t = t_fit), map(fits_years, unclass)))
    # new   <- as_tibble(c(list(t = t_fit), fits_years))
    # 1. curve fitting data
    # pdat1 <- dplyr::full_join(fit$INPUT, new, by = "t") %>%
    #     gather(meth, val, -t, -y) %>% group_by(meth)
    # print(head(pdat1))

    # 2. growing season breaks
    # 3. phenology data
    pdat2 <- fit$pheno$date %>% melt_list("meth") %>% as_tibble() %>%
        gather(index, date, -flag, -origin, -meth) %>%
        mutate(pmeth = str_extract(index, "\\w{1,}"))
    # pdat2 <- pdat2[grep("TRS2", pdat2$index), ]
    # 3. growing season breaks
    # try to add whittaker smoother here
    seasons <- fit$seasons
    seasons$pos$type %<>% factor(labels = c("min", "max"))

    p <- ggplot(pdat1, aes(t, val, color = iters)) +
        geom_line (data = seasons$whit, aes(t, iter3), color = "black", size = 0.6) +
        geom_vline(data = seasons$dt, aes(xintercept = as.numeric(beg)), size = 0.3, linetype=2, color = "blue") +
        geom_vline(data = seasons$dt, aes(xintercept = as.numeric(end)), size = 0.3, linetype=2, color = "red") +
        geom_point(data = subset(seasons$pos, type == "max"), aes(t, val), color= "red") +
        geom_point(data = subset(seasons$pos, type == "min"), aes(t, val), color= "blue") +
        # geom_line(size = 0.4) +
        facet_grid(meth~.) +
        scale_x_date(breaks = fit$seasons$dt$beg, date_labels = "%Y/%m")

    if ('SummaryQA' %in% colnames(d)){
        # p <- p + geom_point(data = x, aes(date, y, shape = SummaryQA, color = SummaryQA), size = 1, alpha = 0.7) +
        #     scale_shape_manual(values = c(21,22, 24:25)) +
        #     scale_fill_manual(values = c("grey40", "#7CAE00", "#F8766D", "#C77CFF")) +
        #     guides(shape = guide_legend(override.aes = list(size = 2)))
        p  <- p + geom_point(data = d, aes(date, y, shape = SummaryQA, color = SummaryQA), size = 2, alpha = 0.7) +
            geom_line(aes(color = iters), size = 1, alpha = 0.7) +
            scale_color_manual(values = c(" good" = "grey60", " margin" = "#00BFC4",
                                          " snow&ice" = "#F8766D", " cloud" = "#C77CFF",
                                          "iter1" = "blue", "iter2" = "red"), drop = F) +
            scale_shape_manual(values = c(19, 15, 4, 17), drop = F) +
            guides(shape = FALSE,
                   color = guide_legend(
                       "legend",
                       override.aes = list(
                           shape = c(c(19, 15, 4, 17)[c(4, 1, 2, 3)], NA, NA),
                           linetype = c(0, 0, 0, 0, 1, 1),
                           # color = c(" good" = "grey60", " margin" = "#00BFC4",
                           #           " snow&ice" = "#F8766D", " cloud" = "#C77CFF",
                           #           "iter1" = "blue", "iter2" = "red"),
                           name = letters[1:6],
                           size = 1.2
                       )
                   ))
    }else{
        p <- p + geom_point(aes(t, y), size = 2, alpha = 0.5, color = "grey60") +
            geom_line(aes(color = iters), size = 1)
    }
    # geom_vline(data = pdat2, aes(xintercept=date, linetype = pmeth, color = pmeth), size = 0.4, alpha = 1) +
    # scale_linetype_manual(values=c(2, 3, 1, 1, 1, 4))
    # p + facet_grid(meth~pmeth)
    #return
    if (plotly){
        plotly::ggplotly(p)
    }else{
        p + ggtitle(title)
    }
}


#' tidyFits_pheno
#'
#' Tidy for every method with multiple years phenology data
#' @export
tidyFits_pheno <- function(pheno, origin){
    tidToDate <- function(datenum)  origin - 1 + unlist(datenum)
    # phenonames <- c('TRS2.sos', 'TRS2.eos', 'TRS5.sos', 'TRS5.eos', 'TRS6.sos', 'TRS6.eos',
    #                 'DER.sos', 'DER.pop', 'DER.eos',
    #                 'GU.UD', 'GU.SD', 'GU.DD', 'GU.RD',
    #                 'ZHANG.Greenup', 'ZHANG.Maturity', 'ZHANG.Senescence', 'ZHANG.Dormancy')
    p_date <- ldply(pheno, tidToDate, .id = "flag") %>%
      mutate(origin = ymd(paste0(substr(flag, 1, 4), "-01-01"))) %>% as_tibble()
    phenonames <- contain(p_date,  "\\.")

    df <- gather(p_date, meth, date, -flag, -origin) %>%
        mutate(doy = as.numeric(date - origin + 1))

    # date <- spread(pheno[, c("flag", "meth", "date")], meth, date)
    p_doy  <- spread(df[, c("flag", "meth", "doy", "origin")], meth, doy)

    vars  <- c("flag", "origin", phenonames)
    list(date = p_date[vars], doy = p_doy[vars]) #quickly return
}

phenonames <- c('TRS2.SOS', 'TRS2.EOS', 'TRS5.SOS', 'TRS5.EOS', 'TRS6.SOS', 'TRS6.EOS',
    'DER.SOS', 'DER.POP', 'DER.EOS',
    'UD', 'SD', 'DD','RD',
    'GreenUp', 'Maturity', 'Senescence', 'Dormancy')

#' curvefit_optimx
#'
#' With the help of `optimx` package, try to find which optimization function
#' is best.
#' @param optimFUN `I_optim` or `I_optimx`
#' @param meth c('BFGS','CG','Nelder-Mead','L-BFGS-B','nlm','nlminb',
#' spg','ucminf','Rcgmin','Rvmmin','newuoa','bobyqa','nmkb','hjkb')
#' @param pfun c(p_nlminb, p_ncminf, p_nlm, p_optim)
#'
#' @export
curvefit_optimx <- function(x, t = index(x), tout = t,
    optimFUN = I_optimx,
    methods = c('spline', 'beck', 'elmore', 'AG'),
    meth, pfun, ...)
{
    if (all(is.na(x))) return(NULL)
    ##2. curve fitting
    if (length(methods) == 1 && methods == 'all')
        methods <- c('spline', 'beck', 'elmore', 'klos', 'AG', 'Gu', 'zhang')

    # failed: 'BFGS', 'Nelder-Mead', 'L-BFGS-B'
    # meth = 'L-BFGS-B'
    # ok: 'L-BFGS-B'; failed: 'BFGS', 'Nelder-Mead'
    params <- list(x, t, tout, optimFUN = optimFUN, pfun = pfun, method = meth, ...)
    # fit.beck   <- FitDL.Beck   #even Nelder-Mead was faster and convergent, but nlminb was better
    # ok: BFGS; failed: 'L-BFGS-B'
    if ('spline' %in% methods) fit.spline <- splinefit(x, t, tout)
    if ('beck'   %in% methods) fit.beck   <- do.call(FitDL.Beck, params)        #nlminb
    if ('elmore' %in% methods) fit.elmore <- do.call(FitDL.Elmore, params)      #nlminb

    # best: BFGS, but its speed lower than other function, i.e. nlm
    if ('Gu'     %in% methods) fit.Gu     <- do.call(FitDL.Gu, params)          #nlm, ucminf
    if ('zhang'  %in% methods) fit.zhang  <- do.call(FitDL.Zhang, params)       #nlm
    if ('AG'     %in% methods) fit.AG     <- do.call(FitAG, params)             #nlm
    if ('klos'   %in% methods) fit.klos   <- do.call(FitDL.Klos, params)        #BFGS, Nelder-Mead, L-BFGS-B

    # test for optimx methods
    # fit   <- FitDL.Zhang  (x, t, tout, optimFUN = optimx_fun, debug = T, method = 'BFGS')
    names <- ls(pattern = "fit\\.") %>% set_names(., .)
    fits  <- lapply(names, get, envir = environment()) %>%
        set_names(toupper(gsub("fit\\.", "", names))) #remove `fit.` and update names
    return(fits)
}
