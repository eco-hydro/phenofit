# get curve fitting of spline, deprecated
getFitVI <- function(fit_years){
    map(fit_years, ~window(.x$pred, .x$data$t)) %>% do.call(c, .) %>% set_names(NULL)
}

getFitValueYear <- function(fit_year, first = FALSE){
    t         <- fit_year$data$t
    I         <- match(t, fit_year$tout)
    d_iters   <- map_df(fit_year$fits, ~.x[I])
    d_iters$t <- t
    if (!first) d_iters <- d_iters[-1, ]
    return(d_iters)
}

#' getFitValueYears
#' @examples
#' origin.date <- fit$INPUT$t[1]
#' out$t %<>% add(origin.date - 1)
getFitValueYears <- function(fit_years){
    # Avoiding overlap, only first year include first point, 
    nyear <- length(fit_years)
    first  <- map_df(fit_years[1], getFitValueYear, first = T)
    others <- map_df(fit_years[2:nyear], getFitValueYear)
    rbind(first, others)
}

#' @export
getFittings <- function(fit){
    #global variables
    origin.date <- fit$INPUT$t[1] #
    I   <- match(fit$tout, fit$INPUT$t)
    org <- as_tibble(fit$INPUT[c('t', 'y', 'w')])[I, ]

    # fit_years: fits for one method, multiple years
    getFitVI_iters <- function(fit_years){
        if (!is.null(fit_years[[1]]$fun)){
        #     out <- getFitVI(fit_years) %>% as.numeric() %>% {tibble(iter1 = .)}
        # }else{
            out <- getFitValueYears(fit_years)
        }
        # out
        # iter1 iter2 t
        # <dbl> <dbl> <date>
        # 1 0.594 0.780 2005-02-06
        out$t %<>% add(origin.date - 1)
        left_join(org, out, by = "t")
    }

    fits_years <- map(fit$fits, getFitVI_iters)
    res <- melt(fits_years, id.vars = colnames(org)) %>%
        set_names(c(colnames(org), "iters", "val", "meth"))
    return(res)
}

#' plot_phenofit
#'
#' @param fit data from phenofit_site
#' @importFrom dplyr left_join
#' @export
plot_phenofit <- function(fit, d, title = NULL, show.legend = T, plotly = F){
    pdat1 <- get_fittings(fit)
    # t_fit <- index(fits_years[[1]]) + t[1] - 1
    # 1. curve fitting data
    # pdat1 <- dplyr::full_join(fit$INPUT, new, by = "t") %>%
    #     gather(meth, val, -t, -y) %>% group_by(meth)
    # print(head(pdat1))
    
    # 3. phenology data
    # pdat2 <- fit$pheno$date %>% melt_list("meth") %>% as_tibble() %>%
    #     gather(index, date, -flag, -origin, -meth) %>%
    #     mutate(pmeth = str_extract(index, "\\w{1,}"))
    # pdat2 <- pdat2[grep("TRS2", pdat2$index), ]
    # 3. growing season breaks
    # try to add whittaker smoother here
    seasons <- fit$seasons
    # seasons$pos$type %<>% factor(labels = c("min", "max"))

    p <- ggplot(pdat1, aes(t, val, color = iters)) +
        geom_line (data = seasons$whit, aes(t, iter2), color = "black", size = 0.8) +
        geom_vline(data = seasons$dt, aes(xintercept = as.numeric(beg)), size = 0.4, linetype=2, color = "blue") +
        geom_vline(data = seasons$dt, aes(xintercept = as.numeric(end)), size = 0.4, linetype=2, color = "red") +
        geom_point(data = seasons$dt, aes(peak, y_peak), color= "red") +
        geom_point(data = seasons$dt, aes(beg , y_beg ), color= "blue") +
        # geom_line(size = 0.4) +
        facet_grid(meth~.) +
        scale_x_date(breaks = fit$seasons$dt$beg, date_labels = "%Y/%m") + ggtitle(title)

    if ('SummaryQA' %in% colnames(d)){
        # p <- p + geom_point(data = x, aes(date, y, shape = SummaryQA, color = SummaryQA), size = 1, alpha = 0.7) +
        #     scale_shape_manual(values = c(21,22, 24:25)) +
        #     scale_fill_manual(values = c("grey40", "#7CAE00", "#F8766D", "#C77CFF")) +
        #     guides(shape = guide_legend(override.aes = list(size = 2)))
        p  <- p + geom_point(data = d, aes(t, y, shape = SummaryQA, color = SummaryQA), size = 2, alpha = 0.7) +
            geom_line(aes(color = iters), size = 0.8, alpha = 0.7) +
            scale_color_manual(values = c(" good" = "grey60", " margin" = "#00BFC4",
                                          " snow&ice" = "#F8766D", " cloud" = "#C77CFF",
                                          "iter1" = "blue", "iter2" = "red"), drop = F) +
            scale_shape_manual(values = c(19, 15, 4, 17), drop = F) +
            ylab('Vegetation Index')
            # guides(shape = FALSE,
            #        color = guide_legend(
            #            "legend",
            #            override.aes = list(
            #                shape = c(c(19, 15, 4, 17)[c(4, 1, 2, 3)], NA, NA),
            #                linetype = c(0, 0, 0, 0, 1, 1),
            #                # color = c(" good" = "grey60", " margin" = "#00BFC4",
            #                #           " snow&ice" = "#F8766D", " cloud" = "#C77CFF",
            #                #           "iter1" = "blue", "iter2" = "red"),
            #                name = letters[1:6],
            #                size = 1.2
            #            )
            #        ))
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
        if (show.legend){
            p <- p + theme(legend.position="none")
            p <- arrangeGrob(p, lgd, nrow = 2, heights = c(15, 1), padding = unit(1, "line")) #return, 
        }
        return(p)
    }
}

#' make_legend
make_legend <- function(){
    labels <- c(" good", " margin", " snow/ice", " cloud", "iter1", "iter2", "whit")
    colors <- c("grey60", "#00BFC4", "#F8766D", "#C77CFF", "blue", "red", "black")
    pch <- c(19, 15, 4, 17, NA, NA, NA)
    lty <- c(0, 0, 1, 0, rep(1, 3))
    lwd <- c(1, 1, 1, 1, 3, 3, 3)

    I <- 1:7
    lgd <- grid::legendGrob(labels[I], pch = pch[I], nrow = 1,
                       # do.lines = T,
                       gp=grid::gpar(lty = lty[I], lwd = lwd[I],
                               col = colors[I], fill = colors[I]))
    return(lgd)
}
lgd <- make_legend()

#' tidyFits_pheno
#'
#' Tidy for every method with multiple years phenology data
#' @export
tidyFitPheno <- function(pheno, origin){
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
