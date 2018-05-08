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
        #     1 0.594 0.780 2005-02-06
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
