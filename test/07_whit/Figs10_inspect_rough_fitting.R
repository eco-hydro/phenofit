check_roughtFitting <- function(x, prefix){
    df_info <- map(x, ~melt_list(., "type")) %>%
        {melt_list(., "meth")[iters == "iter1", ]}
    df_info <- df_info[n_sim > 40, ] # at least 10% good values

    # df_info  <- melt_list(lst, "meth")[iters == "iter1", ]
    df_info2 <- df_info[, .(site, iters, type, meth, R2, Bias, RMSE, Rg, Rg_norm_by_obs, Rg_norm_by_pred)] %>%
        melt(.(site, iters, type, meth) %>% names(), variable.name = "index") %>%
        spread(meth, value) %>%
        melt(.(site, iters, type, index, wWH) %>% names(), variable.name = "meth")

    ## 1. set critical values for different index
    DEFAULT.VALUE <- .02
    index <- .(R2, Bias, RMSE, Rg, Rg_norm_by_obs, Rg_norm_by_pred) %>% names()
    dmax  <- c(.05, .02, .02, .005, .005, .005)

    index_all   <- levels(df_info2$index) %>% union(index, .) %>% factor(., .)
    dmax_all    <- rep(DEFAULT.VALUE, length(index_all))
    I           <- match(index, index_all) %>% rm_empty()
    dmax_all[I] <- dmax

    d_diff0 <- data.table(index = index_all, dmax = dmax_all)

    ## 2. visualization
    d_diff <- merge(df_info2, d_diff0)
    d_diff[, diff := value - wWH]

    # a <- d_diff[type == "part" & meth == "wWH2" & index == "R2" & diff > 0.05]
    # write_rds(a, "wWHd_inspect_sites.RDS")

    d_diff$kind <- 0
    d_diff[wWH - value >  dmax, kind := 1]
    d_diff[wWH - value < -dmax, kind := -1]
    d_diff$kind %<>% factor(levels = c(1, 0, -1),
                        labels = c("Bigger", "Similar", "Smaller"))
    d_diff[, label:= sprintf("%s~%s", meth, index)]
    d_diff$label %<>% factor()

    # ps <- list()
    # for (i in seq_along(indice)){
    #     indexi = indice_label[i]
    d_dominant <- ddply_dt(d_diff, .(table(kind)), .(meth, index, type, label))
    d_dominant[, text := sprintf("%2d vs %2d", Bigger, Smaller)]
    d_dominant[, `:=`(
        text_big   = sprintf("bold('%s' * phantom(' vs ' * '%s'))", Bigger, Smaller),
        text_vs    = sprintf("bold(phantom('%s') * ' vs ' * phantom('%s'))", Bigger, Smaller),
        text_small = sprintf("bold(phantom('%s vs ') * '%s')", Bigger, Smaller))]

    pdat <- d_diff
    pdat$value %<>% clamp(c(-2, 2))

    p1 <- vis_fun(pdat, d_dominant, "all")
    p2 <- vis_fun(pdat, d_dominant, "part")

    # itersI <- 0
    # file <- sprintf("Fig10_compare_with_wWH2_%s.pdf", paste0(prefix, "_all"))
    # file_tiff <- gsub(".pdf", ".tif", file)
    # width = 14; height = 7
    # write_fig(p1, file_tiff, width, height, T, res = 300)

    file <- sprintf("Fig10_compare_with_wWH2_%s.pdf", paste0(prefix, "_good"))
    file_tiff <- gsub(".pdf", ".tif", file)
    width = 14; height = 7
    write_fig(p2, file_tiff, width, height, T, res = 300)
}

vis_fun <- function(pdat, d_dominant, typename = "all"){
    d_dom <- d_dominant[type == typename]
    pdati  <- pdat[type == typename]

    ggplot(pdati, aes(wWH, value, color = kind)) +
        geom_point(data = pdati[kind != "Similar"], alpha = 0.4) +
        geom_point(data = pdati[kind == "Similar"], alpha = 0.1) +
        # facet_wrap(~label, scales = "free") +
        # facet_wrap(.~index, scales = "free", labeller = label_parsed) + #, ncol = 4, labeller = label_value
        facet_grid(meth~index, scales = "free", labeller = label_parsed) + #, ncol = 4, labeller = label_value
        geom_abline(slope = 1, color ="red", size = 0.5) +
        scale_color_manual(values = colors) +

        geom_text(data = d_dom, aes(label = text_big, color = NULL), parse = T,
                  color = colors[1], fontface = "bold", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, show.legend = F, size = text_size) +
        geom_text(data = d_dom, aes(label = text_vs, color = NULL), parse = T,
                  color = "black", fontface = "bold", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, show.legend = F, size = text_size) +
        geom_text(data = d_dom, aes(label = text_small, color = NULL), parse = T,
                  color = colors[3], fontface = "bold", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, show.legend = F, size = text_size) +
        # ylab("wWH2") +
        labs(x = "wWHd", y = "wWH2") +
        theme_gray(base_size = 12) +
        theme(legend.title = element_blank(),
              legend.position = "bottom",
              legend.text = element_text(size=12),
              legend.margin = margin(),
              legend.spacing.x = unit(0.2, 'cm'),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size = 14, face = "bold"),
              strip.text = element_text( margin = margin(0.9, 1, 1, 1, "pt")*1.5, size = 12, face = "bold"),
              plot.margin = margin(2, 3, 0, 0, "pt"))
}
################################################################################


res  <- lst %>% {mapply(check_roughtFitting, ., names(.))}
