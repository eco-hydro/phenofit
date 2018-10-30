source('test/stable/load_pkgs.R')

## 1. GLOBAL parameters #######################################################
# To compare with different sites, Roughtness is normalized by `Yobs`.
typename     <- "part"

indice       <- c("R2", "Bias", "RMSE", "Rg_norm_by_obs")
indice_label <- c("bold((a)*' '*R^2)", "bold((b)*' '*Bias)",
                  "bold((c)*' '*RMSE)", "bold((d)*' '*Roughness)")

indir <- "V:/result/gee_whittaker/valid_v3_adjparam"
prefix <- "valid_v3_adjparam"
# prefix <- "valid_v3_noadjparam"

files <- dir("V:/result/gee_whittaker", "*.RDS", full.names = T) %>% last() %>%
    set_names(basename(.) %>% gsub("0.RDS", "", .))
# lst <- get_sbatch(dirname(indir),
#            paste0(prefix, ".*RDS"), Save = F)
lst   <- llply(files, readRDS, .progress = "text")
x     <- lst$valid_v3_noadjparam

## 2. Visualization ###########################################################
text_size = 5

colors <- scales::hue_pal()(3)
colors <- c(colors[2], "grey20", colors[1])

nmax <- 424 # max length of good values for a point
# check_roughtFitting <- function(x, prefix){
    df_info <- map(x, ~melt_list(., "type")) %>%
        {melt_list(., "meth")[iters == "iter1", ]}
    df_info <- df_info[n_sim > 40, ] # at least 10% good values

    # df_info  <- melt_list(lst, "meth")[iters == "iter1", ]
    df_info2 <- df_info[, .(site, iters, type, meth, R2, Bias, RMSE,
                            Rg, Rg_norm_by_obs, Rg_norm_by_pred)] %>%
        melt(.(site, iters, type, meth) %>% names(), variable.name = "index") %>%
        spread(meth, value) %>%
        melt(.(site, iters, type, index, wWH) %>% names(), variable.name = "meth")

    df_info2$index %<>% factor(indice, indice_label)
    df_info2 <- df_info2[!is.na(index) & type == typename,]
    ## 1. set critical values for different index
    DEFAULT.VALUE <- .02
    index <- .(R2, Bias, RMSE, Rg, Rg_norm_by_obs, Rg_norm_by_pred) %>% names() %>%
        factor(indice, indice_label)
    dmax  <- c(.05, .02, .02, .005, .005, .005)

    # index_all   <- levels(df_info2$index) %>% union(index, .) %>% factor(., .)
    # dmax_all    <- rep(DEFAULT.VALUE, length(index_all))
    # I           <- match(index, index_all) %>% rm_empty()
    # dmax_all[I] <- dmax

    # d_diff0 <- data.table(index = index_all, dmax = dmax_all)
    d_diff0 <- data.table(index = index, dmax = dmax) %>% .[!is.na(index), ]
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
    d_diff$value %<>% clamp(c(-2, 2))
    d_diff[type == typename & !is.na(index), ]
    # ps <- list()
    # for (i in seq_along(indice)){
    #     indexi = indice_label[i]
    d_dominant <- ddply_dt(d_diff, .(table(kind)), .(meth, index, type, label))
    d_dominant[, text := sprintf("%2d vs %2d", Bigger, Smaller)]
    d_dominant[, `:=`(
        text_big   = sprintf("bold('%s' * phantom(' vs ' * '%s'))", Bigger, Smaller),
        text_vs    = sprintf("bold(phantom('%s') * ' vs ' * phantom('%s'))", Bigger, Smaller),
        text_small = sprintf("bold(phantom('%s vs ') * '%s')", Bigger, Smaller))]

#     list(d_diff, d_dominant = d_dominant)
# }

# fun_Figure9(pdat, d_dominant)
# fun_Figure9 <- function(pdat, d_dominant, typename = "part", prefix = Sys.Date()){
    d_dom  <- d_dominant[meth == "wWH2"]
    pdat_i <- d_diff[meth == "wWH2"]

    p <- ggplot(pdat_i, aes(wWH, value, color = kind)) +
        geom_point(data = pdat_i[kind != "Similar"], alpha = 0.4) +
        geom_point(data = pdat_i[kind == "Similar"], alpha = 0.1) +
        # facet_wrap(~label, scales = "free") +
        facet_wrap(.~index, scales = "free", labeller = label_parsed) + #, ncol = 4, labeller = label_value
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

    file <- sprintf("Fig10_compare_with_wWH2_%s.pdf", paste0(prefix, "_" ,typename))
    file_tiff <- gsub(".pdf", ".tif", file)
    width = 14; height = 7
    write_fig(p, file_tiff, width, height, T, res = 300)
# }s

# fun_Figure11 <- function(d_diff, d_dom){
    text_size = 5
    ps <- list()
    for (i in seq_along(indice_label)){
        index_i = indice_label[i]

        d_dom    <- d_dominant[index == index_i]
        pdat_i   <- d_diff[index == index_i]

        p <- ggplot(pdat_i, aes(wWH, value, color = kind)) +
            geom_point(data = pdat_i[kind != "Similar"], alpha = 0.2) +
            geom_point(data = pdat_i[kind == "Similar"], alpha = 0.2) +
            # facet_wrap(~label, scales = "free") +
            facet_grid(meth~index, scales = "free", labeller = label_parsed) + #, ncol = 4, labeller = label_value
            geom_abline(slope = 1, color ="red", size = 0.5) +
            scale_color_manual(values = colors) +

            geom_text(data = d_dom, aes(label = text_big, color = NULL), parse = T,
                      color = colors[1], fontface = "bold", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, show.legend = F, size = text_size) +
            geom_text(data = d_dom, aes(label = text_vs, color = NULL), parse = T,
                      color = "black", fontface = "bold", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, show.legend = F, size = text_size) +
            geom_text(data = d_dom, aes(label = text_small, color = NULL), parse = T,
                      color = colors[3], fontface = "bold", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, show.legend = F, size = text_size) +
            ylab("Other methods") +
            theme_gray(base_size = 14) +
            theme(legend.title = element_blank(),
                  legend.position = "top",
                  legend.margin = margin(),
                  legend.spacing.x = unit(0.2, 'cm'),
                  panel.grid.minor = element_blank(),
                  axis.title = element_blank(),
                  strip.text = element_text( margin = margin(3, 3, 3, 3, "pt"), size = 14, face = "bold"),
                  plot.margin = margin(2, 3, 0, 0, "pt"))
        # guides(color = guide_legend(keywidth = 2))
        if (i < 4){
            p <- p + theme(strip.text.y = element_blank(),
                           legend.position = "none")
        }
        # if (i > 1){
        #     p <- p + theme(axis.title.y = element_blank())
        # }
        ps[[i]] <- p
    }
    ps

    lgd <- g_legend(p)
    p   <- p + theme(legend.position = "none")
    ps[[i]] <- p

    # bottom <- arrangeGrob(,
    #                       lgd)
    xtitle <- textGrob("wWHd", gp=gpar(fontsize=14, fontface = "bold"))
    g <- arrangeGrob(grobs = ps, nrow = 1, widths = c(1, 1, 1, 1.1),
                     bottom = xtitle,
                     left = textGrob("Other methods", gp=gpar(fontsize=14, fontface = "bold"), rot = 90))
    g <- arrangeGrob(g, bottom = lgd)
    file <- sprintf("Fig11_compare_with_other_methods_%s.pdf", typename)
    file_tiff <- gsub(".pdf", ".tif", file)
    width = 11; height = 8
    write_fig(g, file_tiff, width, height, T, res = 300)
    # write_fig(g, file, width, height, T)
# }
