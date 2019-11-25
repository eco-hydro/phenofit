source("test/load_pkgs.R")

# load("data_test/whit_lambda/MOD13A1_st_1e3_20180731.rda")
# st$site %<>% as.character()

## 1. GLOBAL parameters #######################################################
# To compare with different sites, Roughtness is normalized by `Yobs`.
typename  <- "all"
str_iter  <- "iter1"
methods   <- c("AG", "ZHANG", "smooth_wHANTS", "smooth_wSG", "wWH2")

indice       <- c("R2", "Bias", "RMSE", "Rg_norm_by_obs")
indice_label <- c("bold((a)*' '*R^2)", "bold((b)*' '*Bias)",
                  "bold((c)*' '*RMSE)", "bold((d)*' '*Roughness)")

indir <- "V:/result/gee_whittaker/valid_v3_adjparam"
# prefix <- "valid_v3_noadjparam"

# phenofit curve fitting methods-098762076212
df_pt <- readRDS_tidy("file:///V:/result/gee_whittaker/phenofit_st1e3 (20180910).RDS") %>%
    do.call(c, .) %>% melt_list("site") %>%
    reorder_name(c("site", "meth", "iters", "type"))

# valid_v3_adjparam: -2

files <- dir("V:/result/gee_whittaker", "*.RDS", full.names = T) %>% nth(-2) %>%
    set_names(basename(.) %>% gsub("0.RDS", "", .))
prefix <- names(files[1]) # "valid_v3_adjparam"
# lst <- get_sbatch(dirname(indir),
#            paste0(prefix, ".*RDS"), Save = F)
lst   <- llply(files, readRDS, .progress = "text")
x     <- lst[[1]]

## 2. Visualization ###########################################################
text_size = 5

colors <- scales::hue_pal()(3)
colors <- c(colors[2], "grey20", colors[1])

nmax <- 424 # max length of good values for a point
# check_roughtFitting <- function(x, prefix){
df_info <- melt_list(x, "meth") %>%
    reorder_name(c("site", "meth", "iters", "type"))
df_info <- rbind(df_info, df_pt)
df_info <- df_info[n_sim > 5, ] # at least 10% good values

# MAIN --------------------------------------------------------------------

# df_info  <- melt_list(lst, "meth")[iters == "iter1", ]
df_info2 <- df_info[iters == "iter1", .(site, iters, type, meth, R2, Bias, RMSE,
                        Rg, Rg_norm_by_obs, Rg_norm_by_pred)] %>%
    melt(.(site, iters, type, meth) %>% names(), variable.name = "index") %>%
    spread(meth, value) %>%
    melt(.(site, iters, type, index, wWH) %>% names(), variable.name = "meth")

df_info2$index %<>% factor(indice, indice_label)
df_info2 <- df_info2[!is.na(index) & type == typename &
                         meth %in% methods & iters == str_iter,]
df_info2$meth %<>% factor(methods)


## 1. set critical values for different index
DEFAULT.VALUE <- .02
index <- .(R2, Bias, RMSE, Rg, Rg_norm_by_obs, Rg_norm_by_pred) %>% names() %>%
    factor(indice, indice_label)
dmax  <- c(.05, .005, .01, rep(.002, 3)) # good: 0.005, all: 0.002

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


# Figure10 ----------------------------------------------------------------

d_dom  <- d_dominant[meth == "wWH2"]
pdat_i <- d_diff[meth == "wWH2"]

kinds <- c("Bigger", "Similar", "Smaller")
kinds_fix <- c("wWHd is larger", "Similar", "wWHd is smaller")
pdat_i$kind %<>% mapvalues(kinds, kinds_fix)

p <- ggplot(pdat_i, aes(wWH, value, color = kind)) +
    scale_color_manual(values = set_names(colors, kinds_fix), breaks = kinds_fix) +
    geom_point(data = pdat_i[kind != "Similar"], alpha = 0.4) +
    geom_point(data = pdat_i[kind == "Similar"], alpha = 0.1) +
    # facet_wrap(~label, scales = "free") +
    facet_wrap(.~index, scales = "free", labeller = label_parsed) + #, ncol = 4, labeller = label_value
    geom_abline(slope = 1, color ="red", size = 0.5) +
    geom_text(data = d_dom, aes(label = text_big, color = NULL), parse = T,
              color = colors[1], fontface = "bold", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, show.legend = F, size = text_size) +
    geom_text(data = d_dom, aes(label = text_vs, color = NULL), parse = T,
              color = "black", fontface = "bold", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, show.legend = F, size = text_size) +
    geom_text(data = d_dom, aes(label = text_small, color = NULL), parse = T,
              color = colors[3], fontface = "bold", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, show.legend = F, size = text_size) +
    # ylab("wWH2") +
    labs(x = expression("Weighted Whittaker with dynamic " * lambda * " (wWHd)"),
         y = expression("Weighted Whittaker with " * lambda * " = 2 (wWH2)")) +
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

file <- sprintf("Figure09_compare_with_wWH2_%s_%s_%s.pdf", prefix, typename, str_iter)
file_tiff <- gsub(".pdf", ".tif", file)
width = 8; height = 6
write_fig(p, file_tiff, width, height, T, res = 300)

# pdat_i[, .(wWH = range(wWH, na.rm = T),
#            value = range(value, na.rm = T))]

# Figure11 ----------------------------------------------------------------

# fun_Figure11 <- function(d_diff, d_dom){
text_size = 5
ps <- list()
for (i in seq_along(indice_label)){
    index_i = indice_label[i]

    pdat_i   <- d_diff[index == index_i & meth != "wWH2"]
    pdat_i$kind %<>% mapvalues(kinds, kinds_fix)

    d_dom      <- d_dominant[index == index_i & meth != "wWH2"]
    d_dom_lab2 <- d_dom[meth %in% levels(meth)[1:2], ]

    lims     <- with(pdat_i, range(c(wWH, value), na.rm = T))
    d_blank  <- data.table(wWH=lims, value=rep(lims, each = 2))

    p <- ggplot(pdat_i, aes(wWH, value, color = kind)) +
        geom_blank(data = d_blank, aes(color = NULL)) +
        geom_point(data = pdat_i[kind != "Similar"], alpha = 0.2) +
        geom_point(data = pdat_i[kind == "Similar"], alpha = 0.2) +
        # facet_wrap(~label, scales = "free") +
        facet_grid(meth~index, scales = "fixed", labeller = label_parsed) + #, ncol = 4, labeller = label_value
        geom_abline(slope = 1, color ="red", size = 0.5) +
        scale_color_manual(values = set_names(colors, kinds_fix), breaks = kinds_fix) +
        geom_text(data = d_dom, aes(label = text_big, color = NULL), parse = T,
                  color = colors[1], fontface = "bold", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, show.legend = F, size = text_size) +
        geom_text(data = d_dom, aes(label = text_vs, color = NULL), parse = T,
                  color = "black", fontface = "bold", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, show.legend = F, size = text_size) +
        geom_text(data = d_dom, aes(label = text_small, color = NULL), parse = T,
                  color = colors[3], fontface = "bold", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, show.legend = F, size = text_size) +
        # ylab("Other methods") +
        # labs(x = expression("Weighted Whittaker with dynamic " * lambda * " (wWHd)"),
        #      y = expression("Other methods"))
        theme_gray(base_size = 14) +
        theme(legend.title = element_blank(),
              legend.position = "top",
              legend.margin = margin(),
              legend.spacing.x = unit(0.2, 'cm'),
              panel.grid.minor = element_blank(),
              axis.title = element_blank(),
              strip.text = element_text( margin = margin(3, 3, 3, 3, "pt"), size = 14, face = "bold"),
              plot.margin = margin(2, 3, 0, 0, "pt"))
    if ((i %in% 1:2) && typename == "all"){
        p <- p +
            geom_text(data = d_dom_lab2, aes(label = text_small, color = NULL),
                      parse = T, color = "white", fontface = "bold", x = -Inf, y = Inf,
                      hjust = -0.1, vjust = 1.5, show.legend = F, size = text_size)
            # geom_text(data = d_dom_lab2, aes(label = text_small, color = NULL), parse = T,
            #            color = "white", fontface = "bold", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, show.legend = F, size = text_size)
    }
    if (i == 2){
        lims <- c(-0.05, 0.05)
        p <- p +
            scale_x_continuous(breaks = seq(-0.05, 0.05, 0.05), limits = lims) +
            scale_y_continuous(breaks = seq(-0.05, 0.05, 0.05), limits = lims)
    }
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
# ps

lgd <- g_legend(p)
p   <- p + theme(legend.position = "none")
ps[[i]] <- p

# bottom <- arrangeGrob(,
#                       lgd)
xtitle <- textGrob(expression("Weighted Whittaker with dynamic " * lambda * " (wWHd)"),
                   gp=gpar(fontsize=14, fontface = "bold"))
g <- arrangeGrob(grobs = ps, nrow = 1, widths = c(1, 1, 1, 1.07),
                 bottom = xtitle,
                 left = textGrob("Other methods", gp=gpar(fontsize=14, fontface = "bold"), rot = 90))
g <- arrangeGrob(g, bottom = lgd)
file <- sprintf("Figure10_compare_with_other_meth_%s_%s_%s.pdf", prefix, typename, str_iter)

width = 10; height = 8
write_fig(g, gsub(".pdf", ".tif", file), width, height, T, res = 300)
# write_fig(g, gsub(".pdf", ".svg", file), width, height, T, res = 300)
# write_fig(g, file, width, height, T)

nrow <- 1e5
ncol <- 100
x <- matrix(rnorm(nrow * ncol), nrow = nrow, ncol = ncol)
v <- data.table(x)

rbenchmark::benchmark(
    q <- colQuantiles(x, probs = probs),
    q_0  <- apply(x, 2, FUN = quantile, probs = probs) %>% t(),
    q_dt <- apply(v, 2, quantile,probs =c(.1,.9,.5),na.rm=TRUE) %>% t(),
    replications = 10
)

# lattice version ---------------------------------------------------------

strip.math <- function(which.given, which.panel, var.name,
                       factor.levels, ...) {
    vn <- var.name[which.given]
    fl <- factor.levels
    expr <- paste(vn, "==", fl, collapse = ",")

    expr <- paste(fl, collapse = ",")
    expr <- paste("expression(", expr, ")", sep = "")

    # expr <- paste("expression(", fl, ")", sep = "")
    # browser()
    fl <- eval(parse(text = expr))
    strip.default(which.given, which.panel, vn, fl, ...)
}

# xyplot( value ~ wWH|index, data = pdat_i,
#         panel = function(x, y, ...){
#             panel.smoothScatter(x, y, ...)
#             panel.abline(a = 0, b = 1, col = "red", size = 1)
#             panel.text(0.5, 0.5, "R2 = 0.5, RMSE = 0.9")
#         },
#         strip = strip.math,
#         scales=list(relation="free"), as.table = T)

