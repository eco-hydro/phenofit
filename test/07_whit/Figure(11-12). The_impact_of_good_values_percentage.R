library(scales)

# FUNCTIONS ---------------------------------------------------------------

which_max <- function(x){
    if (all(is.na(x))){
        NA
    } else {
        which.max(x)
    }
}

which_min <- function(x){
    if (all(is.na(x))){
        NA
    } else {
        which.min(x)
    }
}

# load df_info from Fig10_.R
df = df_info[type == "all", ]

levels  <- seq(0, 1, 0.01)
xmid0   <- c((levels[-1] + levels[-length(levels)])/2, 1)

###############################################################################
## 03. Fig10. The influence of good values percentage

df[, grp_perc := findInterval(perc_good, levels, rightmost.closed = T)]
df[, xmid := xmid0[grp_perc]]

d <- df[, .(site, meth, type, iters, RMSE, R2, Bias, Rg = Rg_norm_by_obs, grp_perc, xmid)] %>% #, "NSE"
    melt(id.vars = c("site", "meth", "type", "iters", "grp_perc", "xmid"), variable.name = "index")

# ggplot(d[index == "R2"], aes(xmid,value)) + geom_point() + geom_density2d() +
#     facet_wrap(~meth)
# d <- df[meth == "wWH" & iter == "iter1"]

alphas <- c(.05, .1, .25, .5) %>% set_names(., .)
d_envelope <- llply(alphas, function(alpha){
    # print(alpha)
    expr <- substitute(quote(res <- ddply_dt(d, .(quantile_envelope(value, alpha)), .(meth, index, iters, grp_perc))),
               list(alpha = alpha))
    # print(eval(expr))
    eval(eval(expr))
})

d_enve <- d_envelope %>% melt_list("alpha")
d_enve[, xmid := xmid0[grp_perc]]

# hue_pal()(4) %>% show_col()
cols <- scales::hue_pal()(4)
colors <- c(alpha(cols[1], 0.5), cols[2:3], "white")

pdat <- d_enve[meth == "wWH" & iters == "iter2"]

indice <- c("R2", "Bias", "RMSE", "Rg")
labels <- c("bold((a)*' '*R^2)", "bold((b)*' '*Bias)",
            "bold((c)*' '*RMSE)", "bold((d)*' '*Roughness)")
# pdat <- d[meth != "whit_R" & iter == "iter2" & index %in% indice]
pdat$index %<>% factor(indice, labels)

# sub_id <- geom_text(data = data.table(index = factor(labels, labels),
#                                        label = sprintf("(%s)", letters[1:length(labels)])),
#                     aes(label = label),
#                     x = -Inf, y = Inf, vjust = 1, hjust = 0)
d_vline <- data.table(index = factor(labels[1:4], labels),
                      x = c(.25, .25, .25, .25))
v_line <- geom_vline(data = d_vline ,
                     aes(xintercept = x), color = "black", linetype = 2, size = 1)
v_line_lab <- geom_text(data = d_vline, hjust = -0.1, vjust = -0.4,
                        aes(x, y = -Inf,label = sprintf("%2d%%", x*100)))

p <- ggplot(pdat, aes(xmid, ymin)) +
    # geom_point() + geom_density2d() +
    geom_ribbon(aes(x = xmid, ymin = ymin, ymax = ymax, fill = alpha)) +
    facet_wrap(~index, scales = "free", labeller=label_parsed) +
    # sub_id +
    geom_line(data = pdat[alpha == 0.5, ], color ="white") +
    v_line + v_line_lab +
    # geom_vline(xintercept = c(0.25), color= "blue", size = 1, linetype = 2) +
    scale_fill_manual(values = colors, labels = c("  5% and 95%", "10% and 90%", "25% and 75%", "50%")) +
    guides(fill = guide_legend(override.aes = list(color = "black"))) +
    labs(x = "The percentage of good points (%)") +
    scale_x_continuous(labels = function(x) sprintf("%d", x*100)) +
    scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
    theme_gray(base_size = 14) +
    theme(
        legend.title = element_blank(),
        legend.position = c(0.98, 0.985),
        legend.justification = c(1, 1),
        legend.spacing.x = unit(0.1, 'cm'),
        strip.text = element_text(face = "bold", size = 13, margin = margin(2, 0, 1, 0, "pt")),  # not work for meth expr
        axis.text = element_text(size = 12),
        axis.title.x = element_text(face = "bold", size = 13)
    )+
    ylab(NULL)


file <- "Figure11 good vlaues pertenage impact.pdf"
width = 9.5; height = 6
write_fig(p, gsub(".pdf", ".svg", file), width, height, F, res = 300)
write_fig(p, gsub(".pdf", ".tif", file), width, height, T, res = 300)
# write_fig(p, file, width, height, T)

# second try --------------------------------------------------------------

methods <- c("AG", "ZHANG", "smooth_wHANTS", "smooth_wSG", "wWHd") %>% factor(., .)

d_perc <- d %>% dcast(site+iters+grp_perc+xmid+index+type~meth, value.var = "value") %>%
    .[iters == "iter2", ]

I    <- d_perc[index != "R2", .(AG, ZHANG, smooth_wHANTS, smooth_wSG, wWHd = wWH)] %>%
    as.matrix %>% aaply(1, which_min) %>% as.numeric()
I_r2 <- d_perc[index == "R2", .(AG, ZHANG, smooth_wHANTS, smooth_wSG, wWHd = wWH)] %>%
    as.matrix %>% aaply(1, which_max) %>% as.numeric()

x    <- cbind(d_perc[index != "R2", 1:5], dominant = methods[I])
x_r2 <- cbind(d_perc[index == "R2", 1:5], dominant = methods[I_r2])


x <- rbind(x, x_r2)
x$index %<>% factor(indice, indice_label)
x[, grp := cut(xmid*100, seq(0, 1, 0.1)*100)]
pdat <- x[!is.na(dominant), .(count = .N), .(grp, dominant, index)]
d_count_grp <- pdat[, .(sum = sum(count)), .(grp, index)]
pdat <- merge(pdat, d_count_grp)
pdat[, perc := count/sum]

p <- ggplot(pdat,
            aes(grp, perc*100, fill = dominant)) +
    geom_bar(stat = "identity") +
    labs(x = "The percentage of good points (%)",
         y = "The pectenage of dominant times (%)") +
    theme_gray(base_size = 14) +
    theme(legend.title = element_blank(),
          legend.position = "bottom",
          legend.spacing.x = unit(0.2, 'cm'),
          legend.margin = margin(-2),
          strip.text = element_text( margin = margin(1, 1, 1, 1, "pt")*1.5, size = 14, face = "bold"),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    ) +
    # coord_flip() +
    facet_wrap(~index, labeller = label_parsed)

file <- "Figure12 dominant times across good_values_percentage.pdf"
width = 9; height = 6
write_fig(p, gsub(".pdf", ".svg", file), width, height, F, res = 300)
write_fig(p, gsub(".pdf", ".tif", file), width, height, T, res = 300)
