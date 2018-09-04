
###############################################################################
## 03. Fig10. The influence of good values percentage
levels <- seq(0, 1, 0.01)
xmid   <- c((levels[-1] + levels[-length(levels)])/2, 1)

df[, grp_perc := findInterval(perc_good, levels, rightmost.closed = T)]
df[, xmid := xmid[grp_perc]]

d <- df[, .(site, meth, type, iter, RMSE, R2, Bias, Rg, grp_perc, xmid)] %>% #, "NSE"
    melt(id.vars = c("site", "meth", "type", "iter", "grp_perc", "xmid"), variable.name = "index")

ggplot(d[index == "R2"], aes(xmid,value)) + geom_point() + geom_density2d() +
    facet_wrap(~meth)
# d <- df[meth == "wWH" & iter == "iter1"]

alphas <- c(.05, .1, .25, .5) %>% set_names(., .)
d_envelope <- llply(alphas, function(alpha){
    # print(alpha)
    expr <- substitute(quote(res <- ddply_dt(d, .(quantile_envelope(value, alpha)), .(meth, index, iter, grp_perc))),
               list(alpha = alpha))
    # print(eval(expr))
    eval(eval(expr))
})

d_enve <- d_envelope %>% melt_list("alpha")
d_enve[, xmid := xmid[grp_perc]]

# hue_pal()(4) %>% show_col()
cols <- hue_pal()(4)
colors <- c(trans_col(cols[1], 0.5), cols[2:3], "white")

pdat <- d_enve[meth == "wWH" & iter == "iter1"]

indice <- c("R2", "Bias", "RMSE", "Rg")
labels <- c("bold((a)*' '*R^2)", "bold((b)*' '*Bias)",
            "bold((c)*' '*RMSE)", "bold((d)*' '*Roughness)")
# pdat <- d[meth != "whit_R" & iter == "iter2" & index %in% indice]
pdat$index %<>% factor(indice, labels)

# sub_id <- geom_text(data = data.table(index = factor(labels, labels),
#                                        label = sprintf("(%s)", letters[1:length(labels)])),
#                     aes(label = label),
#                     x = -Inf, y = Inf, vjust = 1, hjust = 0)
d_vline <- data.table(index = factor(labels[1:3], labels),
                      x = c(.3, .25, .25))
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
    labs(x = "The percentage of good values (%)") +
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


file <- "Fig13_good_values_percentage_impact.pdf"
file_tiff <- gsub(".pdf", ".tif", file)
width = 9.5; height = 6
write_fig(p, file_tiff, width, height, T, res = 300)
# write_fig(p, file, width, height, T)
# 