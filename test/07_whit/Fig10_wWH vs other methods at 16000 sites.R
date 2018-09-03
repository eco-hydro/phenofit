## second solution
library(grid)
library(gridExtra)

label_value <- function (labels, multi_line = TRUE, sep = "*'~'*")
{
    out <- do.call("Map", c(list(paste, sep = sep), labels))
    out <- list(unname(unlist(out)))

    lapply(out, function(values) {
        values <- paste0("list(", values, ")")
        lapply(values, function(expr) c(parse(text = expr)))
    })
}


d2 <- d[iter == "iter2"] %>% {dcast(., site+IGBPname+iter+index~meth, value.var = "value")[, c(1:5, 8, 10:13)]} %>%
    melt(c("site", "IGBPname", "iter", "index", "wWH"), variable.name = "meth")
d2[, kind:=0]

d_diff <- data.table(index = levels(d2$index) %>% factor(., .),
                     dmax = c(.05, .02, 0.02, .005))
d2 %<>% merge(d_diff)
d2[, label:= sprintf("%s~%s", meth, index)]

d2[wWH - value > dmax, kind := 1]
d2[wWH - value < -dmax, kind := -1]
d2$kind %<>% factor(levels = c(1, 0, -1),
                    labels = c("Bigger", "Similar", "Smaller"))

colors <- scales::hue_pal()(3)
colors <- c(colors[2], "grey70", colors[1])

ps <- list()
for (i in seq_along(indice)){
    indexi = indice[i]
    pdat <- d2[index == indexi]

    d_dominant <- ddply_dt(pdat, .(table(kind)), .(meth, index, label))
    d_dominant[, text := sprintf("%2d vs %2d", Bigger, Smaller)]
    d_dominant[, `:=`(
        text_big   = sprintf("bold('%s' * phantom(' vs ' * '%s'))", Bigger, Smaller),
        text_vs    = sprintf("bold(phantom('%s') * ' vs ' * phantom('%s'))", Bigger, Smaller),
        text_small = sprintf("bold(phantom('%s vs ') * '%s')", Bigger, Smaller))]

    text_size = 5
    p <-
        ggplot(pdat, aes(wWH, value, color = kind)) +
        geom_point(data = pdat[kind != "Similar"], alpha = 0.4) +
        geom_point(data = pdat[kind == "Similar"], alpha = 0.2) +
        # facet_wrap(~label, scales = "free") +
        facet_grid(meth~index, scales = "free", labeller = label_parsed) + #, ncol = 4, labeller = label_value
        geom_abline(slope = 1, color ="red", size = 0.5) +
        scale_color_manual(values = colors) +

        geom_text(data = d_dominant, aes(label = text_big, color = NULL), parse = T,
                  color = colors[1], fontface = "bold", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, show.legend = F, size = text_size) +
        geom_text(data = d_dominant, aes(label = text_vs, color = NULL), parse = T,
                  color = "black", fontface = "bold", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, show.legend = F, size = text_size) +
        geom_text(data = d_dominant, aes(label = text_small, color = NULL), parse = T,
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


lgd <- g_legend(p)
p   <- p + theme(legend.position = "none")
ps[[i]] <- p

# bottom <- arrangeGrob(,
#                       lgd)
xtitle <- textGrob("wWH", gp=gpar(fontsize=14, fontface = "bold"))
g <- arrangeGrob(grobs = ps, nrow = 1, widths = c(1, 1, 1, 1.1),
                 bottom = xtitle,
                 left = textGrob("Other methods", gp=gpar(fontsize=14, fontface = "bold"), rot = 90))
g <- arrangeGrob(g, bottom = lgd)
file <- "Fig8_compare_with_other_methods_v2.pdf"
file_tiff <- gsub(".pdf", ".tif", file)
width = 11; height = 9
write_fig(g, file_tiff, width, height, T, res = 300)
# write_fig(g, file, width, height, T)
