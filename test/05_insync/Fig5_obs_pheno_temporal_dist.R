

# stations212 <- fread("F:/Github/MATLAB/PML/MATLAB/LAI/flux-212.txt")
# df <- merge(df_obs, stations212)

# d_obs.sos <- select(df_obs, matches("sos|UD|SD|Greenup|Maturity"))
# d_obs.eos <- select(df_obs, matches("eos|DD|RD|Senescence|Dormancy"))
pdat <- gather(df_obs, phenophase, value, -flag, -origin, -meth, -site, -IGBP, -lat) %>% data.table()

t1 <- pdat[, .(median = median(value, na.rm = T),
               sd = sd(value, na.rm = T)), .(phenophase)] %>%
    plyr::mutate(label = sprintf("%.1f±%.1f", median, sd))
t1$phenophase %<>% fix_level()
t1 <- t1[order(phenophase), ]

pdat$phenophase %<>% fix_level()
pdat <- pdat[!is.na(phenophase), ]




p <- ggplot(pdat, aes(phenophase, value)) +
    # geom_violin() +
    stat_summary(fun.data = box_qtl, geom = "errorbar", width = 0.5) +
    geom_boxplot2(aes(), coef = 0, width = 0.8, notch = T) +
    # geom_vline(xintercept = c(8.5, 11.5, 15.5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(size = 0.5)
        # panel.grid.major.x = element_blank()
        ) +
    labs(x = NULL, y = "Doy (day of year)", color = 'Phenophase')

doy_trans <- function(x){
    date <- as.Date(sprintf("2010-%03d", as.integer(x)), "%Y-%j")
    format(date, '%m/01')
}
brks <- ymd(sprintf("2010-%02d-01", seq(2, 12, 2))) %>% yday()
p_final <- p + scale_y_continuous("Date", breaks = brks, label = doy_trans,
    sec.axis = sec_axis(~., breaks = brks, name = "DOY (day of year)")) +
    geom_vline(xintercept = c(4, 7, 10, 12, 15) + 0.5, linetype = 2, size = 0.5, col = "red")

# df_text <- data_frame(y = rep(Inf, 6), x = c(1, 5, 8, 11, 13, 16), label = letters[1:6])
p_final

cairo_pdf("Fig6_GPPobs_phenology_temporal_boxplot_dist.pdf", 9, 5.5)
print(p_final)
dev.off()

file.show("Fig6_GPPobs_phenology_temporal_boxplot_dist.pdf")
