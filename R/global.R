# Other levels will be ignored.
qc_levels <- c("good", "marginal", "snow", "cloud", "aerosol", "shadow")
qc_colors <- c("grey60", "#00BFC4", "#F8766D", "#C77CFF", "#B79F00", "#C77CFF") %>% set_names(qc_levels)
qc_shapes <- c(19, 15, 4, 25, 25, 17) %>% set_names(qc_levels)

date.origin = as.Date("2000-01-01")

# color function of i_th iteration smoothed time-series
iter_colors <- colorRampPalette(c("blue", "green", "red")) #(4) %>% scales::show_col()
