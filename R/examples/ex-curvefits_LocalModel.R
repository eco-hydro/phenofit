library(phenofit)
library(Ipaper)

data("CA_NS6")
d = CA_NS6

nptperyear <- 23
INPUT <- check_input(d$t, d$y, d$w, QC_flag = d$QC_flag,
     nptperyear = nptperyear, south = FALSE,
     maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
# plot_input(INPUT)

# Rough fitting and growing season dividing
wFUN <- "wTSM"
brks2 <- season_mov(INPUT,
    options = list(
        rFUN = "smooth_wWHIT", wFUN = wFUN,
        r_min = 0.05, ypeak_min = 0.05,
        lambda = 10,
        verbose = FALSE
    ))
# plot_season(INPUT, brks2, d)

# Fine fitting
fits <- curvefits_LocalModel(
    INPUT, brks2,
    options = list(
        methods = c("AG", "Beck", "Elmore", "Zhang", "Gu"), #,"klos", "Gu"
        wFUN = wFUN,
        nextend = 2, maxExtendMonth = 2, minExtendMonth = 1, minPercValid = 0.2
    ),
    constrain = TRUE
)
# merge local model function into global model function
fits_merged = merge_LocalModels(fits) 

## Visualization ---------------------------------------------------------------
\dontrun{
l_fitting = map(fits %>% guess_names, get_fitting) #%>% melt_list("period")

d_merged = get_fitting(fits_merged[[2]]) %>% cbind(type = "Merged")
d_raw = l_fitting[2:4] %>% set_names(c("Left", "Central", "Right")) %>%
    melt_list("type")
d_obs = d_raw[, .(t, y, QC_flag)] %>% unique()
d_fit = rbind(d_merged, d_raw)[meth == "Zhang"]

levs = c("Left", "Central", "Right", "Merged")
levs_new = glue("({letters[1:4]}) {levs}") %>% as.character()
d_fit$type %<>% factor(levs, levs_new)

p = ggplot(d_obs, aes(t, y)) +
    geom_point() +
    geom_line(data = d_fit, aes(t, ziter2, color = type)) +
    facet_wrap(~type) +
    labs(x = "Date", y = "EVI") +
    scale_x_date(date_labels = "%b %Y", expand = c(1, 1)*0.08) +
    theme_bw(base_size = 13) +
    theme(legend.position = "none",
          strip.text = element_text(size = 14))
p
}
