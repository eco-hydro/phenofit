
d_obs <- fread("file:///Z:/dailyforing_93sites/INPUTS/PML_inputs_8d_dynamicLC_006_CO2_v3.csv")
d_obs$date %<>% ymd()

# 2. zhang yao, 2017, scientific data, VPMGPP, 8day -------------------------------
files <- dir("C:/Users/kon055/Desktop/VPMGPP", full.names = T) %>%
    set_names(gsub(".csv","", basename(.)))
d_vpm <- ldply(files, fread, .id = "site")[, c(1, 3, 4)] %>%
    set_names(c("site", "date", "GPP")) %>% as.data.table()
d_vpm[, date := as.Date(date, "%Y%j")]
d_vpm %<>% set_colnames(c("site", "date", "GPP_vpm"))
d_vpm$w <- 1

df <- merge(d_obs, d_vpm)
fwrite(df, "flux112_GPP_obs&GPP_vpm.csv")

p <- ggplot(mpg, aes(class, hwy))

# how to control not show outliers in boxplot
p +
    stat_summary(fun.data = whits, geom = "errorbar", width = 0.5) +
    geom_boxplot(notch = TRUE, outlier.shape = NA, coef = 0)
    # stat_summary(fun.y=mean, geom="point", size=2) +

