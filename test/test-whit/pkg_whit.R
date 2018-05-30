library(phenofit)
library(Cairo)
library(tidyverse)
library(e1071)

load("Y:/R/phenofit/data/phenofit_MultipleINPUT_flux212.rda")
stations212 <- fread("C:/Users/kon055/Google Drive/Github/data/phenology/station/flux-212.txt")

dt <- lst$MOD13A1
dt <- merge(dt, stations212[, .(site, lat, IGBP)], by = "site")

nptperyear = 23
sites <- unique(dt$site)

IsPlot = TRUE

x <- map_df(res, ~.x[c("site", "IGBP", "lambda")]) #, "lat"
x$lambda %<>% log10()

info <- ddply(dt, .(site), function(d){
    y <- d$y[1:(23*3)]
    c(mean = mean(y, na.rm = T),
     sd = sd(y, na.rm = T),
     kurtosis = kurtosis(y, na.rm = T, type = 2),
     skewness = skewness(y, na.rm = T, type = 2))
})
info <- dt[, .(mean = mean(y, na.rm = T),
               sd = sd(y, na.rm = T),
               kurtosis = kurtosis(y, na.rm = T, type = 2),
               skewness = skewness(y, na.rm = T, type = 2)), .(site)]
# info2 <- dt[, kurtosis(x, w), site]
stat <- merge(data.table(x), info, by = "site")

# ggplot(x, aes(IGBP, lambda)) + geom_boxplot() +
#     geom_jitter(width = 0.2)

xs <- melt(stat, c("site", "IGBP", "lambda"))

ggplot(xs, aes(value, lambda)) + geom_point() + geom_smooth() +
    facet_wrap(~variable, scales = "free")
    facet_grid(IGBP~variable, scales = "free")


summary(slm1)
# slm1$anova

par(mfrow = c(4, 4))
lms <- dlply(x, .(IGBP), function(d){
    lm_model <- lm(lambda~kurtosis+mean+sd,data = d)

    plot(lm_model$fitted.values, d$lambda); grid()
    title(sprintf("%s", d$IGBP[1]))
    abline(a = 0, b = 1, col = "red", lwd = 1.5)
    lm_model
})

for (i in seq_along(lms)){
    lm_model = lms[i]
    plot(lm1$fitted.values, x$lambda); grid()
}
fitted_models = x %>% group_by(IGBP) %>%
    do(model = lm(lambda~kurtosis+mean+sd,data = .))

summary(lm1)
par(mfrow = c(2, 2))
plot(lm1)

color=rgb(0,0,0,alpha=0.3)
# dplot = data_table(fitted = slm1$fitted.values, value = x$lambda)
plot(slm1$fitted.values, x$lambda, col = color); grid()
abline(a = 0, b = 1, col = "red", lwd = 1.5)

d_plot <- data.table(fitted = lm1$fitted.values, x[, .(lambda, site, IGBP)])
d_label <- d_plot[abs(fitted - lambda) >= 0.8]

ggplot(d_plot, aes(lambda, fitted)) +
    geom_point() +
    geom_label_repel(data = d_label, aes(label = site)) +
    geom_abline(slope = 1)
