library(Ipaper)
library(magrittr)
library(plyr)
library(data.table)
source('test/07_whit/whit_lambda/smooth_whit_lambda.R')

lst <- get_sbatch("result/whit_lambda/whit2_grp1/")
# lst <- get_sbatch("result/whit_lambda/wBisquare_grp1/")

I <- sapply(lst, length) %>% {which(. <= 1 | names(.) == "")} %>%
    sort() #%T>% print
lst_good <- lst[-I]


df <- llply(lst_good, function(x) {
    data.table(x$coef, lambda = x$lambda, yearid = 1:length(x$lambda))
}, .progress = "text")
d <- melt_list(df, "site")
setkeyv(d, c("site", "yearid"))
d[, cv := mean/sd]
d %<>% na.omit()

## get formula now

st <- sf::read_sf("F:/Github/phenology/phenofit/data_test/whit_lambda/shp/st_1e3_gee.shp") %>%
    data.table() %>% .[, 1:3]
st$site %<>% as.character()
st$IGBPcode_0 %<>% factor(levels = 1:16, labels = IGBPnames_006[1:16])


temp <- merge(d, st[, c(1, 3)], by = "site")
source("test/07_whit/whit_lambda/main_lambda.R")
colnames(temp)[9] <- "IGBP"
temp$lambda %<>% log10()

# temp$lambda <- 10^(temp$lambda)
a <- select_model(temp, index = NULL, IsPlot = T, file = "lambda_TSM.png", type = "coef",
                  optim = T, robust = F)


# temp[lambda > 100, lambda := 100]
# ggplot(temp, aes(IGBPcode_0, lambda, color = IGBPcode_0)) + geom_boxplot()
