rm(list = ls())
source('test/stable/load_pkgs.R')
source('test/07_whit/whit_lambda/smooth_whit_lambda.R')
source("test/07_whit/whit_lambda/main_lambda.R")


check_lambda <- function(lambda){
    n <- length(lambda)
    lambda <- sort(lambda)
    if (n >= 10){
        lambda <- lambda[3:(n-2)]
    } else if (n >= 5){
        lambda <- lambda[2:(n-1)]
    }
    c(mean = mean(lambda), sd = sd(lambda))
}


library(sf)
indir <- "Y:/Github/phenofit_cluster/"
## 1. station info
st   <- sf::read_sf("data_test/whit_lambda/shp/st_1e3_mask.shp")
coor <- st_geometry(st) %>% do.call(rbind, .) %>% data.table() %>%  set_colnames(c("lon", "lat"))

st <- as.data.table(st)[, 1:3] %>% cbind(coor)

st$site %<>% as.character()
st$IGBPcode %<>% factor(levels = 1:16, labels = IGBPnames_006[1:16])
colnames(st)[3] <- "IGBP"

dirs <- dir("Y:/Github/phenofit_cluster/result/whit_lambda/", full.names = T)[-1] %>%
    set_names(gsub("whit2", "", basename(.)))

## 2. lambda info

res_formula <- list()
for (k in 1:10){ #length(dirs)
    indir <- dirs[k]
    runningId(k, prefix = indir)

    file <- sprintf('Figure/Fig1_lambda_%s.jpg', basename(indir))

    lst <- get_sbatch(indir)
    # lst <- get_sbatch("result/whit_lambda/wBisquare_grp1/")
    I <- sapply(lst, length) %>% {which(. <= 1 | names(.) == "")} %>%
        sort() #%T>% print
    lst_good <- lst[-I]

    d <- llply(lst_good, function(x) {
        data.table(x$coef, lambda = x$lambda, yearid = 1:length(x$lambda))
    }, .progress = "text") %>% melt_list("site") %>% na.omit() %>% # rm NA values
        setkeyv(c("site", "yearid"))

    ## 1. avoid V-curve failure, rm outlier lambda in every site
    #1. rm outlier lambda
    # d$site %<>% as.numeric()
    if (k <= 10){
        d_sd <- ddply_dt(d, .(check_lambda(lambda)), .(site)) %>%
            set_names(c("site", "l_mean", "l_sd"))
        d <- merge(d, d_sd)
        d <- d[lambda < (l_mean + 3*l_sd) & lambda > (l_mean - 3*l_sd), 1:7]
    }
    d[, cv := mean/sd]

    ## get formula now
    temp <- merge(d, st[, c(1, 3:5)], by = "site")
    temp$lambda %<>% log10()

    predictors <- c("mean", "sd", "cv", "skewness", "kurtosis", "lat")
    response   <- c("lambda")
    vars       <- c(predictors, response)

    # dt <- zscore.data.frame(temp[, ..vars])$data %>% as.data.table() %>%
    #     cbind(temp[, .(IGBP)])
    dt <- temp # if you want to use zscore, just comment this line

    # temp$lambda <- 10^(temp$lambda)
    # temp[lambda > 100, lambda := 100]
    # ggplot(temp, aes(IGBPcode_0, lambda, color = IGBPcode_0)) + geom_boxplot()

    n <- nrow(dt)
    res <- list()

    for (i in 1:100){
        runningId(i)
        I <- sample(1:n, n*0.4)
        d <- dt[-I, ]
        res[[i]] <- select_model(d, index = NULL, IsPlot = F, file = file, type = "coef",
                             optim = T, robust = F)
    }
    temp <- select_model(dt, index = NULL, IsPlot = T, file = file, type = "coef",
                             optim = T, robust = F)

    d <- melt_list(res, "Id") %>% data.table()
    info <- d[, .(Id, term, estimate)] %>% dcast(Id~term, value.var = "estimate")
    info <- info[,2:ncol(info)] %>% as.matrix() %>% {
        list(n = aaply(., 2, function(x)sum(!is.na(x))),
             mean = aaply(., 2, mean, na.rm = T),
             sd   = aaply(., 2, sd  , na.rm = T))
    } #%T>% print
    # info$sd/info$mean

    formula <- info$mean[-c(2:3)] %>% {
        c(sprintf("+%.4f", .[[1]]),
          sprintf("%+.4f*b('%s')", ., names(.))[-1])
    } %>% paste(collapse = " ") %>%
        gsub("^[+-]", "var formula = ", .) %>%
        gsub("\\*\\(Intercept\\)", "", .) %T>%
        cat %T>% writeLines("clipboard")
    res_formula[[k]] <- list(info = info, formula = formula)
}

names(res_formula) <- names(dirs)

a <- transpose(res_formula)
coef <- a$info %>% transpose() %>% map(~do.call(rbind, .))

coef$perc <- with(coef, sd/mean*100)
save(res_formula, coef, file = "data_test/lambda_formula.rda")

# var formula = 0.8214 +1.5025*b('mean') -4.0315*b('sd') -0.1018*b('skewness')

# 0.4 part:
# ---------
# $`n`
# cv kurtosis     mean       sd skewness
# 8       94      100      100      100
#
# $mean
# cv     kurtosis         mean           sd     skewness
# -0.003392106 -0.004536728  1.502439739 -4.032131238 -0.101649273
#
# $sd
# cv    kurtosis        mean          sd    skewness
# 0.006084250 0.001470645 0.022139485 0.052780635 0.003400217

get_formula <- function(){
    I <- sample(1:n, n*0.4)
    d <- dt[-I, ]

    source("test/07_whit/whit_lambda/main_lambda.R")
    file <- sprintf("Fig1_lambda_%s_mask.png", grp)
    info <- select_model(d, index = NULL, IsPlot = F, file = file, type = "coef",
                         optim = T, robust = F)
    formula <- info %$% sprintf("%+.4f*%s", estimate, term) %>%
        paste(collapse = " ") %>%
        gsub("^[+-]", "lambda <- ", .) %>%
        gsub("\\*\\(Intercept\\)", "", .) %T>%
        cat %T>% writeLines("clipboard")

    formula_gee <- info %$% sprintf("%+.4f*b('%s')", estimate, term)[-1] %>%
        c(sprintf("+%.4f", info$estimate[1]), .) %>%
        paste(collapse = " ") %>%
        gsub("^[+-]", "var formula = ", .) %>%
        gsub("\\*\\(Intercept\\)", "", .) %T>%
        cat %T>% writeLines("clipboard")
}


# yearly: lambda <- 0.7835 +1.5959*mean -4.0371*sd +0.0048*cv -0.1032*skewness +0.0036*kurtosis
# 4-year: lambda <- 0.8209 +1.5008*mean -4.0286*sd -0.1017*skewness -0.0041*kurtosis
