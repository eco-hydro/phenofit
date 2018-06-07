source('test/stable/load_pkgs.R')

## 0.1 load whittaker lambdas from cluster
lst <- get_slurm_out("Y:/Github/phenofit/result") %>% do.call(c, .) %>% rm_empty()

x <- transpose(lst) %>% llply(unlist) %>% do.call(cbind, .) %>% data.table()
x$lambda %<>% log10()
x <- x[!is.na(lambda), ]
x[, grp := (0:(.N-1)), site]

# fwrite(x, "st1e4_lambda2.csv")
## 02. read point statistics (i.e. mean, sd, cv, ...)
df <- fread('data_whit_lambda_01.csv')
df[, grp := floor(1:.N/(23*3)), site]
df <- df[grp <= 5]

info <- df[!is.na(y), .(mean = mean(y),
                        sd = sd(y),
                        kurtosis = kurtosis(y, type = 2),
                        skewness = skewness(y, type = 2)), .(site, grp)]
info[, cv:=mean/sd]
stat <- merge(data.table(x), info, by = c("site", "grp"))
# stat_s <- znorm(stat)$d # scaled value

## 1.1 check over all gof
# select_model(stat_s, IsPlot = T)
select_model(stat, IsPlot = T)

## 1.2 resample multiple group
lst  <- split(1:nrow(stat), stat$IGBP)# %>% set_names(NULL)
# resample 10 times
times  <- 1:20 %>% set_names(., .)
I_lst  <- llply(times, function(i){ llply(lst, f_sample) %>% unlist })
I_lst2 <- llply(times, function(i){ llply(lst, f_sample) })


## 1.3 get group regression coef
optim <- FALSE
res2 <- llply(I_lst2, get_coef, optim = optim) %>% ldply(self, .id = 'sgrp') %>% data.table()
res  <- llply(I_lst, get_coef, optim = optim, .progress = "text") %>% ldply(self, .id = 'sgrp') %>% data.table()
# res
# table(res$term)
# res[, .(mean = mean(estimate), sd = sd(estimate)), .(term)]
r    <- rbind(res, res2)
coef <- r[, .(sgrp, IGBP, term, estimate)] %>% dcast(sgrp+IGBP~term, value.var = "estimate")

## 2.0 tranform and tidy lm coef
data <- melt.data.table(coef, c("sgrp", "IGBP"))
data <- merge(data[IGBP != "0", ], data[IGBP == "0", c(1, 3:4)], c("sgrp", "variable")) %>%
    melt(measure.vars = c("value.x", "value.y"), variable.name = "type")

# data[, whole := IGBP == "0"]
data[, IGBP  := as.integer(as.character(IGBP))]

## 2.2 aov test
test <- ldply(set_names(1:17, 1:17), function(i){
    d <- data[IGBP == i, ]
    formula <- value ~ variable + type
    # formula <- value ~ type
    aov(formula, d) %>% tidy()
})
a <- test %>% data.table() %>% {.[term == 'type', ]}
a[p.value < 0.05]

################################################################################
# fit <- aov(estimate~IGBP, info[term == 'mean']) %T>% summary()

# dt[, .(lambda2 = init_lambda(y)), .(site, grp)]
# ## 02. visualization
# xs <- melt(x, c("site", "IGBP", "lambda"))
# ggplot(xs, aes(value, lambda)) + geom_point(alpha = 0.05) + geom_smooth() +
#     facet_wrap(~variable, scales = "free")
# # facet_grid(IGBP~variable, scales = "free")
# ggplot(stat, aes(lambda, fitted, color = IGBP)) + geom_point()
