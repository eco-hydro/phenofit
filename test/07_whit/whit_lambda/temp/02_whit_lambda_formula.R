source("test/load_pkgs.R")

stat    <- load_data(indir = "Y:/Github/phenofit/result")
stat_3y <- load_data(indir = "Y:/Github/phenofit/result/grp/")


stat$group <- F
stat_3y$group <- T

df <- rbind(stat, stat_grp)
# test lambda difference between 3y and 16year
aov(lambda ~ IGBP + group, df) %>% summary()

ggplot(df, aes(as.factor(IGBP), lambda, color = group)) +
    geom_boxplot()
# stat_s <- znorm(stat)$d # scaled value

## 1.1 check over all gof
# select_model(stat_s, IsPlot = T)
select_model(stat, IsPlot = T, file = "whit_lambda_IGBP.png")
select_model(stat_3y, IsPlot = T, file = "whit_lambda_IGBP_3y.png")


tuk <- glht(lm(lambda ~ IGBP, stat), linfct = mcp(IGBP = "Tukey"))
plot(cld(tuk, level = .05), col = "lightgrey")

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
