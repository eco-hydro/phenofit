# res <- list()
# for (i in seq_along(sites)){
#     runningId(i)
#     res[[i]] <- optim_lambda(sites[i])
# }

source('test/stable/load_pkgs.R')

lst <- get_slurm_out("result")
a <- do.call(c, lst) %>% rm_empty()

x <- transpose(a) %>% llply(unlist) %>% do.call(cbind, .) %>% data.table()
x$lambda %<>% log10()
x <- x[!is.na(lambda), ]
x[, grp := (0:(.N-1)), site]

# ~llply(.x[c("site", "IGBP", "lambda")], first, default = NA)
# x <- map_df(a, ~.x) %>% data.table()

fwrite(x, "st1e4_lambda2.csv")
## statistics
df[, grp := floor(1:.N/(23*3)), site]
df <- df[grp <= 5]


# info <- ddply(dt, .(site), function(d){
#     y <- d$y[1:(23*3)]
#     c(mean = mean(y, na.rm = T),
#       sd = sd(y, na.rm = T),
#       kurtosis = kurtosis(y, na.rm = T, type = 2),
#       skewness = skewness(y, na.rm = T, type = 2))
# })
info <- df[!is.na(y), .(mean = mean(y),
               sd = sd(y),
               kurtosis = kurtosis(y, type = 2),
               skewness = skewness(y, type = 2)), .(site, grp)]
# info[, cv := sd/mean]

ggplot(x, aes(as.factor(IGBP), lambda)) +
    geom_violin() +
    geom_boxplot(width =0.2) +
    geom_jitter(width = 0.2, alpha = 0.01)


## 01. construct connection
stat <- merge(data.table(x), info, by = c("site", "grp"))
slm(stat)


lms <- dlply(stat, .(IGBP), function(d){
    lm_model    <- rlm(lambda~kurtosis+mean+sd+skewness+I(sd/mean), d, na.action = na.exclude) #

    slm_model   <- stepAIC(lm_model)
    stat$fitted <- predict(slm_model)

    list(glance = glance(slm_model), tidy = tidy(slm_model))
})
lms %<>% transpose()


#' before regression, scale x to 0-1
#' Also return sd, and mean, for the future reconstructure
znorm <- function(d){
    zscore <- function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
    vars   <- c("mean", "sd", "cv", "skewness", "kurtosis")

    sd   <- sapply(d[, ..vars],   sd, na.rm = T)
    mean <- sapply(d[, ..vars], mean, na.rm = T)
    newd <- d[, (vars) := lapply(.SD, zscore), .SDcols = vars]
    list(d = newd, mean = mean, sd = sd)
}

stat[, cv := mean/sd]
stat_s <- znorm(stat) # scaled value


lst  <- split(1:nrow(stat), stat$IGBP) %>% set_names(NULL)
f_sample <- function(x) sort(sample(x, length(x)*0.1))
# resample 10 times
times <- 1:20 %>% set_names(., .)
I_lst <- llply(times, function(i){
    llply(lst, f_sample) %>% do.call(c, .)
})

res <- llply(I_lst, select_model, .progress = "text") %>%
    ldply(function(x) x)
res
table(res$term)


select_model <- function(index){
    robust <- FALSE
    lm_fun   <- ifelse(robust, lmRob, lm)
    step_fun <- ifelse(robust, step.lmRob, stepAIC)

    d <- stat_s$d[index, ]
    l <- lm_fun(lambda~(kurtosis+mean+sd+skewness+I(sd/mean)), d, na.action = na.exclude) #
    l_opt <- step_fun(l, trace = 0)
    tidy(l_opt) # return
}
# set.seed(42); sample(1:100, 10)

#' stepwise regression for whittaker lambda relationship constructure
slm  <- function(stat){
    # robust
    stat$fitted <- predict(l_opt)
    l    <- glm(lambda~kurtosis+mean+sd+skewness+I(sd/mean), stat,
                       family = gaussian(link = "logit"),
                       na.action = na.exclude)
    print(summary(l_opt))
    d2 <- df[, .(lambda2 = init_lambda(y)), .(site, grp)]

    file <- "a.png"
    CairoPNG(file, 10, 7, units = "in", dpi = 300)
    p <- ggplot(stat, aes(lambda, fitted)) + #, color = IGBP
        geom_point(alpha = 0.07) + #, color = "blue"
        # geom_density_2d() +
        geom_abline(slope = 1, intercept = 0, color = "red", size = 1) +
        facet_wrap(~IGBP)
    print(p)
    dev.off()
    file.show(file)
}
slm(stat)

dt[, .(lambda2 = init_lambda(y)), .(site, grp)]
## 02. visualization
xs <- melt(x, c("site", "IGBP", "lambda"))
ggplot(xs, aes(value, lambda)) + geom_point(alpha = 0.05) + geom_smooth() +
    facet_wrap(~variable, scales = "free")
# facet_grid(IGBP~variable, scales = "free")

ggplot(stat, aes(lambda, fitted, color = IGBP)) + geom_point()


