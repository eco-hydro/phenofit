#' before regression, scale x to 0-1
#' Also return sd, and mean, for the future reconstructure
zscore   <- function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
f_sample <- function(x) sort(sample(x, length(x)*0.1))
self <- function(x, ...) x

#' @param ... additional parameters will be passed to select_model (i.e. optim)
get_coef <- function(lst, ...){
    if (is.list(lst)){
        res <- llply(lst, select_model, d = stat, .progress = "text", ...) %>%
            ldply(function(x) x, .id = "IGBP") %>% data.table() # coef of diff IGBP
    }else{
        res <- select_model(d = stat, lst, ...) %>% data.table()
        res <- res[, IGBP := "0"] %>% reorder_name("IGBP")
    }
    res # return
}

znorm <- function(d){
    vars   <- c("mean", "sd", "cv", "skewness", "kurtosis")

    sd   <- sapply(d[, ..vars],   sd, na.rm = T)
    mean <- sapply(d[, ..vars], mean, na.rm = T)
    newd <- d[, (vars) := lapply(.SD, zscore), .SDcols = vars]
    list(d = newd, mean = mean, sd = sd)
}

select_model <- function(d, index = NULL, IsPlot = F, file = "a.png", type = "coef",
                         optim = T, robust = F){
    # robust   <- FALSE
    lm_fun   <- ifelse(robust, lmRob, lm)
    step_fun <- ifelse(robust, step.lmRob, stepAIC)

    if(!is.null(index)) d <- d[index, ]

    l     <- lm_fun(lambda~(mean+sd+cv+skewness+kurtosis), d, na.action = na.exclude) #
    if (optim){
        l_opt <- step_fun(l, trace = 0)
    }else{
        l_opt <- l
    }
    d$fitted <- predict(l_opt)

    if (IsPlot){
        print(summary(l_opt))
        # file <- "b.png"
        CairoPNG(file, 10, 7, units = "in", dpi = 300)
        p <- ggplot(d, aes(lambda, fitted)) + #, color = IGBP
            geom_point(alpha = 0.07) + #, color = "blue"
            # geom_density_2d() +
            geom_abline(slope = 1, intercept = 0, color = "red", size = 1) +
            facet_wrap(~IGBP)
        print(p)
        dev.off()
        file.show(file)
    }

    if (type == "coef"){
        tidy(l_opt)
    }else{
        glance(l_opt)
    }
}
# set.seed(42); sample(1:100, 10)
# l    <- glm(lambda~kurtosis+mean+sd+skewness+I(sd/mean), stat,
#             family = gaussian(link = "logit"),
#             na.action = na.exclude)


# test function
test <- function(){
    ggplot(x, aes(as.factor(IGBP), lambda)) +
        geom_violin() +
        geom_boxplot(width =0.2) +
        geom_jitter(width = 0.2, alpha = 0.01)

    ## 01. construct connection
    # list(glance = glance(slm_model), tidy = tidy(slm_model))
    # lms %<>% transpose()

    melt(stat[IGBP == 1], measure.vars = vars) %>% ggplot(aes(value, lambda)) +
        geom_point(alpha = 0.05) +
        # geom_smooth(method = "rlm", formula = 10^y ~ x) +
        geom_smooth(method = "lm", color = "red", formula = y~poly(x, 2)) +
        facet_wrap(~variable, scales = "free")
}
