#' before regression, scale x to 0-1
#' Also return sd, and mean, for the future reconstructure
zscore <- function(x) UseMethod("zscore")

zscore   <- function(x) {
    x2 <- x[!is.na(x)]
    mean <- mean(x2)
    sd <- sd(x2)
    return(list(data = (x - mean)/sd,
        mean = mean, sd = sd))
}

zscore.data.frame <- function(d){
    res <- llply(d, zscore) %>% transpose()
    res$data %<>% do.call(cbind, .) %>% set_colnames(colnames(d))
    res
}

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

load_data <- function(indir = "Y:/Github/phenofit/result"){
    ## 0.1 load whittaker lambdas from cluster
    lst <- get_slurm_out(indir)
    group = is.list(lst[[1]][[1]])

    if (group){ lst <- do.call(c, lst) }
    lst %<>% rm_empty() %>% transpose()

    I_del <- which(sapply(lst$lambda, length) == 0)
    x <- llply(lst, function(x){
        unlist(x[-I_del])
    }) %>% do.call(cbind, .) %>% data.table()

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
    stat
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

    info <- d[, .(R2 = GOF(lambda, fitted)[["R2"]]), IGBP]
    # info$R2 %<>% round(3)
    info[, label := sprintf("R^2=='%4.2f'", R2)]
    if (IsPlot){
        print(summary(l_opt))
        # file <- "b.png"
        CairoPNG(file, 10, 7, units = "in", dpi = 300)
        p <- ggplot(d, aes(lambda, fitted)) + #, color = IGBP
            geom_point(alpha = 0.05) + #, color = "blue"
            # geom_density_2d() +
            geom_abline(slope = 1, intercept = 0, color = "red", size = 1) +
            facet_wrap(~IGBP) +
            geom_text(data = info, aes(x = Inf, y = Inf, label = label),
                      vjust = 1.2, hjust = 1, parse = T, size = 5) +
            labs(x = expression("Optimized "* log(lambda)),
                 y = expression("Simulated "* log(lambda))) +
            theme_grey(base_size = 15)
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
