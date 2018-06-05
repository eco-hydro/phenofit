source('test/stable/load_pkgs.R')
# source('test/GEE/V-pack.r')

dt = fread('test/GEE/data/MOD13A1_st_1e3.csv')
dt[, `:=`(y    = EVI/1e4,
          t    = ymd(date),
          w    = qc_summary(dt$SummaryQA))]
dt[, per := sum(!is.na(EVI))/.N, site]
df <- dt[per > 0.3, .(site, y, t, w, IGBPcode)]

sites      <- unique(df$site)
nptperyear <- 23

res <- list()

for (i in seq_along(sites)){
    runningId(i)
    sitename <- sites[i]#; grp = 1
    d    <- df[site == sitename]
    IGBP <- d$IGBPcode[1]

    INPUT <- check_input(d$t, d$y, d$w, trim = T, maxgap = nptperyear / 4, alpha = 0.02)

    temp <- list()
    for (j in 1:6){
        temp[[j]] <- tryCatch({
            I  <- ((j-1)*3*23+1):(j*3*23)
            input <- lapply(INPUT[1:3], `[`, I) %>% c(INPUT[5])
            vc    <- v_curve(input$y, w = input$w, llas = seq(-2, 2, by = 0.01), d = 2, show = F)

            listk(site = sitename, IGBP, lambda = vc$lambda, vc)
        }, error = function(e){
            message(sprintf("[%05d]: %s", i, e$message))
        })
    }
    res[[i]] <- temp
}

a <- do.call(c, res) %>% rm_empty()
x <- map_df(a, ~llply(.x[c("site", "IGBP", "lambda")], first, default = NA)) %>% data.table()

x$lambda %<>% log10()
x <- x[!is.na(lambda), ]
x[, grp := (0:(.N-1)), site]

fwrite(x, "st1e4_lambda.csv")
## statistics
dt[, grp := floor(1:.N/(23*3)), site]
dt <- dt[grp <= 5]

# info <- ddply(dt, .(site), function(d){
#     y <- d$y[1:(23*3)]
#     c(mean = mean(y, na.rm = T),
#       sd = sd(y, na.rm = T),
#       kurtosis = kurtosis(y, na.rm = T, type = 2),
#       skewness = skewness(y, na.rm = T, type = 2))
# })
info <- dt[, .(mean = mean(y, na.rm = T),
               sd = sd(y, na.rm = T),
               kurtosis = kurtosis(y, na.rm = T, type = 2),
               skewness = skewness(y, na.rm = T, type = 2)), .(site, grp)]
info[, cv := sd/mean]

## 01. construct connection
stat <- merge(data.table(x), info, by = c("site", "grp"))

slm  <- function(d){
    lm_model    <- lm(lambda~kurtosis+mean+sd+skewness+I(sd/mean), d, na.action = na.exclude) #
    slm_model   <- stepAIC(lm_model)
    d$fitted    <- predict(slm_model)
    print(summary(slm_model))

    d2 <- dt[, .(lambda2 = init_lambda(y)), .(site, grp)]

    ggplot(d, aes(lambda, fitted)) + #, color = IGBP
        geom_point(alpha = 0.3) + #, color = "blue"
        geom_abline(slope = 1, intercept = 0, color = "red", size = 1) +
        facet_wrap(~IGBP)
}
slm(stat)

dt[, .(lambda2 = init_lambda(y)), .(site, grp)]
## 02. visualization
xs <- melt(x, c("site", "IGBP", "lambda"))
ggplot(xs, aes(value, lambda)) + geom_point(alpha = 0.2) + geom_smooth() +
    facet_wrap(~variable, scales = "free")
# facet_grid(IGBP~variable, scales = "free")

ggplot(stat, aes(lambda, fitted, color = IGBP)) + geom_point()


