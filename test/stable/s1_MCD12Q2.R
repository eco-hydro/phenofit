source('F:/Github/phenology/China_Phenology/R/phenofit/test/stable/geom_boxplot_no_outlier.R')

site_rm2 <- c("AU-Cum", "AU-Cpr", "AU-Gin", "AU-RDF", "AU-Stp", "AU-Tum","AU-Wac",
    "IT-Noe","BR-Sa3", "CH-Oe2","DE-Kli", "FR-Pue","GF-Guy", "IT-BCi", "US-Blo", "US-SRG","US-Whs")

## modis phenology product
# pheno[, flag := sprintf("%d_1", year(date))]
# df <- merge(pheno, df_OBS)
#' @param df A data.table, with columns as the below
# # A tibble: 3,100 x 28
#    site   flag   date       Decrease Increase Maximum Minimum origin     TRS1.sos TRS1.eos TRS2.sos TRS2.eos TRS5.sos
#    <chr>  <fct>  <date>        <int>    <int>   <int>   <int> <date>        <dbl>    <dbl>    <dbl>    <dbl>    <dbl>
#  1 AT-Neu 2003_1 2003-01-01      203       69     145     314 2003-01-01     69.0      295     81.0      277      103
#  2 AT-Neu 2003_1 2003-01-01      203       69     145     314 2003-01-01     70.0      290     83.0      274      103
# # ... with 3,090 more rows, and 15 more variables: TRS5.eos <dbl>, TRS6.sos <dbl>, TRS6.eos <dbl>, DER.sos <dbl>,
# #   DER.pop <dbl>, DER.eos <dbl>, GU.UD <dbl>, GU.SD <dbl>, GU.DD <dbl>, GU.RD <dbl>, ZHANG.Greenup <dbl>,
# #   ZHANG.Maturity <dbl>, ZHANG.Senescence <dbl>, ZHANG.Dormancy <dbl>, meth <chr>
gof_prod_MCD12Q2 <- function(df, type = c('SOS', 'EOS'), trim = TRUE, IsPlot = TRUE, ...){
    if(type[1] == 'SOS'){
        d <- select(df, matches("sos|UD|SD|Greenup|Maturity")) - df$Increase
    }else{
        d <- select(df, matches("eos|DD|RD|Senescence|Dormancy")) - df$Minimum
    }
    d %<>% cbind(df[, .(flag, origin, meth, site)], .) #cbind header columns

    diff_boxplot(d, trim = trim, title = type[1], IsPlot, ...)
}

# d.sos <- cbind(header[, .(meth, site)], diff) %>%
#     .[, lapply(.SD, mean, na.rm = T), .(meth, site)] %>%
#     .[order(meth, site)]

gof_prod_phenofit <- function(df_sim, df_obs,
    type = c('SOS', 'EOS'), real = 'TRS1',
    # aggregateCurveMethod = FALSE,
    trim = FALSE, IsPlot = TRUE, ...)
{
    # 1. match df_obs, and df_sim
    I <- match(with(df_obs, paste(meth, site, flag)), with(df_sim, paste(meth, site, flag)))
    df_sim <- df_sim[I]

    if (type[1] == 'SOS'){
        pattern_sim <- "sos|UD|SD|Greenup|Maturity"
        pattern_sim %<>% paste0("|pop") #if (real != "TRS1") 

        pattern_obs <- ifelse(real == 'TRS1', 'TRS1.sos', pattern_sim)
    }else{
        pattern_sim <- "eos|DD|RD|Senescence|Dormancy"
        pattern_obs <- ifelse(real == 'TRS1', 'TRS1.eos', pattern_sim)
    }

    d_obs <- select(df_obs, matches(pattern_obs))
    if (real == 'TRS1') d_obs <- d_obs[[1]]

    d <- select(df_sim, matches(pattern_sim)) - d_obs
    d %<>% cbind(df_sim[, .(flag, origin, meth, site)], .) #cbind header columns

    title <- paste(type[1], real, sep = '_')
    diff_boxplot(d, trim, title, IsPlot, ...)
}


diff_boxplot <- function(d, trim = trim, title, IsPlot = TRUE, ...){
    d <-  d[!is.na(site)] #remove na values in d

    df  <- gather(d, index, val, -flag, -origin, -meth, -site)
    gof <- ddply(df, .(site, index, meth), function(x) GOF2(x$val)) %>% as.data.table()
    gof.all <- ddply(df, .(index, meth), function(x) GOF2(x$val)) %>% as.data.table()

    p <- NULL
    if (IsPlot){
        ## data for plot, if trim is true, points RMSE large than 100 will be removed.
        if (trim){
            pdat <- gather(subset(gof[, 1:6], RMSE < 100), GOFIndex, val, -site, -index, -meth)
        }else{
            pdat <- gather(gof[, 1:6], GOFIndex, val, -site, -index, -meth)
        }

        p <- ggplot(pdat, aes(index, val, color = index)) +
            stat_summary(fun.data = box_qtl, geom = "errorbar", width = 0.5) +
            geom_boxplot2(notch = TRUE, outlier.shape = NA, coef = 0, width = 0.8) +
            facet_grid(GOFIndex~meth, scales = "free_y") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
            ggtitle(title) +
            labs(x = 'Phenophase', y = NULL, color = 'Phenophase') +
            guides(color=FALSE)
            # geom_boxplot() +
            # geom_jitter(width = 0.5) +
        # print(p)
    }
    return(list(d = d, gof = gof, gof.all = gof.all, plot = p))
}

