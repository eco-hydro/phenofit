# source('test/stable/geom_boxplot_no_outlier.R')
# stat_summary(fun.data = box_qtl, geom = "errorbar", width = 0.5) +
# geom_boxplot2(notch = TRUE, outlier.shape = NA, coef = 0, width = 0.8)
        
box_qtl <- function(x){
    x <- stats::na.omit(x)
    quantile(x, c(0.1, 0.9)) %>% set_names(c("ymin", "ymax"))
}
