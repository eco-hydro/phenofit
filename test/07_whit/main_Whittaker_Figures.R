
fix_name <- function(x) {
    names(x) %<>% str_extract(".*(?=_)")
    # names(x)
    x
}

tidy_info <- function(file){
    x <- readRDS(file)
    a <- llply(x, fix_name) %>% purrr::transpose() %>%
        llply(function(l) {
            info <- map(l, "info") %>% do.call(rbind, .)
            rough <- map(l, "rough") %>% do.call(rbind, .)
            list(info = info, rough = rough)
        })
    names <- names(a)

    for (i in 1:length(a)){
        name <- names[i]
        d_info  <- a[[i]]$info
        d_rough <- a[[i]]$rough

        if (name != "phenofit"){
            d_info$meth <- name
            d_rough$meth <- name
        }
        d_info  %<>% reorder_name(c("site", "meth"))
        d_rough %<>% reorder_name(c("site", "meth"))

        a[[i]] <- merge(d_info, d_rough) %>%
            reorder_name(c("site", "meth", "type", "iters"))#list(info = d_info, rough = d_rough)
    }
    a %<>% do.call(rbind, .) #transpose() %>% map(~
    # a <- fix_name(a)
    # d <- .[NSE > 0, ] %>%
    #     melt(id.vars = c("site", "meth", "type"), variable.name = "index") %>%
    #     .[index %in% c("R2", "NSE", "RMSE")]

    # d$meth %<>% factor(methods)
    # a$Rg %<>% unlist()
    return(a)
}

boxplot2 <- function(p, width = 0.95, size = 0.7){
    # width  <- 0.95
    width2 <- width - 0.15
    dodge <- position_dodge(width = width)

    p +
        stat_summary(fun.data = box_qtl,
                     position = dodge, size = size,
                     geom = "errorbar", width = width2) +
        geom_boxplot2(coef = 0,
                      width = width2,
                      lwd = size - 0.2,
                      notch = F, outlier.shape = NA, position=dodge) +
        grid_x +
        geom_text(data = d_lab, aes(x = "ENF",
                                    y = Inf, color = NULL, label = label),
                  vjust = 1.5, hjust = 1.1, fontface = "bold", size =5, show.legend = F)
}