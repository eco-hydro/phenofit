fontsize = 12

# save pdf just like `ggsave`
save_pdf <- function(file = "Rplot.pdf", width = 10, height = 5, p, open = F){
    if (missing(p)) p <- last_plot()

    FUN <- print
    if ("grob" %in% class(p)) FUN <- grid::grid.draw

    Cairo::CairoPDF(file, width = width, height = height)
    FUN(p)
    dev.off()
    if (open) file.show(file)
}

# get curve fitting results from phenofit object
getFittings2 <- function(fit){
    df_fit <- getFittings(fit) %>% data.table()

    whit <- melt(fit$seasons$whit, measure.vars = c("iter1", "iter2"),
                 variable.name = "iters")
    whit$meth <- "whit_R"
    df_fit <- rbind(df_fit, whit)
    df_fit <- unique(df_fit) # remove duplicated value

    return(df_fit)
}

stat_fun <- function(Y_obs, Y_sim){
    R      <- NA_real_
    pvalue <- NA_real_
    tryCatch({
        cor.obj <- cor.test(Y_obs, Y_sim, use = "complete.obs")
        R       <- cor.obj$estimate[[1]] # statistic
        pvalue  <- cor.obj$p.value
    }, error = function(e){
        message(e$message)
    })
    # c(R = R, pvalue = pvalue)[1]
    # pvalue
    R
}

table_count <- function(x, levels){
    t <- table(x)
    I <- match(names(t), levels)

    res <- rep(NA, length(levels))
    res[I] <- t
    set_names(res, levels)
}

# agreement index
agree_index <- function(Y_obs, Y_sim){
    I <- which(!(is.na(Y_sim) | is.na(Y_obs))) # | is.na(w)))
    # n_obs <- length(Y_obs)
    n_sim <- length(I)

    Y_sim <- Y_sim[I]
    Y_obs <- Y_obs[I]

    u <- mean(Y_sim)

    d = 100 - sum((Y_sim - Y_obs)^2) / sum( (abs(Y_sim - u) + abs(Y_obs - u))^2 )*100
    d # agreement index
}

# boxplot for over all correlation and agreement index
boxplot <- function(p, width = 0.8){
    p + stat_summary(fun.data = box_qtl,
                 position = position_dodge(width = width),
                 geom = "errorbar", width = width) +
    geom_boxplot(coef = 0, width = width,
                 lwd = 0.3,
                 notch = F, outlier.shape = NA) +
    theme_light(base_size = fontsize, base_family = "Arial") +
    theme(legend.position = c(1-0.01, 0.01), legend.justification = c(1, 0),
          panel.grid.major = element_line(linetype = 2),
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"))
}

################################################################################
# merge phenofit INPUT, OUTPUT, GEE whittaker result and validation data
get_phenofit_result <- function(infile){
    prefix  <- str_extract(infile, "\\w*(?=_MOD)")

    ## 1. INPUT
    df_in   <- fread(infile, strip.white = T)                               # phenofit INPUT
    df_in$t    %<>% ymd()
    df_in$date %<>% ymd()

    ## 2.1 phenofit curve fitting
    lst     <- get_slurm_out(paste0("Y:/github/phenofit_cluster/result/", prefix, "/"))
    df_temp <- llply(lst, getFittings2, .progress = "text")
    df_out  <- melt_list(df_temp, "site") %>% data.table() # phenofit OUTPUT

    ## 2.2 whit curve fitting
    df_whit  <- readwhitMAT(dir_gdrive, prefix)            # GEE whit smoothing
    measure.vars <- colnames(df_whit) %>% .[grep("iter", .)]
    # unify curve fitting output format as phenofit
    df_whit <- merge(df_in[, .(site, t, date)], df_whit[, -1], by = c("site", "date")) #rm raw data in gee
    df_whit <- melt(df_whit, id.vars = c("site", "t"), measure.vars = measure.vars,
         variable.name = "iters")
    df_whit$meth <- "whit_gee"

    df_fits <- rbind(df_out[, .(site, t, iters, value, meth)],
                     df_whit[, .(site, t, iters, value, meth)])

    ## 3. validation data
    valid_file <- sprintf("%svalid_%s_16day.csv", dir_data, prefix)
    df_valid <- fread(valid_file)
    df_valid$date %<>% ymd()

    # remove NA values
    varnames <- colnames(df_valid)[4:5]
    eval(parse(text = sprintf("df_valid = df_valid[!(is.na(%s) & is.na(%s))]",
        varnames[1], varnames[2])))

    res <- merge(df_in, df_fits, c("site", "t"), all.x = T) %>%
        .[date < ymd(20180101) & date >= ymd(20000218)] %>%
        merge(df_valid[, c(1, 6, 4:5)], c("site", "date"), all.x = T)

    ## 4. station info
    # st_file <- sprintf("%sst_%s.csv", dir_data, prefix)
    # st <- fread(st_file)   # station info
    res
}

methods <- c('AG', 'BECK', 'ELMORE', 'ZHANG', 'whit_R', 'whit_gee')
#' @examples
#' over_perform(df, formula, prefix)
over_perform <- function(df, formula, prefix){
    # only period when all curve fitting methods have result is kept.
    df_trim <- dcast(df, formula, value.var = "value", fun.aggregate = mean) %>% na.omit()
    df_trim <- melt(df_trim, measure.vars = methods, variable.name = "meth")

    # visualization
    info_ai <- df_trim[SummaryQA == "good", .(ai = agree_index(y, value)), .(site, meth)] %>% merge(st)
    if(i == 1){
        info_r  <- df_trim[, .(R = stat_fun(value, GPP_DT)), .(site, meth)] %>% merge(st)
    }else{
        info_r  <- df_trim[, .(R = stat_fun(value, vci)), .(site, meth)] %>% merge(st)
    }

    # 1. show correlation
    p1 <- ggplot(info_r, aes(IGBPname, R, colour = meth), position = "dodge") %>%
        boxplot() %>% `+`(labs(x = "IGBP", y = "Correlation (r)"))

    save_pdf(sprintf("valid_%s_R.pdf", prefix), 12, 5, p1)

    p2 <- ggplot(info_ai, aes(IGBPname, ai, colour = meth), position = "dodge") %>%
        boxplot() %>% `+`(labs(x = "IGBP", y = "Agreement Index (AI)"))
    save_pdf(sprintf("valid_%s_AI.pdf", prefix), 12, 5, p2)
}

# geom_point(aes(fill = meth), pch = 21,
#            position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9))
# ggplot(info_df, aes(meth, R), position = "dodge") +
#     stat_summary(fun.data = box_qtl,
#                  position = position_dodge(width = 0.9),
#                  geom = "errorbar", width = 0.9) +
#     geom_boxplot2(coef = 0, width = 0.9, notch = F, outlier.shape = NA)

# get dominant method occurred times
stat_dominant <- function(){
    info <- dcast(info_r, site+lat+lon+IGBPname~meth, value.var = "R")
    cols_del <- c(1:4, 9) # 9:whit_R
    methods <- colnames(info)[-cols_del]

    mat <- info[, -cols_del, with = F] %>% as.matrix()
    I <- rowSums(is.na(mat)) == 0

    best <- mat[I, ] %>% {
        data.table(min = methods[apply(., 1, which.min)],
                   max = methods[apply(., 1, which.max)])
    } %>% cbind(info[I, 1:4], .)
    best
    info %<>% cbind(best)

    a_min <- ddply(info, .(IGBP), function(d) table_count(d$min, methods))
    a_max <- ddply(info, .(IGBP), function(d) table_count(d$max, methods))

    t_min <- table(best$min)
    t_max <- table(best$max)
    listk(info, a_min, a_max, t_min, t_max)
    # writelist_ToXlsx(listk(c_min, c_max), "gee_info_count.xlsx")
}

get_range <- function(d, alpha = c(0, 1)){
    if (length(alpha) == 1) alpha %<>% rep(2)
    res <- d[, .(min = quantile(value, alpha[1], na.rm = T),
                 max = quantile(value, alpha[2], na.rm = T))]
    unlist(res)
}

lgd_vci <- phenofit:::make_legend(linename = c("Curve fitting", "VCI"),
                   linecolor = c("black", "blue"))
lgd_gpp <- phenofit:::make_legend(linename = c("Curve fitting", "GPP"),
                   linecolor = c("black", "blue"))

#' @examples
#' plot_whit(sitename, df_trim, st, prefix_fig = "whit")
plot_whit <- function(sitename, df_trim, st, prefix_fig = "whit"){
    ## figure title and filename
    sp    <- st[site == sitename, ] # station point
    titlestr <- with(sp, sprintf('[%03d,%s] %s, lat = %5.2f, lon = %6.2f',
                                     ID, site, IGBPname, lat, lon))
    file_pdf <- sprintf('Figure/%s_[%03d]_%s.pdf', prefix_fig, sp$ID[1], sp$site[1])

    ##
    x <- df_trim[site == sitename, ]
    if (all(is.na(x$whit_gee))) return()
    # d <- melt(x, measure.vars = methods, variable.name = "meth")
    d <- melt(x, measure.vars = c(colnames(x)[c(4, 5:6)], methods), variable.name = "meth")
    ## scale validation variable (e.g. GPP or VCI), to keep the same range as
    # `whit_gee`
    lim_fit   <- get_range(d[ grep("whit_gee", meth)], alpha = c(0.01, 0.99))
    lim_valid <- get_range(d[ grep("GPP|vci", meth)], alpha = c(0.01, 0.99))

    lim_raw   <- get_range(d[ grep("y", meth)]) # ylim

    # r <- phenofit:::.normalize(y_fit, y_raw)
    # lim_valid_adj <- lm(lim_valid~r) %>% predict(data.frame(r = c(0, 1)))
    coef <- lm(lim_valid~lim_fit) %>% coef() # coefs used to rescale validation data

    if ("vci" %in% colnames(x)){
        d_valid <- x[, .(valid = (vci - coef[1])/coef[2]), .(site, date, t)]
        ylab_r  <- "VCI"
        lgd     <- lgd_vci
    }else{
        d_valid <- x[, .(valid = (GPP_DT - coef[1])/coef[2]), .(site, date, t)]
        ylab_r  <- expression("GPP ( gC "*mm^-1*d^-1*" )")
        lgd     <- lgd_gpp
    }

    d_raw   <- x[, .(date, y, SummaryQA)]

    ## ggplot, not only whit_gee, I also need to know all curve fitting methods
    #  performance
    p1 <- ggplot(d[-grep("y|GPP|vci|gcc", meth)], aes(date, value)) +
        geom_line(color = "black", size = 0.9) +
        geom_line(data = d_valid, aes(date, valid), size = 0.9, color = "blue") +
        geom_point(data = d_raw, aes(date, y, shape = SummaryQA, color = SummaryQA), size = 1.2) +
        labs(y = "EVI") +
        theme(legend.position = "none",
              axis.text.y.right = element_text(color = "blue"),
              axis.title.y.right = element_text(color = "blue"),
              plot.margin = margin(t = 4, r = 2, b = 0, l = 2, unit = "pt"),
              legend.margin = margin(),
              # axis.ticks.y.right = element_text(color = "blue"),
              ) +
        scale_y_continuous(lim = lim_raw,
                           sec.axis = sec_axis(~.*coef[2]+coef[1], name = ylab_r)) +
        facet_wrap(~meth, ncol = 1) +
        scale_color_manual(values = qc_colors, drop = F) +
        scale_shape_manual(values = qc_shapes, drop = F) +
        ggtitle(titlestr)

    df_lab <- data.frame(meth = methods,
                         lab = sprintf("(%s) %-8s", letters[1:length(methods)], methods))

    p1 <- p1 + geom_text(data = df_lab, x = -Inf, y =Inf, vjust = 1.5, hjust = -0.08,
                   aes(label = lab), fontface = "bold") +
        theme(strip.text = element_blank())

    #
    p1 <- gridExtra::arrangeGrob(p1, lgd, nrow = 2, heights = c(20, 1), padding = unit(0.5, "line")) #return,

    save_pdf(file_pdf, 11, 7, p = p1)
}
