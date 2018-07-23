windowsFonts(Times = windowsFont("Times New Roman"),
             Arial = windowsFont("Arial"))
fontsize = 14

# save pdf just like `ggsave`
save_pdf <- function(file = "Rplot.pdf", width = 10, height = 5, p, open = F){
    if (missing(p)) p <- last_plot()


    if ("grob" %in% class(p)) {
        FUN <- grid::grid.draw
    } else{
        FUN <- print
    }
    print(FUN)
    # Cairo::CairoPDF, if only one figure cairo_pdf is the best
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
boxplot <- function(p, width = 0.95){
    # width  <- 0.95
    width2 <- width - 0.15
    dodge <- position_dodge(width = width)

    p + stat_summary(fun.data = box_qtl,
                     position = dodge,
                     geom = "errorbar", width = width2) +
        geom_boxplot2(coef = 0,
                  width = width2,
                 lwd = 0.3,
                 notch = F, outlier.shape = NA, position=dodge) +
    theme_light(base_size = fontsize, base_family = "Arial") +
    theme(legend.position = c(1-0.01, 0.01), legend.justification = c(1, 0),
          panel.grid.major = element_line(linetype = 2),
          panel.grid.minor = element_blank(),
          legend.title=element_blank(),
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

methods <- c('AG', 'BECK', 'ELMORE', 'ZHANG', 'whit_R', 'whit_gee')[-5]
#' @examples
#' over_perform(df, formula, prefix)
over_perform <- function(d, formula, prefix, ylim2, IGBP.all = F){
    # only period when all curve fitting methods have result is kept.
    df_trim <- dcast(d, formula, value.var = "value", fun.aggregate = mean) %>% na.omit()
    df_trim[, raw := y]
    df_trim <- melt(df_trim, measure.vars = c( "raw", methods), variable.name = "meth")

    new_levels <- c("raw ", "AG ", "Beck ", "Elmore ", "Zhang ", "WH")
    cols <- hue_pal()(5) %>% c("black", .)%>% set_names(new_levels)

    df_trim$meth %<>% mapvalues(c("raw", methods), new_levels)
    # visualization
    info_ai <- df_trim[SummaryQA == "good" & meth != "raw ",
                       .(ai = agree_index(y, value)), .(site, meth)] %>% merge(st)
    if(i == 1){
        info_r  <- df_trim[, .(R = stat_fun(value, GPP_DT)), .(site, meth)] %>% merge(st)
    }else{
        info_r  <- df_trim[, .(R = stat_fun(value, vci)), .(site, meth)] %>% merge(st)
    }

    # add a column 'all', independent of IGBP
    add_IGBPall <- . %>% {
        temp <- .; temp$IGBPname <- "all"
        rbind(., temp)
    }
    if (IGBP.all){
        info_ai %<>% add_IGBPall()
        info_r  %<>% add_IGBPall()
    }

    # 1. show correlation, , alpha = 0.6, fill = meth
    p1 <- { ggplot(info_r, aes(IGBPname, R, color = meth), position = "dodge") +
                scale_colour_manual(values = cols,
                                    guide = guide_legend(direction = "horizontal", nrow = 1, keywidth = 1)) +
                # scale_fill_manual(values = cols) +
                labs(x = "IGBP", y = "Correlation (r)")
          } %>% boxplot()
    p1 <- p1 +
        # geom_text(data = info_r[1, ], aes(fontface = "bold"), x = -Inf, y = -Inf,
        #           label = c("(a)"), hjust = -1, vjust = -2.6,
        #           show.legend = F, size = 4) +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              plot.margin = margin(3, 3, 2, 3))

    p2 <- {ggplot(info_ai, aes(IGBPname, ai, colour = meth), position = "dodge") +
        scale_colour_manual(values = cols) } %>%
        boxplot() %>% `+`(labs(x = "IGBP", y = "Agreement Index (AI)"))
    p2 <-p2 + theme(legend.position = 'none',
                    plot.margin = margin(3, 3, 2, 3))

    if (!missing(ylim2)) p2 <- p2 + coord_cartesian(ylim = ylim2)
    p <- gridExtra::arrangeGrob(p1, p2, nrow = 2, heights = c(0.9, 1), padding = unit(1, "line"))
    # grid.draw(p)
    save_pdf(sprintf("valid_%s.pdf", prefix), 9, 8, p, open = T)
    # save_pdf(sprintf("valid_%s_AI.pdf", prefix), 12, 5, p2)
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

# re-select colors
cols <- colors()[c(74, 134, 50)]
cols <- c("blue", "red", cols[3]) #"#B2B302"
color_valid <- cols[3] #"green4"

lgd_vci <- phenofit:::make_legend(linename = c("iter1", "iter2", "VCI"),
                   linecolor = cols)
lgd_gpp <- phenofit:::make_legend(linename = c("iter1", "iter2", "GPP"),
                   linecolor = cols)

#' plot_methods
#'
#' plot curve fitting series and validation data (i.e.GPP or VCI ) to check the
#' smoothing performance
#'
#' @param df_trim A data.table of two iters curve fitting result.
#' For MODIS product, the current year last value maybe same as
#' the coming year first value. Need to remove the duplicated data. Besides,
#' We don't constrain the equal length of different curve fitting series as
#' performance index part.
#' @param st A dataframe of station information, ID, site, IGBPname, lat, lon.
#' @param methods one of 'AG', 'BECK', 'ELMORE', 'ZHANG', 'whit_R' and 'whit_gee'.
#'
#' @examples
#' \dontrun{
#' plot_whit(sitename, df_trim, st, prefix_fig = "whit")
#' }
plot_methods <- function(sitename, df_trim, st, prefix_fig = "whit", methods, show.legend = T){
    ## figure title and filename
    sp    <- st[site == sitename, ] # station point
    # titlestr <- with(sp, sprintf('[%03d,%s] %s, lat = %5.2f, lon = %6.2f',
    #                                  ID, site, IGBPname, lat, lon))
    titlestr <- sp$titlestr
    if ( length(titlestr) == 0 || is.na(titlestr) ){
        titlestr <- with(sp, sprintf('[%03d,%s] %s, lat = %5.2f, lon = %6.2f',
                                     ID, site, IGBPname, lat, lon))
    }
    file_pdf <- sprintf('Figure/%s_[%03d]_%s.pdf', prefix_fig, sp$ID[1], sp$site[1])

    ##
    x <- df_trim[site == sitename , ]
    if (all(is.na(x$whit_gee))) return()
    # d <- melt(x, measure.vars = methods, variable.name = "meth")
    d <- melt(x,
              measure.vars = c(contain(x, "^y$|GPP_NT|GPP_DT|vci|gcc"), methods),
              variable.name = "meth")
    pdat  <- d[meth %in% methods]

    ## scale validation variable (e.g. GPP or VCI)
    # The range of validation variable and \code{whit_gee} should be equal.
    d_valid_scale <- d[meth %in% c("whit_gee", "GPP_DT", "vci") & iters == "iter2"] %>%
        dcast(site+date~meth, value.var = "value") %>% na.omit() %>%
        melt(id.vars = c("site", "date"), variable.name = "meth")
    lim_fit   <- get_range(d_valid_scale[ grep("whit_gee", meth)], alpha = c(0.01, 0.99))
    lim_valid <- get_range(d_valid_scale[ grep("GPP|vci", meth)], alpha = c(0.01, 0.99))
    lim_raw   <- get_range(d[ grep("y", meth)]) # ylim

    # r <- phenofit:::.normalize(y_fit, y_raw)
    # lim_valid_adj <- lm(lim_valid~r) %>% predict(data.frame(r = c(0, 1)))
    coef <- lm(lim_valid~lim_fit) %>% coef() # coefs used to rescale validation data

    ## only keep iter2
    x <- df_trim[site == sitename & iters == "iter2", ]
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
    # y|GPP|vci|gcc|

    IsSingle <- length(methods) == 1
    if (IsSingle){
        d_lab <- data.frame(meth = methods, lab = titlestr)
    }else{
        d_lab <- data.frame(meth = methods,
                         lab = sprintf("(%s) %-8s", letters[1:length(methods)], methods))
    }
    lwd <- 0.65
    # color_valid <- "green4" # set as global variable
    p1 <- ggplot(pdat, aes(date, value)) +
        geom_vline(xintercept = ymd(0101 + (2001:2017)*1e4),
                   color = "white", linetype = 3, size = 0.4) +
        geom_point(data = d_raw, aes(date, y, shape = SummaryQA, color = SummaryQA), size = 1.4) +
        geom_line(data = pdat[iters == "iter1"], color = "blue", size = lwd) +
        geom_line(data = pdat[iters == "iter2"], color = "red", size = lwd) +
        geom_line(data = d_valid, aes(date, valid), size = lwd, color = color_valid) +
        labs(y = "EVI") +
        theme(legend.position = "none",
              # axis.text.y.left = element_text(color = cols[2]), #iter2 color
              # axis.title.y.left = element_text(color = cols[2]),
              axis.text.y.right = element_text(color = color_valid),
              axis.title.y.right = element_text(color = color_valid),
              plot.margin = margin(t = 4, r = 2, b = 0, l = 2, unit = "pt"),
              legend.margin = margin(),
              strip.text = element_blank(),
              panel.grid.minor = element_blank(),
              # panel.grid.major.x = element_blank(),
              # panel.grid.major.y = element_line(size = 0.2),
              panel.grid = element_line(size = 0.4) #linetype = 2
              # axis.ticks.y.right = element_text(color = "blue"),
              ) +
        scale_y_continuous(lim = lim_raw,
                           sec.axis = sec_axis(~.*coef[2]+coef[1], name = ylab_r)) +
        scale_x_date(breaks = ymd(0101 + seq(2001, 2017, 4)*1e4)) +
        facet_wrap(~meth, ncol = 1) +
        scale_color_manual(values = qc_colors, drop = F) +
        scale_shape_manual(values = qc_shapes, drop = F) +
        geom_text(data = d_lab, aes(label = lab), fontface = "bold",
            x = -Inf, y =Inf, vjust = 1.5, hjust = -0.08)

    if (!IsSingle) p1 <- p1 + ggtitle(titlestr)

    if (show.legend)
        p1 <- gridExtra::arrangeGrob(p1, lgd, nrow = 2,
            heights = c(20, 1), padding = unit(0.5, "line")) #return,

    if (!IsSingle) save_pdf(file_pdf, 11, 7, p = p1)
    p1
}




