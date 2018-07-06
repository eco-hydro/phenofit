# save pdf just like `ggsave`
save_pdf <- function(file = "Rplot.pdf", width = 10, height = 5, p){
    if (missing(p)) p <- last_plot()
    CairoPDF(file, width = 12, height = 5)
    print(p)
    dev.off()
    file.show(file)
}

# get curve fitting results from phenofit object
getFittings2 <- function(fit){
    df_fit <- getFittings(fit) %>% data.table()

    whit <- melt(fit$seasons$whit, measure.vars = c("iter1", "iter2"),
                 variable.name = "iters", value.name = "val")
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
agr_index <- function(Y_obs, Y_sim){
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
         variable.name = "iters", value.name = "val")
    df_whit$meth <- "whit_gee"

    df_fits <- rbind(df_out[, .(site, t, iters, val, meth)],
                     df_whit[, .(site, t, iters, val, meth)])

    ## 3. validation data
    valid_file <- sprintf("%svalid_%s_16day.csv", dir_data, prefix)
    df_valid <- fread(valid_file)
    df_valid$date %<>% ymd()

    # remove NA values
    varnames <- colnames(df_valid)[4:5]
    eval(parse(text = sprintf("df_valid = df_valid[!(is.na(%s) & is.na(%s))]",
        varnames[1], varnames[2])))

    res <- merge(df_in, df_fits, c("site", "t")) %>%
        .[date < ymd(20180101) & date >= ymd(20000218)] %>%
        merge(df_valid[, c(1, 6, 4:5)], c("site", "date"))

    ## 4. station info
    # st_file <- sprintf("%sst_%s.csv", dir_data, prefix)
    # st <- fread(st_file)   # station info
    res
}