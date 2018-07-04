source('test/stable/load_pkgs.R')
# source("R/plot_phenofit.R")

# main functions
stations212 <- fread("F:/Github/MATLAB/PML/data/flux-212.txt")
# get curve fitting results
getFit <- function(fit){
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

################################################################################

## 00. raw MOD13A1
df_raw     <- fread("data_test/phenofit166_MOD13A1_INPUT.csv", strip.white = F)
df_raw[, `:=`( t = ymd(t), SummaryQA = factor(SummaryQA, qc_levels))]

## 01. load MOD13A1 curve fitting
lst <- get_slurm_out("Y:/github/phenofit_cluster/result/flux166_MOD13A1/")
df_lst <- llply(lst, getFit, .progress = "text")

df_R   <- melt_list(df_lst, "site")
df_R   <- merge(df_raw[, .(site, t, date)], df_R) %>% # SummaryQA, IGBPname
    {.[date < ymd(20180101)]}
x <- merge(df_raw[, .(site, t, date, SummaryQA, lat, lon, IGBPname)], df_R, c("site", "t"))

# df_R[iters == "iter1", .N, .(site, meth)]
# df_fit[iters == "iter1", .N, .(meth)]
# df_R[, .(begin = first(date), end = last(date), n = .N), .(site, meth)]

## 02. gee_whit
indir <- "D:/Document/GoogleDrive/phenofit/data/gee_phenofit/v2/"
files <- dir(indir, "*.geojson", full.names = T)
df_gee <- read_whits.gee(files)
df_gee$date %<>% format()

df_gee <- merge(df_raw[, .(site, t, date)], df_gee[, -1], by = c("site", "date")) #rm raw data in gee
df_gee <- melt(df_gee,  measure.vars = colnames(df_gee)[4:ncol(df_gee)],
     variable.name = "iters", value.name = "val")
df_gee$meth <- "whit_gee"

df_fits <- rbind(df_R[, .(site, t, date, iters, val, meth)],
                 df_gee[, .(site, t, date, iters, val, meth)])
## 03. GPP data
df_gpp = fread("data_test/flux166_16d_ForWhittaker.csv")
df_gpp[, date := format(parse_date_time(paste0(year, (d16-1)*16+1), "%Y%j"))]
df_gpp <- df_gpp[!(is.na(GPP_NT) & is.na(GPP_DT)), ] # remove gpp NA value

res <- merge(df_fits, df_gpp[, .(site, date, GPP_NT, GPP_DT)], by = c("site", "date")) %>%
    merge(df_raw, ., by = c("site", "date", "t"))

# fwrite(res, "data_test/gee_phenofit_multi.csv")

## agreement index
agr_index <- function(Y_obs, Y_sim){
    I <- which(!(is.na(Y_sim) | is.na(Y_obs))) # | is.na(w)))
    # n_obs <- length(Y_obs)
    n_sim <- length(I)

    Y_sim <- Y_sim[I]
    Y_obs <- Y_obs[I]

    u <- mean(Y_sim)

    d = 100 - sum((Y_sim - Y_obs)^2) / sum( (abs(Y_sim - u) + abs(Y_obs - u))^2 )
    d # agreement index
}

a <- res[SummaryQA == " good", .(d = agr_index(y,val)), .(site, meth, IGBPname)]
ggplot(a, aes(IGBPname, d, color = meth)) + geom_boxplot()

info_df <- res[iters == "iter2", .(R = stat_fun(val, GPP_DT)), .(site, meth)] %>%
    merge(stations212[, .(site, lat, lon = long, IGBP)])

info <- dcast(info_df, site+lat+lon+IGBP~meth, value.var = "R")

cols_del <- c(1:4, 9) # 9:whit_R
methods <- colnames(info)[-cols_del]

best <- info[, -c(1:4, 9)] %>% {
    data.table(min = methods[apply(., 1, which.min)],
               max = methods[apply(., 1, which.max)])
}

info %<>% cbind(best)
table(info$max); table(info$min)

c_min <- ddply(info, .(IGBP), function(d) table_count(d$min, methods))
c_max <- ddply(info, .(IGBP), function(d) table_count(d$max, methods))
writelist_ToXlsx(listk(c_min, c_max), "gee_info_count.xlsx")

# ggplot(info_df, aes(IGBP, R, color = meth)) + geom_boxplot()
