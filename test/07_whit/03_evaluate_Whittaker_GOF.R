source('test/stable/load_pkgs.R')
source("test/GEE/main_phenofit_test.R")
# source("R/plot_phenofit.R")

stations212 <- fread("F:/Github/MATLAB/PML/data/flux-212.txt")

infile     <- file_cam
dir_gdrive <- "D:/Document/GoogleDrive/"

df_cam  <- get_phenofit_result(file_cam)
df_flux <- get_phenofit_result(file_flux)

st_cam  <- fread(file_st_cam)
st_flux <- fread(file_st_flux)

## agreement index
i <- 1
prefix <- c("phenoflux", "phenocam")[i]
df <- if(i == 1) df_flux else df_cam
st <- if(i == 1) st_flux else st_cam

df      <- df[iters == "iter2"]
st$IGBPname %<>% factor(IGBPnames_006)

methods <- c('AG', 'BECK', 'ELMORE', 'ZHANG', 'whit_R', 'whit_gee')

# make sure different curve fitting methods have the same length fitting
formula <- if(i == 1) formula(site+date+t+y+GPP_NT+GPP_DT+SummaryQA~meth) else
    formula(site+date+t+y+gcc+vci+SummaryQA~meth)

# only period when all curve fitting methods have result is kept.
df_trim <- dcast(df, formula, value.var = "val", fun.aggregate = mean) %>% na.omit()
df_trim <- melt(df_trim, measure.vars = methods, variable.name = "meth")

# visualization
info_ai <- df_trim[SummaryQA == "good", .(ai = agr_index(y, value)), .(site, meth)] %>% merge(st)
if(i == 1){
    info_r  <- df_trim[, .(R = stat_fun(value, GPP_DT)), .(site, meth)] %>% merge(st)
}else{
    info_r  <- df_trim[, .(R = stat_fun(value, vci)), .(site, meth)] %>% merge(st)
}

source('test/stable/geom_boxplot_no_outlier.R')

# 1. show correlation
p1 <- ggplot(info_r, aes(IGBPname, R, colour = meth), position = "dodge") %>%
    boxplot() %>% `+`(labs(x = "IGBP", y = "Correlation (r)"))

save_pdf(sprintf("valid_%s_R.pdf", prefix), 12, 5, p1)

p2 <- ggplot(info_ai, aes(IGBPname, ai, colour = meth), position = "dodge") %>%
    boxplot() %>% `+`(labs(x = "IGBP", y = "Agreement Index (AI)"))
save_pdf(sprintf("valid_%s_AI.pdf", prefix), 12, 5, p2)

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
}
# writelist_ToXlsx(listk(c_min, c_max), "gee_info_count.xlsx")

sites <- unique(df_trim$site)
sitename <- sites[100]

################################################################################
colnames(df_trim)[4] <- "raw"
qc_levels <- c("good", "margin", "snow/ice", "cloud")
df_trim$SummaryQA %<>% factor(qc_levels)

get_range <- function(d, grp = "b1", alpha = 0){
    if (length(alpha) == 1) alpha %<>% rep(2)
    res <- d[, .(min = quantile(value, alpha[1], na.rm = T),
                 max = quantile(value, alpha[2], na.rm = T))]
    unlist(res)
}

plot_whit <- function(sitename, df, st){
    ## figure title and filename
    sp    <- st[site == sitename, ] # station point
    titlestr <- with(sp, sprintf('[%03d,%s] %s, lat = %5.2f, lon = %6.2f',
                                     ID, site, IGBPname, lat, lon))
    file_pdf <- sprintf('Figure/%s_[%03d]_%s.pdf', prefix_fig, sp$ID[1], sp$site[1])

    ##
    x <- df_trim[site == sitename, ]

    d <- melt(x, measure.vars = colnames(x)[c(4, 5:6, 8:13)], variable.name = "meth")
    d_valid <- x[, .(valid = (GPP_DT - coef[1])/coef[2]), .(site, date, t)]

    ## scale validation variable (e.g. GPP or VCI), to keep the same range as 
    # `whit_gee`
    lim_fit <- get_range(d[ grep("whit_gee", meth)], "fit", alpha = [0.01, 0.99]) 
    lim_valid <- get_range(d[ grep("GPP", meth)],  "valid", alpha = [0.01, 0.99])
    
    lim_raw <- get_range(d[ grep("raw", meth)], "raw") # ylim
    
    # r <- phenofit:::.normalize(y_fit, y_raw)
    # lim_valid_adj <- lm(lim_valid~r) %>% predict(data.frame(r = c(0, 1)))
    coef <- lm(lim_valid~lim_fit) %>% coef() # coefs used to rescale validation data

    ## ggplot, not only whit_gee, I also need to know all curve fitting methods 
    #  performance
    p1 <- ggplot(d[-grep("raw|GPP", meth)], aes(date, value)) +
        geom_point(aes(shape = SummaryQA, color = SummaryQA)) +
        geom_line(data = d_valid, aes(date, valid), size = 0.9, color = "blue") +
        geom_line(data = d[meth == "whit_gee", .(date, value, SummaryQA)],
                  size =1, color = "black") +
        labs(y = "EVI") +
        theme(legend.position = "none") +
        scale_y_continuous(lim = lim_raw,
                           sec.axis = sec_axis(~.*coef[2]+coef[1],
                                               name = expression("GPP ( gC "*mm^-1*d^-1*" )"))) +
        facet_wrap(~meth, ncol = 1) +
        scale_color_manual(values = c("good" = "grey60", "margin" = "#00BFC4",
                                      "snow/ice" = "#F8766D", "cloud" = "#C77CFF"), drop = F) +
        scale_shape_manual(values = c(19, 15, 4, 17), drop = F) +
        ggtitle(titlestr)

}



# p0 <- ggplot(d[meth == "raw"], aes(date, value, shape = SummaryQA, color = SummaryQA)) +
#     geom_point() +
#     theme(legend.position = "none") +
#     scale_color_manual(values = c("good" = "grey60", "margin" = "#00BFC4",
#                                   "snow/ice" = "#F8766D", "cloud" = "#C77CFF"), drop = F) +
#     scale_shape_manual(values = c(19, 15, 4, 17), drop = F) +
#     scale_y_continuous(lim = lim_raw)

# p3 <- ggplot_dual_axis(p1, p2) #%>% as.ggplot()
# p <- ggplot_dual_axis(p3, p0, add_yaxis_r = F)
