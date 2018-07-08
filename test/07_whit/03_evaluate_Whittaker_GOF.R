source('test/stable/load_pkgs.R')
source("test/07_whit/main_phenofit_test.R")
source('test/stable/ggplot/geom_boxplot_no_outlier.R')
# source("R/plot_phenofit.R")

# stations212 <- fread("F:/Github/MATLAB/PML/data/flux-212.txt")

infile     <- file_cam
dir_gdrive <- "D:/Document/GoogleDrive/"

df_cam  <- get_phenofit_result(file_cam)
df_flux <- get_phenofit_result(file_flux)

st_cam  <- fread(file_st_cam)
st_flux <- fread(file_st_flux)

methods <- c('AG', 'BECK', 'ELMORE', 'ZHANG', 'whit_R', 'whit_gee')[-5]

##
i <- 2
prefix <- c("phenoflux", "phenocam")[i]
df <- if(i == 1) df_flux else df_cam
st <- if(i == 1) st_flux else st_cam

df      <- df[iters == "iter2"]
st$IGBPname %<>% factor(IGBPnames_006)

# make sure different curve fitting methods have the same length fitting
formula <- if(i == 1) formula(site+date+t+y+GPP_NT+GPP_DT+SummaryQA~meth) else
    formula(site+date+t+y+gcc+vci+SummaryQA~meth)

over_perform(df, formula, prefix)

# site figure data input
qc_levels <- c("good", "margin", "snow/ice", "cloud")
qc_colors <- c("grey60", "#00BFC4", "#F8766D", "#C77CFF") %>% set_names(qc_levels)
qc_shapes <- c(19, 15, 4, 17) %>% set_names(qc_levels)

df_trim <- dcast(df, formula, value.var = "value", fun.aggregate = mean) # %>% na.omit()
df_trim$SummaryQA %<>% factor(qc_levels)
# df_trim <- melt(df_trim, measure.vars = methods, variable.name = "meth")

sites <- unique(df$site)
sitename <- sites[100]

vars <- c("get_range", "save_pdf", "lgd",
          "qc_levels", "qc_colors", "qc_shapes", "methods")
cl <- cluster_Init(pkgs = c("data.table", "ggplot2", "magrittr"),
                   vars = vars)
clusterExport(cl, vars)
res <- parLapplyLB(cl, sites, plot_whit,
                   df_trim, st, prefix_fig = paste0("whit_", prefix))
for (i in seq_along(sites)){
    runningId(i)
    sitename <- sites[i]
    plot_whit(sitename, df_trim, st, prefix_fig = paste0("whit_", prefix))
}

# merge_pdf('../whit_phenoflux166.pdf', indir = "./")
# merge_pdf('whit_phenoflux166.pdf', indir = "Figure/")
# merge_pdf('whit_phenocam133.pdf', indir = "Figure/", del = T)
################################################################################
# colnames(df)[4] <- "raw"


get_range <- function(d, alpha = c(0, 1)){
    if (length(alpha) == 1) alpha %<>% rep(2)
    res <- d[, .(min = quantile(value, alpha[1], na.rm = T),
                 max = quantile(value, alpha[2], na.rm = T))]
    unlist(res)
}

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
        ylab_r <- "VCI"
    }else{
        d_valid <- x[, .(valid = (GPP_DT - coef[1])/coef[2]), .(site, date, t)]
        ylab_r <- expression("GPP ( gC "*mm^-1*d^-1*" )")
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

lgd <- phenofit:::make_legend(linename = c("Curve fitting", "VCI"),
                   linecolor = c("black", "blue"))
grid.newpage(); grid.draw(lgd)
# arrangeGrob(p, lgd, nrow = 2, heights = c(15, 1), padding = unit(1, "line")) #return,

# p0 <- ggplot(d[meth == "raw"], aes(date, value, shape = SummaryQA, color = SummaryQA)) +
#     geom_point() +
#     theme(legend.position = "none") +
#     scale_color_manual(values = c("good" = "grey60", "margin" = "#00BFC4",
#                                   "snow/ice" = "#F8766D", "cloud" = "#C77CFF"), drop = F) +
#     scale_shape_manual(values = c(19, 15, 4, 17), drop = F) +
#     scale_y_continuous(lim = lim_raw)

# p3 <- ggplot_dual_axis(p1, p2) #%>% as.ggplot()
# p <- ggplot_dual_axis(p3, p0, add_yaxis_r = F)
