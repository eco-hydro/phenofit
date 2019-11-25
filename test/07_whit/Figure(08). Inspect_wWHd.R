# Update 20181024
# ---------------
# This function is used to test the rough curve fitting methods, i.e. (smooth_wSG,
# smooth_wHANTS, wWHd, wWH2) at 16000 points, sampled at global scale.
source("test/load_pkgs.R")
source("R/season_mov.R")
source("R/curvefits.R")
source('R/smooth_wWHIT_lambda.R')
source('test/07_whit/main_gee_Whittaker.R')

## MAKE A legend
make_legend <- function(linename = c("iter1", "iter2", "whit"),
        linecolor = c("blue", "red", "black")){
    qc_levels <- c("good", "margin", "snow", "cloud")
    qc_colors <- c("grey60", "#00BFC4", "#F8766D", "#C77CFF") %>% set_names(qc_levels)
    qc_shapes <- c(19, 15, 4, 17) %>% set_names(qc_levels)
    npoints   <- length(qc_levels)

    labels <- c(qc_levels,linename)
    colors <- c(qc_colors, linecolor)

    # labels <- c(" good", " margin", " snow/ice", " cloud", linename)
    # colors <- c("grey60", "#00BFC4", "#F8766D", "#C77CFF", linecolor)
    nline <- length(linename)
    pch <- c(qc_shapes, rep(NA, nline))

    lty <- rep(0, npoints);  lty[3] <- 1
    lty <- c(lty, rep(1, nline))
    lwd <- c(rep(1, npoints), rep(3, nline))

    I   <- 1:length(colors)
    lgd <- grid::legendGrob(labels[I], pch = pch[I], nrow = 1,
                       # do.lines = T,
                       gp=grid::gpar(lty = lty[I], lwd = lwd[I],
                                     cex = 0.65, fontsize = 20,
                               col = colors[I], fill = colors[I]))
    lgd$children[[5]]$children[[1]]$children %<>% .[2] # fix cross point type
    return(lgd)
}

lgd <- make_legend(linename = c("raw", "smoothed"), linecolor = c("black", "red"))
# grid.draw(lgd)

stat_season <- function(INPUT, brks){
    d_org <- as.data.table(INPUT[c("t", "y", "w")])
    d_fit <- brks$whit %>% .[,.SD,.SDcols=c(1, ncol(.))] %>% set_colnames(c("t", "ypred"))
    d <- merge(d_org, d_fit, by = "t")

    stat <- with(d, GOF_extra2(y, ypred))# %>% as.list()
    stat$nseason <- nrow(brks$dt)

    stat$Roughness = stat$Rg_norm_by_obs
    # str_title <- sprintf("[%s] IGBP = %s, %s, lat = %.2f", sitename, IGBP_name, stat_str, lat)
    # str_title <- paste(titlestr, stat_txt)
    # NSE <- stat$NSE
    # cv  <- stat$cv
    stat
}

plot_season <- function(INPUT, brks, plotdat, ylu, IsOnlyPlotbad = FALSE){
    stat <- stat_season(INPUT, brks)
    stat_txt  <-
        # stat[c("R2", "NSE", "sim.cv", "obs.cv")] %>%
        stat[c("R2", "Bias", "RMSE", "Roughness")] %>% unlist() %>%
        {paste(names(.), round(., 3), sep = "=", collapse = ", ")}

    # if (NSE < 0 | (cv < 0.1 & NSE < 0.1)) {}
    # if(IsPlot && (NSE < 0 && cv < 0.2)){
    if (IsOnlyPlotbad && stat['NSE'] < 0.3) return()

    t  <- brks$whit$t
    dt <- brks$dt
    zs <- dplyr::select(brks$whit,dplyr::matches("ziter*"))
    ypred <- last(zs)

    # if (missing(xlim))
    xlim <- c(first(brks$dt$beg), last(brks$dt$end))

    xlim <- c("2000-01-01", "2017-12-31") %>% ymd()
    at   <- (seq(2000, 2018, 3)*10000+0101) %>% ymd()
    # 7.1 PLOT CURVE FITTING TIME-SERIES
    # need to plot outside, because y, w have been changed.
    plot_input(plotdat, xlim = xlim, xaxt="n")
    axis(1, at = at, labels = seq(2000, 2018, 3))

    colors <- c("red", "blue", "green")
    iters  <- ncol(zs)
    if (iters < 3) colors <- c("red", "green")

    for (i in 1:iters){
        lines(t, zs[[i]], col = colors[i], lwd = 2)
    }

    # 7.2 plot break points
    # points(dt$peak, dt$y_peak, pch=20, cex = 1.8, col="red")
    # points(dt$beg , dt$y_beg , pch=20, cex = 1.8, col="blue")
    # points(dt$end , dt$y_end , pch=20, cex = 1.8, col="blue")

    if (!missing(ylu)) abline(h=ylu, col="red", lty=2) # show ylims
    legend('topleft', stat_txt, adj = c(0.05, 0.8), bty='n',
        cex = 1, text.col = "black")
}

# FUNCTION END ------------------------------------------------------------

adj.param = TRUE
version   = paste0("v3_",
                   ifelse(adj.param, "adjparam", "noadjparam"))

is_flux  <- F

if (is_flux){
    # about 2 minutes
    source("test/07_whit/dat_flux&cam_phenofit.R") # load data at flux sites
    df <- df_org # Get data
}else{
    # about 1 hour
    load("data_test/whit_lambda/MOD13A1_st_1e3_20180731.rda")
}
dir_flux <- ifelse(is_flux, "flux/", "")

## examples
sites    <- unique(df$site) %>% set_names(., .)
sitename <- sites[1]

# 1.1 lambda formula coefs
noise_percs = c(0.1, 0.3, 0.5)
noise_perc  = 0 # default is zero

methods  <- c("smooth_wHANTS", "smooth_wSG", "smooth_wWHIT", "smooth_wWHIT")
methods2 <- c("smooth_wHANTS", "smooth_wSG", "wWH", "wWH2")

lst <- list()

# df <- df[t < ymd("2018-01-01")]
################################################################################

date <- fill_missdate()

# for (k in 4){ # 1:nrow(coefs)
# noise_perc <- noise_percs[k]
# df <- select_valid(df, noise_perc = noise_perc)[, 1:10]

k = 4 # k = 4 is corresponding to `grp01_Extend`
pattern <- rownames(coefs)[k]
param   <- as.list(coefs[k, ])

info  <- readRDS("data_test/wWHd_inspect_sites.RDS")

sites <- info$site
sitename <- sites[1]
sitename <- 9971

sp <- st[site == sitename]
k = 4
iters = 1
param   <- as.list(coefs[k, ]); print(str(param))


svg("Fig8_wWH2 vs wWHd one site.svg", 9, 6)
par(mfrow = c(2, 1), mar = c(1.5, 2, 2, 1), mgp = c(3, 0.6, 0),
    oma = c(3.5, 0, 0, 0),
    cex = 1)
a <- rough_fitting(sitename, df, st, .FUN = smooth_wWHIT, lambda = 2, IsPlot = T, iters = iters);
title("(a) wWH2", adj = 0)

b1 <- rough_fitting(sitename, df, st, .FUN = smooth_wWHIT, lambda = NULL, IsPlot = T, iters = iters)
title("(b) wWHd", adj = 0)
mtext("Time", side =1, line = 1.6, font = 2, cex = 1.2)
pushViewport(viewport(x = 0.5, y = 0.04, width=0.8, height=0.1))
# grid.rect(gp=gpar(fill="white"))
grid.draw(lgd)
dev.off()

# b2 <- rough_fitting(sitename, df, st, .FUN = smooth_wWHIT, lambda = NULL, T, IsOptim_lambda = T, iters = iters)
# title("wWHd_optim")
# info <- list(wWH2 = a$GOF,
#              wWHd_v0 = b1$GOF,
#              wWHd_latest = b2$GOF) %>%
#     melt_list("meth") %>%
#     # .[iter == "iter1"] %>%
#     .[order(type, iter), ]
# print(info)
