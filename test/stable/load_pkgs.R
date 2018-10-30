# source('test/stable/load_pkgs.R')
suppressMessages({
    library(Matrix)
    library(plyr)
    library(data.table)
    library(tidyverse)
    library(broom)

    library(magrittr)
    library(lubridate)
    library(purrr)
    library(zoo)

    library(devtools)
    library(jsonlite)
    library(openxlsx)
    # library(pbmcapply)
    # library(MASS)

    ## visualization pkgs
    library(grid)
    library(gridExtra)
    library(Cairo)
    # library(plotly)

    ## self pkgs
    library(Ipaper)
    library(phenofit)
})

fontsize  <- 14

subl <- function() system('subl .')

qc_values <- c("0", "1", "2", "3")
qc_levels <- c("good", "margin", "snow/ice", "cloud")
qc_colors <- c("grey60", "#00BFC4", "#F8766D", "#C77CFF") %>% set_names(qc_levels)
qc_shapes <- c(19, 15, 4, 17) %>% set_names(qc_levels)

## fluxnet2015 and phenocam dataset were used to test phenofit
dir_data <- "data_test/"
dir_flux <- "F:/Github/data/flux/"

file_st_cam  <- paste0(dir_data, "st_phenocam133.csv")
file_st_flux <- paste0(dir_data, "st_phenoflux166.csv")

file_flux <- paste0(dir_data, "phenoflux166_MOD13A1_INPUT.csv")
file_cam  <- paste0(dir_data, "phenocam133_MOD13A1_INPUT.csv")

##
if (.Platform$OS.type == "unix"){
    dir_climate <- "/OSM/CBR/CoRE/working/timeseries/Climate/"
    dir_flush   <- "/flush1/kon055/"
} else{
    dir_climate <- "//clw-03-cdc.it.csiro.au/OSM_CBR_CoRE_working/timeseries/Climate/"
    dir_flush   <- "//braggflush1/flush1/kon055/"
}

# MCD12Q1.006 land cover 1-17 and 255, IGBP scheme
IGBPnames_006 <- c("ENF", "EBF", "DNF", "DBF", "MF" , "CSH",
                   "OSH", "WSA", "SAV", "GRA", "WET", "CRO",
                   "URB", "CNV", "SNO", "BSV", "water", "UNC")
# MCD12Q1.005 land cover 0-16 and 254, IGBP scheme
# https://code.earthengine.google.com/dataset/MODIS/051/MCD12Q1
IGBPnames_005 <- c("water", "ENF", "EBF", "DNF", "DBF", "MF" , "CSH",
                   "OSH", "WSA", "SAV", "GRA", "WET", "CRO",
                   "URB", "CNV", "SNO", "BSV", "UNC")


metric_spring <- c('TRS1.sos', 'TRS2.sos', 'TRS5.sos', 'TRS6.sos', 'DER.sos',
                   'GU.UD', 'GU.SD', 'ZHANG.Greenup', 'ZHANG.Maturity')
metric_autumn <- c('DER.eos', 'TRS6.eos', 'TRS5.eos', 'TRS2.eos', 'TRS1.eos',
                   'GU.DD', 'GU.RD', 'ZHANG.Senescence', 'ZHANG.Dormancy')
metrics <- c(metric_spring, "DER.pop", metric_autumn)
#' phase: 'spring', 'pop', 'autumn'
metric_phase <- function(metric){
    phase <- rep("pop", length(metric))
    phase[metric %in% metric_spring] <- "spring"
    phase[metric %in% metric_autumn] <- "autumn"
    phase %>% factor(c("spring", "pop", "autumn"))
}

#' nth_max
#'
#' The nth maximum value
#' @examples
#' x <- c(12.45,34,4,0,-234,45.6,4)
#' nth_max(x)
nth_max <- function(x, n = 2){
    len <- length(x)

    if (sum(!is.na(x)) <= n){
        min(x, na.rm = T)
    } else {
        sort(x, decreasing = T)[n]
    }
    # i   <- len-n+1
    # sort(x, partial=i, na.last=T)[i]
}

fix_null <- function(x, default = NA){
    I <- sapply(x, is.null)
    x[I] <- default
    x
}

clamp_min <- function(x, value = 0){
    x[x < value] <- value
    x
}

#' clamp
#' clamp values in the range of `lims`
clamp <- function(x, lims = c(0, 1)){
    x[x < lims[1]] <- lims[1]
    x[x > lims[2]] <- lims[2]
    x
}

# save pdf just like `ggsave`
write_fig <- function(p, file = "Rplot.pdf", width = 10, height = 5, show = T, res = 300){
    if (missing(p)) p <- last_plot()

    if ("grob" %in% class(p)) {
        FUN <- grid::grid.draw
    } else{
        FUN <- base::print
    }

    file_ext <- str_extract(basename(file), "(?<=\\.).{1,4}$")

    param <- list(file, width = width, height = height)
    if (file_ext == "pdf"){
        devicefun <- cairo_pdf # Cairo::CairoPDF #
    } else if (file_ext == "svg"){
        devicefun <- svg
    } else {
        if (file_ext %in% c("tif", "tiff")){
            devicefun <- tiff
        } else if (file_ext == "png") {
            devicefun <- Cairo::CairoPNG
        }
        param %<>% c(list(units = "in", res = res, compression = "lzw")) #, dpi = 300
    }

    # print(FUN)
    # Cairo::CairoPDF, if only one figure cairo_pdf is the best
    do.call(devicefun, param)
    FUN(p)
    dev.off()
    if (show) file.show(file)
}

g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}


readRDS_tidy <- function(file){
    file <- gsub("file:///", "", file)
    readRDS(file)
}

ddply_dt <- function(d, j, by){
    if (is.quoted(by)) by <- names(by)

    byname  <- paste(by, collapse = ", ")
    operate <- j[[1]] %>% deparse()

    eval(parse(text = sprintf("res <- d[, .(res = list(%s)), .(%s)]", operate, byname)))
    # print(res)
    res$res %>% do.call(rbind, .) %>% data.table() %>% cbind(res[, ..by], .)
}


#' @return
#' Rg   : normalized
#' Rg_0 : no normalized
GOF_extra <- function(Y_obs, Y_pred){
   
    # roughness
    Roughness <- function(Y_pred){
        temp = diff(Y_pred)^2 %>% .[!is.na(.)]
        if (is_empty(temp)) return(NA_real_)
        Rg = sqrt(sum(temp)/length(temp))
        Rg
    }

    znorm <- function(Y_pred, Y_norm){
        min <- min(Y_norm, na.rm = T)
        max <- max(Y_norm, na.rm = T)
        A   <- max - min

        if (is.finite(A)){
            Yz <- (Y_pred  - min) / A  # normalized
        } else {
            Yz <- Y_pred * NA
        }
        return(Yz)
    }

    ## roughness second definition
    if (sum(is.finite(Y_pred)) <= 1){
        Rg              <- NA_real_
        Rg_norm_by_pred <- NA_real_
        Rg_norm_by_obs  <- NA_real_
        cv              <- NA_real_
        ACF <- NULL
    } else {
        # for different methods, using the same `min` and `max` value
        Y_norm_by_pred  <- znorm(Y_pred, Y_pred)
        Y_norm_by_obs   <- znorm(Y_pred, Y_obs)

        Rg              <- Roughness(Y_pred)
        Rg_norm_by_pred <- Roughness(Y_norm_by_pred)
        Rg_norm_by_obs  <- Roughness(Y_norm_by_obs)
        cv <- cv_coef(Y_pred)[['cv']]

        # the autocorrelation of residuals
        ACF = acf(Y_pred - Y_obs, lag.max = 10, plot = F, na.action = na.pass)$acf[,,1][-1]
    }

    c(Rg = Rg,
        Rg_norm_by_obs  = Rg_norm_by_obs,
        Rg_norm_by_pred = Rg_norm_by_pred, 
        cv = cv, 
        acf = list(list(ACF))) #GOF(Y_obs, Y_pred),
}

GOF_extra2 <- function(Y_obs, Y_pred){
    gof_1 <- GOF(Y_obs, Y_pred)
    gof_2 <- GOF_extra(Y_obs, Y_pred)

    res <- c(gof_1, gof_2)
    res
}

# function to separate data to steps of x, obtain 95 quantile value for smooth
upper_envelope <- function(x, y, step = 0.2, alpha = 0.95){
    xrange <- range(x, na.rm = T)

    brks <- seq(xrange[1], xrange[2], by = step)
    n    <- length(brks)
    xmid <- (brks[-n] + brks[-1])/2

    brks[n] <- Inf

    res <- numeric(n-1)*NA_real_

    for (i in 1:(n-1)){
        val_min <- brks[i]
        val_max <- brks[i+1]

        I <- x >= val_min & x < val_max
        res[i] <- quantile(y[I], alpha, na.rm = T)
    }

    data.table(x = xmid, y = res)
}

siteorder <- function(sites){ factor(sites) %>% as.numeric() }

# Get statistics information of phenofit
getPheno_phenofit <- function(fit, INPUT, brks, d, file_pdf, titlestr){
    ## check the curve fitting parameters
    params <- getparam(fit)
    # print(str(params, 1))
    # print(params$AG)

    ## Get GOF information
    stat  <- ldply(fit$fits, function(fits_meth){
        ldply(fits_meth, statistic.phenofit, .id = "flag")
    }, .id = "meth")
    fit$stat <- stat
    # print(head(stat))

    # pheno: list(p_date, p_doy)
    p <- lapply(fit$fits, ExtractPheno)
    pheno  <- map(p, tidyFitPheno, origin = INPUT$t[1]) %>% purrr::transpose()
    fit$pheno  <- pheno

    fit$INPUT   <- INPUT
    fit$seasons <- brks

    if (IsPlot){
        # svg("Figure1_phenofit_curve_fitting.svg", 11, 7)
        Cairo::CairoPDF(file_pdf, 11, 6) #
        # grid::grid.newpage()
        plot_phenofit(fit, d, titlestr) %>% grid::grid.draw()# plot to check the curve fitting
        dev.off()
    }
    return(fit)
}

#' Get GPPobs phenology dat
#' @examples
#' get_phenofit(df, st, brks_lst, sites, wFUN = 'wTSM')
get_phenofit_GPPobs <- function(sitename,
    df, st, brks_lst, sites, wFUN = 'wTSM',
    prefix_fig = 'phenofit_0.1.6', IsPlot = T)
{
    # sitename <- sites[i]
    i <- grep(sitename, sites)
    brks2 <- brks_lst[[i]]

    ## 1. prepare inputs
    d   <- df[site == sitename, .(t = date, GPP_DT, GPP_NT, w = 1)] #%T>% plotdata(365)
    d$y <- rowMeans(d[, .(GPP_DT, GPP_NT)], na.rm = T)
    d[y < 0, y := 0] # for GPP_NT

    sp       <- st[site == sitename, ]
    titlestr <- with(sp, sprintf('[%03d,%s] %s, lat = %5.2f, lon = %6.2f',
                                     ID, site, IGBPname, lat, lon))
    file_pdf <- sprintf('Figure/%s_[%03d]_%s.pdf', prefix_fig, sp$ID[1], sp$site[1])

    # parameters for season_3y
    INPUT <- getINPUT_GPPobs(df, st, sitename)

    ## 3. Get daily curve fitting result
    wFUN <- get(wFUN)
    fit  <- curvefits(INPUT, brks2,
                      methods = c("AG", "zhang", "beck", "elmore"), #,"klos",, 'Gu'
                      debug = F,
                      wFUN = wFUN,
                      nextent = 5, maxExtendMonth = 2, minExtendMonth = 1/3,
                      # qc = as.numeric(dnew$SummaryQA),
                      minPercValid = 0.2,
                      print = FALSE)

    fit <- getPheno_phenofit(fit, INPUT, brks2, d, file_pdf, titlestr)

    return(fit)
}

#' get_phenofit
#'
#' @param nextent Suggest as ceiling(nptperyear/12*1.5). Daily GPPObs, 5; 16-day
#' VI, 2.
#' @param isVarLambda If `isVarLambda` = true, lambda will be automatically assigned
#' with `init_lambda` in different years.
get_phenofit <- function(sitename, df, st, prefix_fig = 'phenofit_v0.1.6', IsPlot = F,
    nptperyear = 23,
    brks = NULL,
    ypeak_min = 0.05,
    wFUN_season = wTSM, wFUN_fit = wFUN_season,
    lambda = NULL, isVarLambda = FALSE,
    nextent = 5, maxExtendMonth = 3, minExtendMonth = 1)
{
    d     <- df[site == sitename, ] # get the first site data
    sp    <- st[site == sitename, ] # station point
    south <- sp$lat < 0
    if (length(south) == 0) south <- F

    lat <- lon <- NA; IGBPname <- ""
    titlestr <- with(sp, sprintf('[%03d,%s] %s, lat = %.2f, lon = %.2f',
                                     ID, site, IGBPname, lat, lon))
    file_pdf <- sprintf('Figure/%s_[%03d]_%s.pdf', prefix_fig, sp$ID[1], sp$site[1])

    tryCatch({
        dnew  <- add_HeadTail(d, south, nptperyear)
        # 1. Check input data and initial parameters for phenofit
        INPUT <- check_input(dnew$t, dnew$y, dnew$w, nptperyear,
            maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
        INPUT$y0    <- dnew$y
        INPUT$south <- south

        IsPlot2   <- FALSE # for brks
        # 2. The detailed information of those parameters can be seen in `season`.
        # If lambda initialized at here, yearly variation will be not include.
        if (is.null(lambda)) lambda <- init_lambda(INPUT$y)#*2
        if (isVarLambda) lambda <- NULL

        # brks   <- season(INPUT, nptperyear,
        #                FUN = whitsmw2, wFUN = wFUN, iters = 2,
        #                lambda = lambda,
        #                IsPlot = IsPlot, plotdat = d,
        #                south = d$lat[1] < 0,
        #                rymin_less = 0.6, ymax_min = ymax_min,
        #                max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5) #, ...
        ## get growing season breaks in a 3-year moving window
        brks2 <- season_3y(INPUT, south = south,
            wFUN = wFUN_season, rFUN = wWHIT, iters = 2,
            lambda = lambda,
            minpeakdistance = nptperyear/6,
            MaxPeaksPerYear = 3,
            MaxTroughsPerYear = 4,
            ypeak_min = ypeak_min,
            IsPlot = F, print = FALSE, IsOnlyPlotbad = F)
        if (!is.null(brks)) brks2$dt <- brks$dt

        # 3. curve fitting
        fit  <- curvefits(INPUT, brks2,
                          methods = c("AG", "zhang", "beck", "elmore"), #,"klos",, 'Gu'
                          debug = F,
                          wFUN = wFUN_fit,
                          nextent = nextent, maxExtendMonth = maxExtendMonth, minExtendMonth = minExtendMonth,
                          qc = as.numeric(dnew$QC_flag),
                          minPercValid = 0.2,
                          print = FALSE)

        fit <- getPheno_phenofit(fit, INPUT, brks2, d, file_pdf, titlestr)
        return(fit)
    }, error = function(e){
        message(sprintf('[e]: %s, %s', titlestr, e$message))
    })
}


# rename phenofit phenology metrics names
fix_level <- function(x){
    phenophase <- c(
        'TRS1.sos', 'TRS2.sos', 'ZHANG.Greenup', 'GU.UD',
        'TRS5.sos', 'TRS6.sos', 'DER.sos',
        'ZHANG.Maturity','GU.SD', 'DER.pop',
        'ZHANG.Senescence', 'GU.DD',
        'TRS5.eos', 'TRS6.eos','DER.eos',
        rev(c('TRS1.eos', 'TRS2.eos', 'ZHANG.Dormancy', 'GU.RD')))
    phenophase_spl <- c(
        'TRS1.SOS', 'TRS2.SOS', 'Greenup', 'UD',
        'TRS5.SOS', 'TRS6.SOS', 'DER.SOS',
        'Maturity', 'SD', 'POP',
        'Senescence', 'DD',
        'TRS5.EOS', 'TRS6.EOS','DER.EOS',
        rev(c('TRS1.EOS', 'TRS2.EOS', 'Dormancy', 'RD')))
    factor(x, phenophase) %>% mapvalues(phenophase, phenophase_spl)#return
}

############################# GEE WHITTAKER ####################################
#' This function is only used to read gee_phenofit whittaker result.
read_whitMat.gee <- function(file){
    # file: phenoflux166_WH_p2_2008_2011.geojson
    years <- str_extract_all(basename(file), "\\d{4}")[[1]] %>% map_dbl(as.numeric)
    years <- years[1]:years[2]
    doy   <- seq(1, 366, 16)
    date  <- sprintf("%4d%03d", rep(years, each = 23), doy) %>% parse_date_time("%Y%j") %>% date()
    if (years[1] == 2000) date <- date[-(1:3)]

    lst   <- read_json(file)$features
    ncol  <- length(lst[[1]]$properties$array[[1]])
    names <- c("raw", paste0("iter", 1:(ncol-1) ) )

    data <- map(lst, function(l){
        l$properties$array %>% purrr::transpose() %>%
            map(unlist) %>%
            do.call(cbind.data.frame, .) %>%
            set_colnames(names) %>%
            data.table()
    })

    sites <- map_chr(lst, ~.x$properties$site)
    # sites <- map_chr(lst, "id") %>% str_extract(".*(?=_)")
    df <- set_names(data, sites) %>% melt_list(var.name = "site")
    df$date <- date
    df
}


#' read gee_phenofit whittaker from multiple json files
#'
#' @examples
#' indir <- "D:/Document/GoogleDrive/phenofit/data/gee_phenofit/v2/"
#' files <- dir(indir, "*.geojson", full.names = T)
#' df_gee <- read_whits.gee(files)
read_whitMats.gee <- function(files){
    lst   <- llply(files, read_whitMat.gee, .progress = "text")
    df    <- do.call(rbind, lst) %>% {.[order(site), ]}
    # 20180705, Simon has fixed image missing
    df
}

#' read gee whit matrix
#'
readwhitMAT <- function(indir, prefix){
    pattern <- paste0(prefix, "_.*.geojson")
    files   <- dir(indir, pattern, full.names = T)
    df      <- read_whitMats.gee(files)
    df
}

#' FigsToPages
#'
#' Subplots xlab and ylab are unified, and only keet singe one.
#' Currently, only support ggplot figures; And only support arrange figures
#' into rows (nrow*1).
#'
#' @param ps A list of ggplot figure objects. And
#' @param lgd A grid grob object, legend to show in the bottom.
#' @param ylab y label title
#' @param width
#' @param height
#'
#' @param
FigsToPages <- function(ps, lgd, ylab.right, file, width = 10, height, nrow = 6){
    npage <- ceiling(length(ps)/nrow)
    if (missing(height)) height = nrow*1.6

    ylab.left        <- ps[[1]]$labels$y
    ylab.left.color <- ps[[1]]$theme$axis.title.y.left$colour
    ylab.right.color <- ps[[1]]$theme$axis.title.y.right$colour

    params <- list(ncol = 1, padding = unit(1, "line"),
        left  = textGrob(ylab.left , rot = 90,
                         gp=gpar(fontsize=14, col=ylab.left.color)) )

    # parameters for arrangeGrob
    if (!missing(ylab.right))
        params$right = textGrob(ylab.right, rot = 270,
                                gp=gpar(fontsize=14, col=ylab.right.color))

    Cairo::CairoPDF(file, width, height)
    for (i in 1:npage){
        runningId(i)
        I_beg <- (i - 1) * nrow + 1
        I_end <- min(i*nrow, length(ps))

        I  <- I_beg:I_end
        n  <- length(I)

        ps_i <- ps[I]
        for (j in seq_along(I)){
            theme_j <- theme(
                axis.text.x = element_blank(),
                axis.title = element_blank(),
                axis.title.y.right = element_blank(),
                axis.title.y.left  = element_blank()
            )
            if (j == n)
                theme_j <- theme(
                    axis.title.y.right = element_blank(),
                    axis.title.y.left  = element_blank() )
            ps_i[[j]] <- ps_i[[j]] + theme_j
        }

        ps_i  <- c(ps_i, list(lgd))
        nx    <- length(ps_i)

        params$grobs <- ps_i
        params$nrow  <- nx

        if (missing(lgd)){
            params$heights <- c(rep(5, nx - 1), 5.5)
        } else{
            params$heights <- c(rep(5, nx - 2), 5.5, 1)
        }

        g <- do.call(gridExtra::arrangeGrob, params)

        if (i != 1) grid.newpage();
        grid::grid.draw(g)
    }
    dev.off()
    file.show(file)
}


############################### VALIDATION FUNCTIONS ###########################
select_validI <- function(Id, perc = 0.2) {
    n <- length(Id)
    set.seed(100)
    I <- sample(1:n, n*perc)
    unique(Id[I])
}

# only select good points (with the percentage of noise_perc)
select_valid <- function(df, noise_perc = 0.3, group = F){
    df     <- unique(df) # sometimes data is duplicated at begin and end.
    if (class(df$SummaryQA[1]) != "factor") {
        df[, `:=`( SummaryQA = factor(SummaryQA, qc_levels))]
    }
    if (class(df$t[1]) != "Date") df[, `:=`( t = ymd(t))]

    ## 1.1 get range for every site and divide grp
    alpha = 0.05
    d_range <- df[SummaryQA == "good", .(ymin = quantile(y, alpha/2, na.rm = T),
                                         ymax = quantile(y, 1-alpha/2, na.rm = T)), .(site)] %>%
        plyr::mutate(A = ymax - ymin,
                     brk_min = ymin + 0.2*A,
                     brk_max = ymin + 0.6*A) %>% .[, c(1, 4:6)]
    if (!group) d_range <- d_range[, .(site, A)]

    df <- merge(df, d_range, all.x = T) # rm points has no good points
    df[, `:=`(Id = 1:.N, y0 = y, w0 = w, I_valid = 0)]

    if (group){
        df[, grp:=1]
        df[y >= brk_max, grp := 2]
        df[y <= brk_min, grp := 0]
    }

    # 1. select cross-validation points ---------------------------------------
    if (noise_perc > 0){
        ## 1.2 select validation points
        # & grp == 1
        I <- df[SummaryQA == "good", select_validI(Id, noise_perc), .(site)]$V1

        set.seed(I[1])
        desc_perc <- runif(length(I), .05, .5)

        ## 1.3 adjust y and w
        # randomly reduced by 5%–50% with of their amplitude; w set to wmin

        # df_val <- df[I, ]
        df[I, `:=`(w=0, y=y-A*desc_perc, I_valid = 1)]
    }
    return(df)
}

# nptperyear = 46
# # df <- fread("data/lc006/PMLv2_flux112_sgfitw&TSM.csv")
# df <- fread("PMLv2_flux112_CV.csv")
# df$date %<>% ymd
# df <- df[order(site, date), ]
#
# # For each site, remove na values at head and tail
# df <- ddply(df, .(site), function(x){
#     if (all(is.na(x$GPPobs))) return(NULL)
#     I <- which(!is.na(x$GPPobs)) %>% {first(.):last(.)}
#     x[I, ]
# }) %>% as.data.table()
# df$YYYY <- df$year
# df[lat < 0, YYYY := year + as.integer(date > ymd(sprintf("%d0701", year))) - 1L];
#
# sites_rm1 <- df[lat < 0, unique(lat<0), .(site)]$site
# sites_rm2 <- c("GF-Guy", "BR-Sa3", "US-Whs")
# sites_rm  <- union(sites_rm1, sites_rm2)
## df <- df[!(site %in% sites_rm), ]
#
# save(df, file = "phenofit_flux90_INPUTS.rda")

# source("R/phenofit.R")
#
# source("R/derivs.R")
# source("R/optimDL.R")
# source("R/doubleLogistics.R")
# source("R/doubleLogistics_fitting.R")
# source("R/PhenoExtract.R")
# source("R/PhenoExtract_main.R")
# source('R/PhenoBrks.R', encoding = "utf-8")
# source('R/pkg_smooth.R', encoding = "utf-8")
# source("F:/Github/PML_v2/fluxsites_tidy/R/mainfunc/load_pkgs.R", encoding = "utf-8")
# stations212 <- fread("C:/Users/kon055/Google Drive/Github/data/phenology/station/flux-212.txt")

# tidy_pheno <- function(RES){
#     id.vars <- colnames(RES[[1]]$pheno$doy$AG)
#     df <- map(rm_empty(RES), ~.x$pheno$doy) %>%
#         rm_empty() %>%
#         melt(id.vars = id.vars) %>%
#         set_names(c(id.vars, "meth", "site")) %>% as.data.table()
#     return(df)
# }

# # re-calculate phenology of every site
# recal_pheno.site <- function(fit){
#     # 3. phenology
#     p <- lapply(fit$fits, getFits_pheno)
#     # pheno: list(p_date, p_doy)
#     fit$pheno  <- map(p, tidyFits_pheno, origin = fit$INPUT$t[1]) %>% purrr::transpose()
#     return(fit)
# }

# plotsites <- function(fits, file = 'Fig3_GPP_phenofit_flux112_v13.pdf'){
#     Cairo::CairoPDF(file, width = 10, height = 7)
#     sites <- names(fits)

#     for (i in seq_along(fits)){
#         runningId(i)
#         site <- sites[i]
#         tryCatch({
#             p <- plot_phenofit(fits[[i]]) + ggtitle(site)
#             print(p)
#         }, error = function(e){
#             message(sprintf("%s:%s", site, e$message))
#         })
#     }
#     dev.off()
# }
