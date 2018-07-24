# source('test/stable/load_pkgs.R')
library(Matrix)
library(plyr)
library(data.table)
library(tidyverse)

library(magrittr)
library(lubridate)
library(purrr)
library(zoo)

library(Ipaper)
library(phenofit)
# library(plotly)
library(devtools)
library(Cairo)
library(jsonlite)
library(openxlsx)
library(pbmcapply)
library(MASS)
library(broom)

fontsize  <- 14

qc_values <- c("0", "1", "2", "3")
qc_levels <- c("good", "margin", "snow/ice", "cloud")
qc_colors <- c("grey60", "#00BFC4", "#F8766D", "#C77CFF") %>% set_names(qc_levels)
qc_shapes <- c(19, 15, 4, 17) %>% set_names(qc_levels)

## fluxnet2015 and phenocam dataset were used to test phenofit
dir_data <- "data_test/"

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

siteorder <- function(sites){ factor(sites) %>% as.numeric() }

get_phenofit <- function(sitename, df, st, prefix_fig = 'phenofit_v3'){
    d     <- df[site == sitename, ] # get the first site data
    sp    <- st[site == sitename, ] # station point
    titlestr <- with(sp, sprintf('[%03d,%s] %s, lat = %5.2f, lon = %6.2f',
                                     ID, site, IGBPname, lat, lon))
    file_pdf <- sprintf('Figure/%s_[%03d]_%s.pdf', prefix_fig, sp$ID[1], sp$site[1])

    tryCatch({
        dnew  <- add_HeadTail(d)
        # 1. Check input data and initial parameters for phenofit
        INPUT <- check_input(dnew$t, dnew$y, dnew$w, maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
        INPUT$y0 <- dnew$y
        IsPlot   <- FALSE # for brks
        # 2. The detailed information of those parameters can be seen in `season`.
        lambda <- init_lambda(INPUT$y)#*2
        # brks   <- season(INPUT, nptperyear,
        #                FUN = whitsmw2, wFUN = wFUN, iters = 2,
        #                lambda = lambda,
        #                IsPlot = IsPlot, plotdat = d,
        #                south = d$lat[1] < 0,
        #                rymin_less = 0.6, ymax_min = ymax_min,
        #                max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5) #, ...
        # get growing season breaks in a 3-year moving window
        brks2 <- season_3y(INPUT, nptperyear, south = sp$lat[1] < 0, FUN = whitsmw2,
                           IsPlot = IsPlot, print = print, partial = F)

        # 3. curve fitting
        fit  <- curvefits(INPUT, brks2, lambda =lambda,
                          methods = c("AG", "zhang", "beck", "elmore"), #,"klos",, 'Gu'
                          nptperyear = nptperyear, debug = F, wFUN = wTSM,
                          nextent = 2, maxExtendMonth = 3, minExtendMonth = 1,
                          qc = as.numeric(dnew$SummaryQA), minPercValid = 0.2,
                          print = print)
        fit$INPUT   <- INPUT
        fit$seasons <- brks2

        # svg("Figure1_phenofit_curve_fitting.svg", 11, 7)
        Cairo::CairoPDF(file_pdf, 11, 6) #
        # grid::grid.newpage()
        plot_phenofit(fit, d, titlestr) %>% grid::grid.draw()# plot to check the curve fitting
        dev.off()

        # temp <- ExtractPheno(fit$fits$ELMORE[1:5], IsPlot = T) # check extracted phenology
        ## 3.2 Get GOF information
        stat  <- ldply(fit$fits, function(fits_meth){
            ldply(fits_meth, statistic.phenofit, .id = "flag")
        }, .id = "meth")
        fit$stat <- stat

        # 4. extract phenology, check extracted phenology extraction for one method.
        # ratio = 1.15
        # file <- "Figure5_Phenology_Extraction_temp.pdf"
        # cairo_pdf(file, 8*ratio, 6*ratio)
        # temp <- ExtractPheno(fit$fits$ELMORE[2:6], IsPlot = T)
        # dev.off()
        # file.show(file)

        # pheno: list(p_date, p_doy)
        p <- lapply(fit$fits, ExtractPheno)
        pheno  <- map(p, tidyFitPheno, origin = INPUT$t[1]) %>% purrr::transpose()

        fit$pheno  <- pheno
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

    years <- 2000:2017
    doy   <- seq(1, 366, 16)
    date    <- sprintf("%4d%03d", rep(years, each = 23), doy)[-(1:3)] %>% parse_date_time("%Y%j") %>% date()
    df$date <- date # 20180705, Simon has fixed image missing
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
FigsToPages <- function(ps, lgd, ylab.right, file, width = 10, height){
    nrow  <- 6
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
