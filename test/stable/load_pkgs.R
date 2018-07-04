# source('test/stable/load_pkgs.R')
library(Matrix)
library(plyr)
library(data.table)
library(tidyverse)

library(magrittr)
library(lubridate)
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
                   "URB", "CNV", "SNOW", "BSV", "water", "UNC")
# MCD12Q1.005 land cover 0-16 and 254, IGBP scheme
# https://code.earthengine.google.com/dataset/MODIS/051/MCD12Q1
IGBPnames_005 <- c("water", "ENF", "EBF", "DNF", "DBF", "MF" , "CSH",
                   "OSH", "WSA", "SAV", "GRA", "WET", "CRO",
                   "URB", "CNV", "SNOW", "BSV", "UNC")

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


siteorder <- function(sites){ factor(sites) %>% as.numeric() }

#' Use exact date not image date
getRealDate <- function(df){
    df[, `:=`(date = ymd(date), year = year(date), doy = as.integer(yday(date)))]
    df[is.na(DayOfYear), DayOfYear := doy] # If DayOfYear is missing

    # In case of last scene of a year, doy of last scene could in the next year
    df[abs(DayOfYear - doy) >= 300, t := as.Date(sprintf("%d-%03d", year+1, DayOfYear), "%Y-%j")] # last scene
    df[abs(DayOfYear - doy) <  300, t := as.Date(sprintf("%d-%03d", year  , DayOfYear), "%Y-%j")]
    df
}

qc_values <- c("0", "1", "2", "3")
qc_levels <- c(" good", " margin", " snow/ice", " cloud")

#' Tidy phenofit INPUT data from raw original data exported from gee
#' This function only suit for MODIS VI products (e.g. MOD13A1, MOD13A2, ...)
tidyMOD13INPUT_gee <- function(infile, outfile){
    df   <- fread(infile)
    df %<>% getRealDate()

    # Initial weights
    df[, w := qc_summary(SummaryQA, wmin = 0.2)]
    # Remap SummaryQA factor level, plot_phenofit use this variable. For other
    # remote sensing data without `SummaryQA`, need to modify `plot_phenofit`
    if ('SummaryQA' %in% colnames(df)){
        df$SummaryQA %<>% factor() %>% mapvalues(qc_values, qc_levels)
    }

    df <- df[, .(site, y = EVI/1e4, t, w, date, SummaryQA
                 # IGBPcode,
                 # IGBPname = as.factor(IGBPnames[IGBPcode]),
                 )]
    df
    # merge coordinate info
    # df <- merge(df, st[, .(Id = site, lat, lon = long, IGBPname = IGBP)], by = "Id")
    # fwrite(df, outfile)
}

############################# GEE WHITTAKER ####################################
#' This function is only used to read gee_phenofit whittaker result.
read_whit.gee <- function(file){
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

    sites <- map_chr(lst, "id") %>% str_extract(".*(?=_)")

    df <- set_names(data, sites) %>% melt_list(var.name = "site")
    df
}

#' read gee_phenofit whittaker from multiple json files
#' 
#' @examples
#' indir <- "D:/Document/GoogleDrive/phenofit/data/gee_phenofit/v2/"
#' files <- dir(indir, "*.geojson", full.names = T)
#' df_gee <- read_whits.gee(files)
read_whits.gee <- function(files){
    lst   <- llply(files, read_whit.gee, .progress = "text")
    df    <- do.call(rbind, lst) %>% {.[order(site), ]}

    years <- 2000:2017
    doy   <- seq(1, 366, 16)
    date    <- sprintf("%4d%03d", rep(years, each = 23), doy)[-(1:3)] %>% parse_date_time("%Y%j") %>% date()
    df$date <- date # 20180705, Simon has fixed image missing
    df
}

get_phenofit <- function(sitename, df, st, prefix_fig = 'phenofit_v3'){
    d     <- df[site == sitename, ] # get the first site data
    sp    <- st[site == sitename, ] # station point
    titlestr <- with(sp, sprintf('[%03d_%s] %s, lat = %5.2f, lon = %6.2f',
                                     ID, site, IGBPname, lat, lon))
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
        file_pdf = sprintf('Figure/%s_%s.pdf', prefix_fig, d$site[1])
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
