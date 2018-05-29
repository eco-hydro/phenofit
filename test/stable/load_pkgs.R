# source('test/load_pkgs.R')

library(Matrix)
library(plyr)
library(data.table)
library(tidyverse)

library(magrittr)
library(lubridate)
library(zoo)

library(Ipaper)
library(phenofit)
library(plotly)
library(devtools)
library(Cairo)
library(jsonlite)
library(openxlsx)

# MCD12Q1.006 land cover 1-17, IGBP scheme
IGBPnames <- c("ENF", "EBF", "DNF", "DBF", "MF" , "CSH", 
              "OSH", "WSA", "SAV", "GRA", "WET", "CRO", 
              "URB", "CNV", "SNOW", "BSV", "water", "UNC")

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
sites_rm2 <- c("GF-Guy", "BR-Sa3", "US-Whs")
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

tidy_pheno <- function(RES){
    id.vars <- colnames(RES[[1]]$pheno$doy$AG)
    df <- map(rm_empty(RES), ~.x$pheno$doy) %>%
        rm_empty() %>%
        melt(id.vars = id.vars) %>%
        set_names(c(id.vars, "meth", "site")) %>% as.data.table()
    return(df)
}

# re-calculate phenology of every site
recal_pheno.site <- function(fit){
    # 3. phenology
    p <- lapply(fit$fits, getFits_pheno)
    # pheno: list(p_date, p_doy)
    fit$pheno  <- map(p, tidyFits_pheno, origin = fit$INPUT$t[1]) %>% purrr::transpose()
    return(fit)
}

plotsites <- function(fits, file = 'Fig3_GPP_phenofit_flux112_v13.pdf'){
    Cairo::CairoPDF(file, width = 10, height = 7)
    sites <- names(fits)

    for (i in seq_along(fits)){
        runningId(i)
        site <- sites[i]
        tryCatch({
            p <- plot_phenofit(fits[[i]]) + ggtitle(site)
            print(p)
        }, error = function(e){
            message(sprintf("%s:%s", site, e$message))
        })
    }
    dev.off()
}

GOF2 <- function(RE){
    RE <- RE[!is.na(RE)] #discard(RE, is.na)
    n_sim <- length(RE)

    if (n_sim <= 1){
        return(c(Bias = NA, MAE = NA,RMSE = NA, n_sim = NA)) #R = R,
    }
    Bias  <- mean(RE)                                        # bias
    MAE   <- mean(abs(RE))                                   # mean absolute error
    RMSE  <- sqrt(sum((RE)^2) / length(RE))               # root mean sqrt error
    # NASH  <- 1  - sum((RE)^2) / sum((Y_obs - mean(Y_obs))^2) # NASH coefficient

    return(c(Bias = Bias, MAE = MAE,RMSE = RMSE, n_sim = n_sim)) #R = R,
}
