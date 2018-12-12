# source('inst/shiny/check_season/global.R')
library(phenofit)
library(shiny)
# library(DT)
library(data.table)
library(magrittr)

library(plyr)
library(purrr)

load('data/phenoflux115.rda')

# load('data/phenoflux115_ET&GPP&VI.rda')
# load('inst/shiny/check_season/data/phenoflux_115.rda')
# load('inst/shiny/check_season/data/ET&GPP&VI_flux115.rda')
# sites <- sort(sites)

# Generate DT::datatable
DT_datatable <- function(
    df,
    pageLength = 10,
    columnDefs = list(list(className = 'dt-center')), ...){

    DT::datatable(df, options = list(
        # autoWidth = TRUE,
        # columnDefs = list(list(width = '10px', targets = c(4:10)))
        searching = FALSE, lengthChange = FALSE,
        pageLength = pageLength,
        columnDefs = columnDefs, ...
    ))
}


# tidy GPP data
tidy_fluxGPP <- function(){
    gpp <- rowMeans(df[, .(GPP_DT, GPP_NT)], na.rm = T)
    gpp[gpp < 0] <- 0

    df <<- data.table(site = df$site, t = df$date, y = gpp) # , w = 1 (optional)
    # df
}

updateINPUT <- function(input){
    if (input$file_type == '.rda | .RData') {
        load(file_rda)
    } else if (input$file_type == 'text'){
        if (file_veg != "") {
            df    <<- fread(file_veg)
            sites <<- unique(df$site) %>% sort()

            if (file_site == ""){
                st <<- data.table(ID = seq_along(sites), site = sites, lat = 30)
            }
        }
        if (file_site != "") st <- fread(file_site)
    }

    list(df = df, st = st, sites = sites)
}


getDf.site  <- function(df, sitename){
    dplyr::select(df[site == sitename, ], dplyr::matches("t|y|w"))
    #%T>% plotdata(365)
}


getINPUT.site <- function(df, st, sitename){
    sp       <- st[site == sitename]
    south    <- sp$lat < 0
    titlestr <- with(sp, sprintf("[%3d] %s, %s, lat = %.2f", ID, site, IGBP, lat))

    d     <-  df[site == sitename, ]#%T>% plotdata(365)
    d_new <- add_HeadTail(d, south = south)
    INPUT <- do.call(check_input, d_new)

    INPUT$south    <- south
    INPUT$titlestr <- titlestr
    # list(INPUT = INPUT, plotdat = d)
    INPUT
}


check_season <- function(INPUT,
                         FUN_season = c("season", "season_3y"),
                         rFUN = "wWHIT",
                         wFUN = "wTSM",
                         lambda = 1000,
                         iters = 3,
                         IsPlot = F, ...) {
    # sitename <- "US-ARM" # "FR-LBr", "ZA-Kru", "US-ARM"

    FUN_season <- get(FUN_season[1])
    wFUN       <- get(wFUN)
    res  <- FUN_season(INPUT, south = INPUT$south,
                        rFUN = get(rFUN),
                        wFUN = wFUN,
                     IsPlot = IsPlot,
                     lambda = lambda,
                     iters = iters,
                     minpeakdistance = 30,
                     MaxPeaksPerYear = 3,
                     MaxTroughsPerYear = 4,
                     ypeak_min = 1, ...,
                     IsOnlyPlotbad = FALSE
    )

    if (IsPlot){
        abline(h = 1, col = "red")
        title(INPUT$titlestr)
    }
    return(res)
}

# plot_data <- function(d, title){
#     par(setting)
#     do.call(check_input, d) %>% plotdata()
#     mtext(title, side = 2, line = 2, cex = 1.3, font = 2)
# }

################################################################################
## global parameters for check_season
nptperyear <- 365
setting    <- list(mar = c(2, 3, 1, 1), mgp = c(1.2, 0.6, 0))
fig.height <- 225
par(setting)

param_step <- 0.1

tidy_fluxGPP() # tidy df_GPP

# https://stackoverflow.com/questions/48592842/show-inf-in-dtdatatable
options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))

# sites_grp  <- list(sites_sg, sites_mg) %>%
#     set_names(c("single season", "multiple season"))

# threshold_max = 0.1,
# threshold_min = 0,
# rytrough_max = 0.6,
