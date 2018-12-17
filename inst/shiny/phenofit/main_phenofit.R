# source('main_phenofit.R')
suppressMessages({
    library(phenofit)
    library(shiny)
    # library(DT)
    library(data.table)
    library(magrittr)

    library(plyr)
    library(purrr)
})

# load('data/phenoflux115_ET&GPP&VI.rda')
# load('inst/shiny/check_season/data/phenoflux_115.rda')
# load('inst/shiny/check_season/data/ET&GPP&VI_flux115.rda')
# sites <- sort(sites)

#' Generate DT::datatable
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

#' tidy GPP data
tidy_fluxGPP <- function(){
    gpp <- rowMeans(df[, .(GPP_DT, GPP_NT)], na.rm = T)
    gpp[gpp < 0] <- 0

    df <<- data.table(site = df$site, t = df$date, y = gpp) # , w = 1 (optional)
    # df
}

#' check_file
#' Check file whether exist. If not, then give a notification.
check_file <- function(file, duration = 10){
    filename <- deparse(substitute(file))
    if (is.null(file)) file <- "NULL"

    if (file.exists(file)) {
        TRUE
    } else {
        showNotification(sprintf("invalid %s: %s", filename, file),
                         duration = duration, type = "warning")
        FALSE
    }
}

#' Make sure date character in \code{df} has been converted to \code{Date} object.
check_datestr <- function(df){
    var_times <-  intersect(c("t", "date"), colnames(df))
    for (i in seq_along(var_times)){
        varname <- var_times[i]
        df[[varname]] %<>% lubridate::ymd()
    }
    df <<- df
}

#' update INPUT of all data
updateINPUT <- function(input){
    status <- FALSE
    if (input$file_type == '.rda | .RData') {
        file_rda  <- input$file_rda$datapath
        if (check_file(file_rda)) {
            load(file_rda)
            check_datestr(df)
            status <- TRUE
        }
    } else if (input$file_type == 'text'){
        file_veg  <- input$file_veg$datapath
        file_site <- input$file_site$datapath

        if (check_file(file_veg)) {
            df    <<- fread(file_veg)
            check_datestr(df)
            sites <<- unique(df$site) %>% sort()

            if (check_file(file_site)){
                st <<- data.table(ID = seq_along(sites), site = sites, lat = 30)
            } else {
                st <<- fread(file_site)
            }
            status <- TRUE
        }
    }
    # list(df = df, st = st, sites = sites)
    return(status)
}

#' update vegetation index variable Y in df
updateY <- function(input){
    if (!is.null(input$txt_varVI) && input$txt_varVI != "") {
        eval(parse(text = sprintf('df[, y := %s]', input$txt_varVI)))
    }
}

#' convert_QC2weight
convert_QC2weight <- function(input){
    qcFUN <- input$qcFUN
    varQC <- input$txt_varQC

    if (varQC %in% colnames(df)){
        warning(sprintf("No QC variable %s in df! ", varQC))
    }

    if (input$check_QC2weight && varQC %in% colnames(df)){
        eval(parse(text = sprintf('df[, c("w", "QC_flag") := %s(%s, wmin = 0.2)]',
            qcFUN, varQC)))
    }
}

################################################################################

getDf.site  <- function(df, sitename){
    dplyr::select(df[site == sitename, ], dplyr::matches("t|y|w"))
    #%T>% plotdata(365)
}

getINPUT.site <- function(df, st, sitename){
    sp       <- st[site == sitename]
    south    <- sp$lat < 0
    titlestr <- with(sp, sprintf("[%3d] %s, %s, lat = %.2f", ID, site, IGBP, lat))

    d     <- df[site == sitename, ]#%T>% plotdata(365)
    d_new <- add_HeadTail(d, south = south, nptperyear)
    INPUT <- do.call(check_input, d_new)

    INPUT$south    <- south
    INPUT$titlestr <- titlestr
    # list(INPUT = INPUT, plotdat = d)
    INPUT
}

#' Cal growing season dividing information
#'
#' @param input Shiny \code{input} variable
#' @param INPUT An object returned by \code{check_season}
#'
cal_season <- function(input, INPUT){
    param <- list(
        FUN_season     = input$FUN_season,
        rFUN           = input$rFUN,
        iters          = input$iters,
        lambda         = input$lambda,
        nf             = input$nf,
        frame          = input$frame,
        wFUN           = input$wFUN,
        maxExtendMonth = input$maxExtendMonth,
        rytrough_max   = input$rytrough_max,
        threshold_max  = input$threshold_max,
        threshold_min  = input$threshold_min
    )
    # param <- lapply(varnames, function(var) input[[var]])
    param <- c(list(INPUT = INPUT), param)
    # print(str(param))
    do.call(check_season, param) # brk return
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

phenofit_all <- function(input, progress = NULL){
    n   <- length(sites)
    res <- list()

    # parameters for Fine Fitting
    params_fineFitting <- list(
        methods      = input$FUN, #c("AG", "zhang", "beck", "elmore", 'Gu'), #,"klos",
        # debug        = FALSE,
        wFUN         = get(input$wFUN2),
        nextent      = 2,
        maxExtendMonth = 3,
        minExtendMonth = 1,
        qc           = 1,
        minPercValid = 0.2,
        print        = TRUE
    )

    showProgress <- !is.null(progress)
    if (showProgress){
        on.exit(progress$close())
        progress$set(message = sprintf("phenofit (n=%d) | running ", n), value = 0)
    }

    # print('debug 1 ...')
    # browser()

    for (i in 1:n){
        # tryCatch({
        # }, error = function(e){
        # })
        if (showProgress){
            progress$set(i, detail = paste("Doing part", i))
        }
        fprintf("phenofit (n = %d) | running %03d ... \n", i, n)

        sitename <- sites[i]
        INPUT    <- getINPUT.site(df, st, sitename)

        # Rough Fitting and gs dividing
        brks   <- cal_season(input, INPUT)

        params <- c(list(INPUT = INPUT, brks = brks), params_fineFitting)
        fit    <- do.call(curvefits, params)

        # Good of fitting of Fine Fitting
        stat <- ldply(fit$fits, function(fits_meth){
            ldply(fits_meth, statistic.phenofit, .id = "flag")
        }, .id = "meth")

        # Phenological Metrics
        p     <- lapply(fit$fits, ExtractPheno)
        pheno <- map(p, tidyFitPheno, origin = INPUT$t[1]) %>% purrr::transpose()

        fit   <- c(fit, list(INPUT = INPUT, seasons = brks, stat = stat, pheno = pheno))
        ############################# CALCULATION FINISHED #####################
        res[[i]] <- fit
    }
    set_names(res, sites)
}

# plot_data <- function(d, title){
#     par(setting)
#     do.call(check_input, d) %>% plotdata()
#     mtext(title, side = 2, line = 2, cex = 1.3, font = 2)
# }

################################################################################
# https://stackoverflow.com/questions/48592842/show-inf-in-dtdatatable
options(
    htmlwidgets.TOJSON_ARGS = list(na = 'string'),
    shiny.maxRequestSize=30*1024^2
)
