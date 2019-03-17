# source('main_phenofit.R')
suppressMessages({
    library(phenofit)
    library(shiny)
    # library(DT)
    library(data.table)
    library(magrittr)

    library(plyr)
    library(purrr)
    library(tidyselect)
})

# load('data/phenoflux115_ET&GPP&VI.rda')
# load('inst/shiny/check_season/data/phenoflux_115.rda')
# load('inst/shiny/check_season/data/ET&GPP&VI_flux115.rda')
# sites <- sort(sites)

################################################################################
# https://stackoverflow.com/questions/48592842/show-inf-in-dtdatatable
options(
    htmlwidgets.TOJSON_ARGS = list(na = 'string'),
    shiny.maxRequestSize=30*1024^2
)

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
    df
}

#' update all INPUT data according to \code{input} file.
updateINPUT <- function(input){
    status <- FALSE
    if (input$file_type == '.rda | .RData') {
        file_veg_rda  <- input$file_veg_rda$datapath
        if (check_file(file_veg_rda)) {
            load(file_veg_rda)
            check_datestr(df)
            status <- TRUE
        }
    } else if (input$file_type == 'text'){
        file_site <- input$file_site$datapath
        file_veg_text  <- input$file_veg_text$datapath

        if (check_file(file_veg_text)) {
            df    <<- fread(file_veg_text)
            check_datestr(df)
            sites <<- unique(df$site) %>% sort()

            if (check_file(file_site)){
                st <<- fread(file_site)
            } else {
                st <<- data.table(ID = seq_along(sites), site = sites, lat = 30)
            }
            status <- TRUE
        }
    }
    # list(df = df, st = st, sites = sites)
    return(status)
}

#' update vegetation index variable Y in df
#'
#' @param rv reactiveValues.
#' @param varname variable name of vegetation index.
update_VI <- function(rv, varname){
    # varname <- input$txt_varVI

    print('\t update_VI ...')
    if (!is.null(varname) && !(varname %in% c("", "y"))) {
        eval(parse(text = sprintf('rv$df$y <- rv$df$%s', varname)))
        df <<- rv$df
    }
}

#' convert_QC2weight
convert_QC2weight <- function(input){
    qcFUN <- input$qcFUN
    varQC <- input$txt_varQC

    if (!(varQC %in% colnames(df))){
        warning(sprintf("No QC variable %s in df! ", varQC))
    }

    if (input$check_QC2weight && varQC %in% colnames(df)){
        eval(parse(text = sprintf('df[, c("w", "QC_flag") := %s(%s, wmin = 0.2)]',
            qcFUN, varQC)))
        df <<- df
    }
}
