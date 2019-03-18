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

# convert shinypath to normal path
tidy_shiny_filepath <- function(path){
    if (is.list(path)) {
        path$datapath
    } else if (is.character(path)) {
        path
    } else {
        ""
    }
}

#' check_file
#' Check file whether exist. If not, then give a notification.
check_file <- function(file, duration = 10){
    filename <- deparse(substitute(file))

    if (length(file) == 0 || !is.character(file)) {
        return (FALSE)
    } else if (file.exists(file)) {
        return(TRUE)
    } else {
        if (duration != NULL){
            showNotification(sprintf("invalid %s: %s", filename, as.character(file)),
                duration = duration, type = "warning")
        }
        return(FALSE)
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
#'
#' @param input should has children of \code{file_site}, and one of
#' \code{file_veg_rda} or \code{file_veg_text}.
#'
#' assign following variables to parent:
#' df, st, sites
load_input <- function(input){
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
            df    <- fread(file_veg_text)
            check_datestr(df)
            sites <- unique(df$site) %>% sort()

            if (check_file(file_site)){
                st <- fread(file_site)
            } else {
                st <- data.table(ID = seq_along(sites), site = sites, lat = 30)
            }
            status <- TRUE
        }
    }

    list2env(listk(df, st, sites), env = parent.frame()) # to calling env
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
    }
}

#' convert_QC2weight
convert_QC2weight <- function(input, rv){
    qcFUN <- input$qcFUN
    varQC <- input$var_qc

    if (!(varQC %in% colnames(rv$df))){
        warning(sprintf("No QC variable %s in df! ", varQC))
    }

    if (input$check_QC2weight && varQC %in% colnames(rv$df)){
        eval(parse(text = sprintf('rv$df[, c("QC_flag", "w") := %s(%s, wmin = 0.2)]',
            qcFUN, varQC)))
    }
}


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

## initial parameters ----------------------------------------------------------
initparam_FUN_rough <- function(nptperyear, session) {
    frame <- nptperyear / 5 * 2 + 1

    if (nptperyear >= 300) {
        lambda <- 1e4
    } else if (nptperyear >= 40){
        lambda <- 15
    } else if (nptperyear >= 20){
        lambda <- 5
    } else {
        lambda <- 2
    }

    if (!missing(session)) {
        # Update wSG parameter
        updateNumericInput(session,
            "frame", "moving window size (frame):",
            floor(nptperyear / 5 * 2 + 1),
            floor(nptperyear / 12),
            floor(nptperyear / 2),
            floor(nptperyear / 12)
        )
        # Update Whittaker parameter
        updateNumericInput(session, "lambda", value = lambda)
    }
}

#' @examples
#' init_options(options, rv, input, session)
load_data <- function(options, ...){
    file_type <- options$file_type
    file_veg_rda  <- options$file_rda
    file_veg_text <- options$file_veg_text
    file_site <- options$file_site

    if (file_type == "text") {
        if (check_file(file_veg_text)) {
            df    <- fread(file_veg_text) %>% check_datestr()
            sites <- unique(df$site) %>% sort()

            if (check_file(file_site)) {
                st <- fread(file_site)
            } else {
                st <- data.table(ID = seq_along(sites), site = sites, lat = 30)
            }
        }
    } else {
        if (check_file(file_veg_rda)) {
            # options$file_veg_rda <- file_veg_rda
            load(file_veg_rda)
            df <- df %>% check_datestr()
            sites <- unique(df$site) %>% sort()
            # rv$st <- st
        }
    }
    
    varname <- options$var_y
    varQC <- options$var_qc
    qcFUN <- options$qcFUN

    if (!is.null(varname) && !(varname %in% c("", "y"))) {
        eval(parse(text = sprintf('df$y <- df$%s', varname)))
    }

    if (!(varQC %in% colnames(df))){
        warning(sprintf("No QC variable %s in df! ", varQC))
    }

    if (varQC %in% colnames(df)){
        eval(parse(text = sprintf('df[, c("QC_flag", "w") := %s(%s, wmin = 0.2)]',
            qcFUN, varQC)))
    }

    # isolate({
    #     rv$st <- st    
    #     rv$sites <- sites
    #     rv$df <- df
    # })
    nptperyear <- options$nptperyear

    listk(df, st, sites, nptperyear)
}
