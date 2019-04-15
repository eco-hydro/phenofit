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

#' select_var_VI
select_var_VI <- function(df){
    tidyselect::vars_select(colnames(df),
        -tidyselect::matches("id|ID|site|scale|date|t|year|doy|DayOfYear|QA|QC|qa|qc|w"))
}

#' updateInput_phenofit
#'
#' Update shiny components when input changes.
#' @param session shiny session
#' @param rv reactiveValues defined in the server, has variables of
#' 'df', 'st' and 'sites'
updateInput_phenofit <- function(session, rv, init = FALSE) {
    ## update `nptperyear`
    colnames <- colnames(rv$df)
    sites <- rv$sites

    var_time   <-  intersect(c("date", "t"), colnames)[1]
    deltaT     <-  isolate( as.numeric(diff(rv$df[[var_time]][c(1, 2)])) )
    nptperyear <<- ceiling(365/deltaT)

    # browser()
    # update_VI(rv, input$txt_varVI) # update_VI right now.

    ## update `var_QC`
    varQCs <- colnames %>% .[grep("QA|QC|qa|qc", .)]
    if (length(varQCs) == 0){
        sel_qc_title <- paste0("vairable of QC: ",
            "No QC variables!")
        seq_qc <- ""; varQCs <- ""
    } else {
        sel_qc_title <- "vairable of QC:"
        sel_qc <- varQCs[1]
    }

    # updateNumericInput(session, 'nptperyear', value = nptperyear)
    # Do not update values when varQC changes. Because, it also rely on qcFUN
    updateSelectInput(session, "var_qc", sel_qc_title,
        choices = varQCs, selected = varQCs[1])
    updateSelectInput(session, "site",
                      choices = sites, selected = sites[1])
    NULL
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

    if (input$is_QC2w && varQC %in% colnames(rv$df)){
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
