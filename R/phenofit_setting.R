#' @name setting
#' @title setting
NULL

# phenofit parameters
options.phenofit <- list(
    file_site              = "", 
    file_veg_text          = "",
    file_veg_rda           = "",
    nptperyear             = NA_real_,

    var_y                  = "y", 
    var_qc                 = "",  # SummaryQA
    qcFUN                  = "",  # qc_summary

    # growing season dividing parameters
    calendarYear           = FALSE, 
    FUN_season             = "season_mov", 
    FUN_rough              = "wWHIT", 
    frame                  = NA_real_, 
    lambda                 = NA_real_,
    nf                     = 3,
    wFUN_rough             = "wKong", 
    iters_rough            = 4, 
    max_extend_month_rough = 2, 
    r_max                  = 0.2, 
    r_min                  = 0.0, 
    rtrough_max            = 0.6, 

    # fine fitting parameters
    FUN_fine               = "Beck", 
    wFUN_fine              = "wKong", 
    iters_fine             = 2, 
    max_extend_month_fine  = 2, 
    nextend_fine           = 2, 
    use.rough              = FALSE,

    # phenology extraction
    meths_pheno            = ""   
)


#' get parameters of phenofit shinyapp
#' 
#' @importFrom jsonlite write_json read_json
#' @export
#' @rdname setting
#' 
#' @examples
#' \dontrun{
#' pars = get_setting(input)
#' print(str(pars))
#' write_json(setting, "phenofit_setting.json", pretty = TRUE)
#' }
setting.get <- function(input){
    list(
        file_site              = input$file_site,
        file_veg_text          = input$file_veg_text,
        file_veg_rda           = input$file_veg_rda,
        nptperyear             = input$nptperyear,

        var_y                  = input$var_y, 
        var_qc                 = input$var_qc,  # SummaryQA
        qcFUN                  = input$qcFUN,   # qc_summary

        # growing season dividing parameters
        calendarYear           = input$calendarYear, 
        FUN_season             = input$FUN_season, 
        FUN_rough              = input$FUN_rough, 
        frame                  = input$frame, 
        lambda                 = input$lambda,
        nf                     = input$nf,
        
        wFUN_rough             = input$wFUN_rough, 
        iters_rough            = input$iters_rough, 
        max_extend_month_rough = input$max_extend_month_rough, 
        r_max                  = input$r_max, 
        r_min                  = input$r_min, 
        rtrough_max            = input$rtrough_max, 

        # fine fitting parameters
        FUN_fine               = input$FUN_fine, 
        wFUN_fine              = input$wFUN_fine, 
        iters_fine             = input$iters_fine, 
        max_extend_month_fine  = input$max_extend_month_fine, 
        nextend_fine           = input$nextend_fine, 
        use.rough              = input$use.rough,

        # phenology extraction
        meths_pheno            = input$meths_pheno
    )
    # TRS_sos
    # TRS_eos
}

#' @export
#' @rdname setting
setting.read <- function(file = "phenofit_setting.json"){
    read_json(file)
}

#' @export
#' @rdname setting
setting.write <- function(options, file = "phenofit_setting.json"){
    write_json(options, file, pretty = TRUE)
}

#' @export
#' @rdname setting
setting.update <- function(options){
    I <- match(names(options), names(options.phenofit))
    options.phenofit[I] <- options
    return(options.phenofit)
}
