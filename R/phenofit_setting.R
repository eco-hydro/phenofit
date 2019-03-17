#' @name setting
#' @title setting
NULL

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

        # growing season dividing parameters
        wFUN_rough             = input$wFUN_rough, 
        iters_rough            = input$iters_rough, 
        max_extend_month_rough = input$max_extend_month_rough, 
        calendarYear           = input$calendarYear, 
        FUN_season             = input$FUN_season, 
        FUN_rough              = input$FUN_rough, 
        r_max                  = input$r_max, 
        r_min                  = input$r_min, 
        rtrough_min            = input$rtrough_min, 

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
setting.read <- function(file){
    read_json(file)
}

#' @export
#' @rdname setting
setting.write <- function(pars, file){
    write_json(setting, "phenofit_setting.json", pretty = TRUE)
}
