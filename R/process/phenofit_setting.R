#' @name setting
#' @title setting
NULL

# phenofit parameters
options.phenofit <- list(
    file_site              = "", 
    file_veg_text          = "",
    file_veg_rda           = "",
    file_type              = "RData",
    nptperyear             = NA_real_,
    ymin                   = 0, 
    var_y                  = "y", 
    is_QC2w                = FALSE,
    var_qc                 = "",  # SummaryQA
    qcFUN                  = "",  # qc_summary
    site                   = "",
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
#' @param options options of phenofit needed to export.
#' @param others options in others will override `options`. 
#' @param ... ignored.
#' 
#' @importFrom jsonlite write_json read_json
#' @export
#' @rdname setting
#' 
#' @examples
#' \dontrun{
#' pars = get_setting(options)
#' print(str(pars))
#' write_json(setting, "phenofit_setting.json", pretty = TRUE)
#' }
setting.get <- function(options, others = NULL, ...){
    options <- list(
        file_site              = options$file_site,
        file_veg_text          = options$file_veg_text,
        file_veg_rda           = options$file_veg_rda,
        file_type              = options$file_type,
        nptperyear             = options$nptperyear,
        ymin                   = options$ymin,

        var_y                  = options$var_y, 
        is_QC2w                = options$is_QC2w,
        var_qc                 = options$var_qc,  # SummaryQA
        qcFUN                  = options$qcFUN,   # qc_summary
        site                   = options$site,

        # growing season dividing parameters
        calendarYear           = options$calendarYear, 
        FUN_season             = options$FUN_season, 
        FUN_rough              = options$FUN_rough, 
        frame                  = options$frame, 
        lambda                 = options$lambda,
        nf                     = options$nf,
        
        wFUN_rough             = options$wFUN_rough, 
        iters_rough            = options$iters_rough, 
        max_extend_month_rough = options$max_extend_month_rough, 
        r_max                  = options$r_max, 
        r_min                  = options$r_min, 
        rtrough_max            = options$rtrough_max, 

        # fine fitting parameters
        FUN_fine               = options$FUN_fine, 
        wFUN_fine              = options$wFUN_fine, 
        iters_fine             = options$iters_fine, 
        max_extend_month_fine  = options$max_extend_month_fine, 
        nextend_fine           = options$nextend_fine, 
        use.rough              = options$use.rough,

        # phenology extraction
        meths_pheno            = options$meths_pheno
    )

    # inorder to replace variables in shinyFiles
    names  <- names(others)
    I_raw  <- match(names(others), names(options))
    I_nona <- which(!is.na(I_raw))

    options[I_raw[I_nona]] <- others[I_nona]
    return(options)
    # TRS_sos
    # TRS_eos
}

#' @param file file path of phenofit setting file (json).
#'  
#' @export
#' @rdname setting
setting.read <- function(file = "phenofit_setting.json"){
    if (file.exists(file)) {
        read_json(file) %>% map(unlist)
    } else {
        warning(sprintf('[w] setting file: %s not exist!', file))
    }
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
