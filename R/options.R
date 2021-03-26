.options <- list(
    verbose_season_mov = TRUE,
    verbose_season     = FALSE,
    verbose_curvefit   = TRUE, 
    wFUN_rough         = "wTSM", 
    wFUN_fine          = "wTSM", 
    FUN_rough          = "wHANTS", 
    calendarYear       = TRUE,
    south              = FALSE
)

#' set and get phenofit option
#'
#' @param options
#' - style: font style, one of "EN" and "CH"
#' - family: font family name, e.g. "Times"
#' - family_CH: Chinese fontname for panel.barchartFreq yaxis title
#' 
#' @examples 
#' set_options(verbose_curvefit = FALSE)
#' get_options("verbose_season")
#' @export
set_options <- function(...) {
    opts = list(...)
    .options <<- modifyList(.options, opts)
}

#' @rdname set_options
#' @export
get_options <- function(names = NULL) {
    if (is.null(names)) return(.options)
    .options[[names]]
}
