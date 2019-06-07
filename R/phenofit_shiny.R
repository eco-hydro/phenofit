#' @name phenofit_shiny
#' @title phenofit shiny app
#' 
#' @description The GUI allows you to interactively visualize curve fitting 
#' time series and phenological metrics.
#' 
#' @importFrom shiny runApp
#' @export
#' 
#' @examples
#' \dontrun{
#' phenofit_shiny
#' }
phenofit_shiny = function(){
    # if(!requireNamespace(c("DT",
    #                      "plotly",
    #                      "shinydashboard",
    #                      "leaflet"), quietly = TRUE)){
    #     stop("Packages \"DT, plotly, shinydashboard and leaflet\" are needed 
    #         for this function to work. Please install it.",
    #         call. = FALSE)
    # }
    appDir = sprintf("%s/shiny/phenofit", path.package("phenofit"))
    suppressWarnings(runApp(appDir,
        display.mode = "normal", launch.browser = TRUE))
}
