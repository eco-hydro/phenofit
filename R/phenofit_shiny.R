#' Starts the phenofit shiny interface
#' 
#' The GUI allows you to interactively visualize curve fitting time series and 
#' phenological metrics.
#' 
#' @keywords GUI, front end, interactive
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
    suppressWarnings(shiny::runApp(appDir,
        display.mode = "normal", launch.browser = TRUE))
}
