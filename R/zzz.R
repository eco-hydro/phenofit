#' @title phenofit
#' @name phenofit
#' @aliases phenofit-package
#' @docType package
#' @keywords Vegetation phenology package
#' @description Vegetation phenology package
#' @import magrittr numDeriv plyr
#' @import tibble ggplot2 
#' @importFrom gridExtra arrangeGrob
#' @importFrom data.table data.table as.data.table := is.data.table fwrite fread
#' @importFrom zoo na.approx index zoo
#' @importFrom dplyr bind_cols bind_rows group_by first last nth
#' @importFrom purrr map map_df map_dbl is_empty
#' @importFrom tidyr gather spread
#' @importFrom lubridate ymd yday year month day dyears is.Date
#' @importFrom stringr str_extract
#' @importFrom utils object.size
#' @importFrom grDevices dev.off cairo_pdf
#' @import stats graphics
#' 
#' @useDynLib phenofit, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

.onLoad <- function (libname, pkgname){
    if(getRversion() >= "2.15.1") {
        utils::globalVariables(
            c(".SD", ".N", 
              "meth", "doy", "origin", # tidyFitPheno
              "DayOfYear", "SummaryQA", "site", "EVI", "w", "QC_flag", # tidy_MOD13.gee
              "beg", "end",  # plot_phenofit
              "val", "type", "flag", "peak" # season
            )
        )
    }
}
# .onUnload <- function (libpath) {
#   library.dynam.unload("phenofit", libpath)
# }
