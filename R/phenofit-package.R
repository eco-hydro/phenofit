#' @title phenofit
#' @name phenofit
#' @aliases phenofit-package
#' @docType package
#' @keywords Vegetation phenology package
#' @description Vegetation phenology package
#' @import magrittr numDeriv plyr
#' @import tibble ggplot2 
#' @import foreach iterators
#' @importFrom gridExtra arrangeGrob
#' @importFrom data.table data.table as.data.table := is.data.table fwrite fread
#' @importFrom zoo na.approx index zoo
#' @importFrom dplyr bind_cols bind_rows group_by first last nth
#' @importFrom purrr map map_df map_dbl is_empty
#' @importFrom tidyr gather spread
#' @importFrom lubridate ymd yday year month day dyears is.Date
#' @importFrom stringr str_extract
#' @importFrom utils object.size
#' @importFrom grDevices dev.off cairo_pdf colorRampPalette
#' @importFrom jsonlite read_json write_json
#' @importFrom shiny getDefaultReactiveDomain showNotification
#' @import stats graphics
#' 
#' @useDynLib phenofit, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

.onLoad <- function (libname, pkgname){
    if(getRversion() >= "2.15.1") {
        utils::globalVariables(
            c(".SD", ".N", 
              "meth", "doy", "origin", # tidyFitPheno
              "DayOfYear", "SummaryQA", "site", "EVI", "w", "QC_flag", # tidy_MOD13.gee
              "beg", "end",  # plot_phenofit
              "val", "type", "flag", "peak", # season, 
              "i", "qc", "y", "sitename" # phenofit_TS.avhrr
            )
        )
    }
}

# .onUnload <- function (libpath) {
#   library.dynam.unload("phenofit", libpath)
# }
