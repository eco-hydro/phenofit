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
#' dcast
#' @importFrom zoo na.approx index zoo
#' @importFrom purrr map is_empty %>%
#' @importFrom lubridate ymd yday year month day dyears is.Date
#' @importFrom utils object.size
#' @importFrom grDevices dev.off cairo_pdf colorRampPalette
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
            c(".SD", ".N", "..vars", 
              "left", "len", "right", "y_peak",
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
