#' @title phenofit
#' @name phenofit
#' @aliases phenofit-package
#' @docType package
#' @keywords Vegetation phenology package
#' @description Vegetation phenology package
#' @import magrittr numDeriv plyr
#' @import tibble ggplot2 
#' @importFrom gridExtra arrangeGrob
#' @importFrom data.table data.table as.data.table := is.data.table
#' @importFrom zoo na.approx index zoo
#' @importFrom dplyr bind_cols bind_rows group_by first last nth
#' @importFrom purrr map map_df is_empty
#' @importFrom tidyr gather spread
#' @importFrom lubridate ymd yday month day dyears is.Date
#' @importFrom stringr str_extract
#' @importFrom utils object.size
#' @importFrom grDevices dev.off
#' @import stats graphics
#' 
#' @useDynLib phenofit
#' @importFrom Rcpp sourceCpp
NULL

# .onLoad <- function (libname, pkgname){
#   # print(search())
# }
# .onUnload <- function (libpath) {
#   library.dynam.unload("phenofit", libpath)
# }
