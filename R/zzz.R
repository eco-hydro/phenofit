#' @title MissInfo
#' @name MissInfo
#' @aliases MissInfo-package
#' @docType package
#' @keywords missing information, interpolate
#' @description detect missing information of meteorological data and interpolate missing values using common interpolation methods
#' @import magrittr numDeriv plyr
#' @import tibble ggplot2
#' @importFrom zoo na.approx index zoo
#' @importFrom dplyr bind_cols bind_rows group_by first last nth
#' @importFrom purrr map map_df is_empty
#' @importFrom tidyr gather spread
#' @importFrom lubridate ymd
#' @importFrom stringr str_extract

#' @useDynLib phenofit
#' @importFrom Rcpp sourceCpp
NULL

# .onLoad <- function (libname, pkgname){
#   # print(search())
# }
# .onUnload <- function (libpath) {
#   library.dynam.unload("phenofit", libpath)
# }
