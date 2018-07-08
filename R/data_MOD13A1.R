#' MOD13A1
#' 
#' A data.table dataset, raw data of MOD13A1 data, clipped in 10 representative points 
#' ('DE-Obe', 'IT-Col', 'CN-Cha', 'AT-Neu', 'ZA-Kru', 'AU-How', 'CA-NS6', 
#'  'US-KS2', 'CH-Oe2', 'CZ-wet').
#' 
#' Variables in MOD13A1:
#' 
#' dt: vegetation index data
#' \describe{
#'   \item{system:index}{image index}
#'   \item{DayOfYear}{Numeric, Julian day of year}
#'   \item{DayOfYear}{corresponding doy of compositing NDVI and EVI}
#'   \item{DetailedQA}{VI quality indicators}
#'   \item{SummaryQA}{Quality reliability of VI pixel}
#'   \item{EVI}{Enhanced Vegetation Index}
#'   \item{NDVI}{Normalized Difference Vegetation Index}
#'   \item{date}{Date, corresponding date}
#'   \item{site}{String, site name}
#'   \item{sur_refl_b01}{Red surface reflectance}
#'   \item{sur_refl_b02}{NIR surface reflectance}
#'   \item{sur_refl_b03}{Blue surface reflectance}
#'   \item{sur_refl_b07}{MIR surface reflectance}
#'   \item{.geo}{geometry}
#' }
#' 
#' st: station info
#' \describe{
#'   \item{ID}{site ID}
#'   \item{site}{site name}
#'   \item{lat}{latitude}
#'   \item{lon}{longitude}
#'   \item{IGBPname}{IGBP land cover type}
#' }
#' 
#' @docType data
#' @usage
#' data('MOD13A1')
#' @references
#' https://code.earthengine.google.com/dataset/MODIS/006/MOD13A1
'MOD13A1'
