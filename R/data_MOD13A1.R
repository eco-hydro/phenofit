#' MOD13A1
#' 
#' A data.table dataset, raw data of MOD13A1 data, clipped in 10 representative points 
#' ('DE-Obe', 'IT-Col', 'CN-Cha', 'AT-Neu', 'ZA-Kru', 'AU-How', 'CA-NS6', 
#'  'US-KS2', 'CH-Oe2', 'CZ-wet').
#' 
#' Variables in MOD13A1:
#' 
#' * __`dt`__: vegetation index data
#'   * `system:index`: image index
#'   * `DayOfYear`: Numeric, Julian day of year
#'   * `DayOfYear`: corresponding doy of compositing NDVI and EVI
#'   * `DetailedQA`: VI quality indicators
#'   * `SummaryQA`: Quality reliability of VI pixel
#'   * `EVI`: Enhanced Vegetation Index
#'   * `NDVI`: Normalized Difference Vegetation Index
#'   * `date`: Date, corresponding date
#'   * `site`: String, site name
#'   * `sur_refl_b01`: Red surface reflectance
#'   * `sur_refl_b02`: NIR surface reflectance
#'   * `sur_refl_b03`: Blue surface reflectance
#'   * `sur_refl_b07`: MIR surface reflectance
#'   * `.geo`: geometry
#' 
#' * __`st`__: station info
#'   * `ID`: site ID
#'   * `site`: site name
#'   * `lat`: latitude
#'   * `lon`: longitude
#'   * `IGBPname`: IGBP land cover type
#' 
#' @docType data
#' @usage
#' data('MOD13A1')
#' @references
#' 1. https://code.earthengine.google.com/dataset/MODIS/006/MOD13A1
'MOD13A1'
