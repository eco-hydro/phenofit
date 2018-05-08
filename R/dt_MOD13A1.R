#' dt_MOD13A1
#' 
#' A data.table dataset, raw data of MOD13A1 data, clipped in 10 representative points 
#' ('DE-Obe', 'IT-Col', 'CN-Cha', 'AT-Neu', 'ZA-Kru', 'AU-How', 'CA-NS6', 
#'  'US-KS2', 'CH-Oe2', 'CZ-wet').
#' 
#' Variables in dt_MOD13A1:
#' 
#' \describe{
#'   \item{site}{String, site name}
#'   \item{date}{Date, corresponding date}
#'   \item{year}{Numeric}
#'   \item{doy }{Numeric, Julian day of year}
#'   \item{d16 }{Numeric}
#'   \item{DetailedQA}{}
#'   \item{SummaryQA}{}
#'   \item{DayOfYear}{corresponding doy of compositing NDVI and EVI}
#'   \item{NDVI}{No buffered NDVI}
#'   \item{EVI}{No buffered EVI}
#'   \item{NDVI_500m}{500m buffered mean NDVI}
#'   \item{EVI_500m }{500m buffered mean EVI}
#'   \item{NDVI}{1000m buffered mean NDVI}
#'   \item{EVI }{1000m buffered mean EVI}
#'   \item{Tn  }{LST_Night_1km (degree) of MOD11A2}
#'   \item{w   }{Initial weights according to `SummaryQA` and \code{\link{qc_summary}}}
#' }
#' 
#' @docType data
#' @usage
#' data('dt_MOD13A1')
'dt_MOD13A1'