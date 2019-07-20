#' getRealDate
#'
#' convert MODIS `DayOfYear` to the exact compositing date.
#' @param date Date vector, the first day of the 16-day composite period.
#' @param DayOfYear Numeric vector, exact composite day of year. 
#' 
#' @keywords internal
#' @return
#' A data.table with a new column `t`, which is the exact compositing date.
#' 
#' @examples
#' library(phenofit)
#' data("MOD13A1")
#' 
#' df  <- MOD13A1$dt
#' df$t <- getRealDate(df$date, df$DayOfYear)
#' @export
getRealDate <- function(date, DayOfYear){
    year = year(date)
    doy  = as.integer(yday(date))

    I_na <- which(is.na(DayOfYear))
    DayOfYear[I_na] <- doy[I_na] # If DayOfYear is missing

    t  <- as.Date(sprintf("%d-%03d", year  , DayOfYear), "%Y-%j")

    # In case of last scene of a year, doy of last scene could in the next year
    t_nextYear <- as.Date(sprintf("%d-%03d", year+1, DayOfYear), "%Y-%j")
    I_nextYear <- abs(DayOfYear - doy) >= 300

    t[I_nextYear] <- t_nextYear[I_nextYear]
    return(t)
}

#' tidy_MOD13.gee
#'
#' Tidy MODIS 'MOD13' VI products' (e.g. MOD13A1, MOD13A2, ...) raw data exported from
#' Google Earth Engine.
#' Tidy contents include: \cr
#' 1. add exact compositing date, see [getRealDate()]. \cr
#' 2. Init weigths according `SummaryQA`, see [qc_summary()]. \cr
#' 
#' @inheritParams check_input
#' @param infile A character csv file path or a data.table
#' @param outfile Output file name. If missing, will not be written to file.
#'
#' @keywords internal
#' @return
#' A tidied data.table, with columns of 'site', 'y', 't', 'w', 'date' and
#' 'SummaryQA'.
#' * `site`: site name
#' * `y`: real value of EVI, `[-1, 1]`
#' * `date`: image date
#' * `t`: exact compositing date constructed from `DayOfYear`
#' * `w`: weights
#' * `SummaryQA`: A factor, QA types, one of "good", "margin", "snow/ice"
#' or "cloud".
#' 
#' @examples
#' library(phenofit)
#' data("MOD13A1")
#' dt <- tidy_MOD13.gee(MOD13A1$dt)
#' @export
tidy_MOD13.gee <- function(infile, outfile, wmin = 0.2){
    df <- infile
    if ("character" %in% class(infile) ) df <- fread(infile)

    df$date %<>% ymd()
    df$t <- getRealDate(df$date, df$DayOfYear)
    
    # Initial weights
    df[, c("QC_flag", "w") := qc_summary(SummaryQA, wmin = 0.2)]
    # Remap SummaryQA factor level, plot_phenofit use this variable. For other
    # remote sensing data without `SummaryQA`, need to modify `plot_phenofit`
    
    df <- df[, .(site, y = EVI/1e4, date, t, w, QC_flag
                 # IGBPcode,
                 # IGBPname = as.factor(IGBPnames[IGBPcode]),
                 )]
    if (!missing(outfile)) fwrite(df, outfile)
    df
    # merge coordinate info
    # df <- merge(df, st[, .(Id = site, lat, lon = long, IGBPname = IGBP)], by = "Id")
}
