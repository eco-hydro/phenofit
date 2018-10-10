#' getRealDate
#'
#' convert MODIS \code{DayOfYear} to the exact compositing date.
#' @param df A data.table, at least with columns of \code{date}, and
#' \code{DayOfYear}.
#'
#' @return
#' A data.table with a new column \code{t}, which is the exact compositing date.
#' @export
getRealDate <- function(df){
    df[, `:=`(date = ymd(date))]
    df[, `:=`(year = year(date), doy = as.integer(yday(date)))]
    df[is.na(DayOfYear), DayOfYear := doy] # If DayOfYear is missing

    # In case of last scene of a year, doy of last scene could in the next year
    df[abs(DayOfYear - doy) >= 300, t := as.Date(sprintf("%d-%03d", year+1, DayOfYear), "%Y-%j")] # last scene
    df[abs(DayOfYear - doy) <  300, t := as.Date(sprintf("%d-%03d", year  , DayOfYear), "%Y-%j")]
    df
}

qc_values <- c("0", "1", "2", "3")
qc_levels <- c("good", "margin", "snow/ice", "cloud")
qc_colors <- c("grey60", "#00BFC4", "#F8766D", "#C77CFF") %>% set_names(qc_levels)
qc_shapes <- c(19, 15, 4, 17) %>% set_names(qc_levels)

#' tidy_MOD13.gee
#'
#' Tidy MODIS 'MOD13' VI products' (e.g. MOD13A1, MOD13A2, ...) raw data exported from
#' Google Earth Engine.
#' Tidy contents include: \cr
#' 1. add exact compositing date, see \code{\link{getRealDate}}. \cr
#' 2. Init weigths according \code{SummaryQA}, see \code{\link{qc_summary}}. \cr
#' 
#' @inheritParams check_input
#' @param infile A character csv file path or a data.table
#' @param outfile Output file name. If missing, will not be written to file.
#'
#' @return
#' A tidied data.table, with columns of 'site', 'y', 't', 'w', 'date' and
#' 'SummaryQA'.
#' \describe{
#' \item{site}{site name}
#' \item{y}{real value of EVI, [-1, 1]}
#' \item{date}{image date}
#' \item{t}{exact compositing date constructed from \code{DayOfYear}}
#' \item{w}{weights}
#' \item{SummaryQA}{A factor, QA types, one of "good", "margin", "snow/ice"
#' or "cloud".}
#' }
#' @export
tidy_MOD13.gee <- function(infile, outfile, wmin = 0.2){
    df <- infile
    if ("character" %in% class(infile) ) df <- fread(infile)
    df %<>% getRealDate()

    # Initial weights
    df[, c("w", "QC_flag") := qc_summary(SummaryQA, wmin = 0.2)]
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
