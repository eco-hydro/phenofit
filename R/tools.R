#' fprintf
#' Print sprintf result into console, just like C style fprintf function
#' @param fmt a character vector of format strings, each of up to 8192 bytes.
#' @param ... other parameters will be passed to \code{sprintf}
#' 
#' @examples
#' cat(fprintf("%s\n", "Hello phenofit!"))
#' @export
fprintf <- function(fmt, ...) cat(sprintf(fmt, ...))

#' print the running ID in the console
#' 
#' @param i the running Id.
#' @param step how long of print step.
#' @param prefix prefix string
#' 
#' @examples
#' for (i in 1:10){
#'     runningId(i, prefix = "phenofit")
#' }
#' @rdname fprintf
#' @export
runningId <- function(i, step = 1, prefix = "") {
    if (mod(i, step) == 0) cat(sprintf("%s running %d ...\n", prefix, i))
}


#' Add n-day flag
#' 
#' To aggregated data into n-day (e.g. 8-day, 16-day) like MODIS product, a 
#' n-day flag is need.
#' 
#' @param d data.frame or data.table
#' @param days Integer number or vector, can't have duplicated value.
#' 
#' @examples
#' date = seq.Date(as.Date("2010-01-01"), as.Date("2010-12-31"), by = "day")
#' d <- data.frame(date)
#' dnew <- add_dn(d, days = c(8, 16))
#' @export
add_dn <- function(d, days = 8){
    if (class(d$date) != 'Date')
        d$date %<>% ymd()
    
    d %<>% plyr::mutate(d, year = year(date), doy = yday(date))
    
    days <- floor(days)
    for (i in seq_along(days)){
        day <- days[i]
        # d$d8  = ceiling(d$doy/8)
        eval(parse(text = sprintf("d$d%d <- ceiling(d$doy/%d)", day, day)))
    }
    return(d)
}

#' rm_empty
#' @param x A vector or list
#' 
#' @examples
#' # numeric
#' x <- c(1:5, NA)
#' rm_empty(x)
#' 
#' # list
#' l <- list(1:5, NULL, NA)
#' rm_empty(l)
#' @rdname tools
#' @export
rm_empty <- function(x){
    if (is.list(x)){
        x[sapply(x, length) > 0]
    }else {
        x[!is.na(x)]
    }
}

#' contain
#' find assigned pattern variable names
#' @param d A data.frame vector, or list
#' @param pattern string used to match \code{names(d)}
#' 
#' @examples
#' df <- data.frame(year = 2010, day = 1:3, month = 1, site = "A")
#' contain(df, "year|month|day")
#' @rdname tools
#' @export
contain <- function(d, pattern = "NDVI|EVI") {
    names(d) %>% .[grep(pattern, .)]
}

#' merge_pdf
#'
#' rely on python pdfmerge package, \code{pip install pdfmerge}
#' 
#' @param outfile String
#' @param indir Directory to search pdfs
#' @param pattern string used to match pdf filename
#' @param del Booolean. If true, after merge, original pdf files will be delete.
#' 
#' @examples
#' \dontrun{
#' merge_pdf("RPlot.pdf", indir = 'Figure', pattern = "*.pdf") 
#' }
#' @export
merge_pdf <- function(outfile = "RPlot.pdf", indir = 'Figure', pattern = "*.pdf", del = FALSE){
    # "Y:/R/phenofit/Figs/"
    files <- dir(indir, pattern, full.names = TRUE)
    cmd <- sprintf("pdfmerge -o %s %s", outfile, paste(files, collapse = " "))

    shell(cmd, wait = del)
    # shell(sprintf('pdfmerge -o %s %s', outfile, pattern) )
    if (del) file.remove(files)
}


#' obj.size
#' 
#' Get object size in \code{unit}
#' @param obj Object
#' @param unit "Kb", "Mb" or "Gb"
#' 
#' @examples
#' obj.size(1:100)
#' @export
obj.size <- function(obj, unit = "Mb"){
    cat(format(object.size(obj), unit), "\n")
}
