#' fprintf
#' Print sprintf result into console, just like C style fprintf function
#' @param fmt a character vector of format strings, each of up to 8192 bytes.
#' @param ... other parameters will be passed to \code{sprintf}
#' @export
fprintf <- function(fmt, ...) cat(sprintf(fmt, ...))

#' print the running ID in the console
#' 
#' @param i the running Id.
#' @param step how long of print step.
#' @param prefix prefix string
#' 
#' @rdname fprintf
#' @export
runningId <- function(i, step = 1, prefix = "") {
    if (mod(i, step) == 0) cat(sprintf("%s running %d ...\n", prefix, i))
}

#' retry
#' retry to execute expr until success before reaches the max try times (maxTimes)
#' 
#' @param expr expression to be evaluated
#' @param maxTimes maximum try times
#' 
#' @export
retry <- function(expr, maxTimes = 3){
    eTimes <- 0
    out    <- NULL
    while (eTimes < maxTimes){
      out <- tryCatch({
          expr
          eTimes <- maxTimes
      }, error = function(e){
          eTimes <<- eTimes + 1
          message(sprintf("[try%d]: %s", eTimes, e))
          NULL #If error return NULL
      })
    }
    return(out)
}

#' melt_list
#' 
#' @param list Data set to melt
#' @param var.name list names convert into var.name
#' @param na.rm Should NA values be removed from the data set? 
#' This will convert explicit missings to implicit missings.
#' @param ... other parameters to melt.
#' 
#' @rdname listk
#' @importFrom reshape2 melt
#' @export
melt_list <- function(list, var.name, na.rm = TRUE, ...){
    if (is.data.table(list[[1]])){
        names <- names(list)
        for (i in seq_along(list)){
            x <- list[[i]]
            eval(parse(text = sprintf("x$%s <- names[i]", var.name)))
            list[[i]] <- x
        }
        res <- do.call(rbind, list)#return
    } else{
        id.vars <- colnames(list[[1]])
        res <- reshape2::melt(list, ..., id.vars = id.vars, na.rm = na.rm)
        colnames(res) <- c(id.vars, var.name)
    }
    return(res)
}
# melt_list <- function(data, var.name, ..., na.rm = FALSE, value.name = "value") {
#   id.vars <- colnames(data[[1]])
#   res <- melt(data, ..., id.vars = id.vars, na.rm = na.rm, value.name = value.name)
#   colnames(res) <- c(id.vars, var.name)
#   return(res)
# }

#' listk
#' @export
listk <- function(...){
  # get variable names from input expressions
  cols <- as.list(substitute(list(...)))[-1]
  vars <- names(cols)
  Id_noname <- if (is.null(vars)) seq_along(cols) else which(vars == "")

  if (length(Id_noname) > 0)
    vars[Id_noname] <- sapply(cols[Id_noname], deparse)
  # ifelse(is.null(vars), Id_noname <- seq_along(cols), Id_noname <- which(vars == ""))
  x <- setNames(list(...), vars)
  return(x)
}

#' @rdname listk
#' @export
list.cbind <- function(x) do.call(cbind.data.frame, x) %>% set_colnames(names(x))

#' @rdname listk
#' @export
list.rbind <- function(x) do.call(rbind.data.frame, x) %>% set_rownames(names(x))#%>% set_rownames(names(x))

#' Add n-day flag
#' 
#' To aggregated data into n-day (e.g. 8-day, 16-day) like MODIS product, a 
#' n-day flag is need.
#' 
#' @param d data.frame or data.table
#' @param days Integer number or vector, can't have duplicated value.
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

#' reorder_name
#' @param headvars headvars will be in the head columns.
#' @param tailvars tailvars will be in the tail columns.
#' @rdname tools
#' @export
reorder_name <- function(d,
                         headvars = c("site", "date", "year", "doy", "d8", "d16"),
                         tailvars = ""){
    headvars %<>% intersect(colnames(d))
    tailvars %<>% intersect(colnames(d))
    varnames <- c(headvars, setdiff(colnames(d), union(headvars, tailvars)), tailvars)
    if (is.data.table(d)){
        # d[, ..varnames]
        d[, varnames, with = F] #return
    }else{
        d[, varnames]
    }
}

#' rm_empty
#' @param x A vector or list
#' 
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
#' @export
merge_pdf <- function(outfile = "RPlot.pdf", indir = 'Figure', pattern = "*.pdf", del = FALSE){
    # "Y:/R/phenofit/Figs/"
    files <- dir(indir, pattern, full.names = T)
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
obj.size <- function(obj, unit = "Mb"){
    cat(format(object.size(obj), unit), "\n")
}
