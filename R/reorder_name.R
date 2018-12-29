#' reorder_name
#'
#' reorder the name of data.frame, date.table or list.
#'
#' @param headvars headvars will be in the head columns.
#' @param tailvars tailvars will be in the tail columns.
#'
#' @examples
#' df <- data.frame(year = 2010, day = 1:3, month = 1, site = "A")
#' dt <- data.table::data.table(year = 2010, day = 1:3, month = 1, site = "A")
#' l <- list(year = 2010, day = 1:3, month = 1, site = "A")
#'
#' newname <- c("site", "year")
#' reorder_name(df, newname)
#' reorder_name(dt, newname)
#' reorder_name(l, newname)
#'
#' @rdname tools
#' @export
reorder_name <- function(
    d,
    headvars = c("site", "date", "year", "doy", "d8", "d16"),
    tailvars = "")
{
    names <- names(d)
    headvars %<>% intersect(names)
    tailvars %<>% intersect(names)
    varnames <- c(headvars, setdiff(names, union(headvars, tailvars)), tailvars)

    if (is.data.table(d)) {
        d[, varnames, with = F]
    } else if (is.data.frame(d)) {
        d[, varnames]
    } else if (is.list(d)){
        d[varnames]
    } else{
        stop("Unknown data type!")
    }
}
