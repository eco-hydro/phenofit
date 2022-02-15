fprintf <- function(fmt, ...) cat(sprintf(fmt, ...))

runningId <- function(i, step = 1, N, prefix = "") {
    perc <- ifelse(missing(N), "", sprintf(", %.1f%% ", i/N*100))
    if (mod(i, step) == 0) cat(sprintf("[%s] running%s %d ...\n", prefix, perc, i))
}

is_empty <- function(x) {
    is.null(x) || (is.data.frame(x) && nrow(x) == 0) || length(x) == 0
}

rm_empty <- function(x){
    if (is.list(x)){
        x[sapply(x, length) > 0]
    }else {
        x[!is.na(x)]
    }
}

last2 <- function(x) {
    if (is.list(x)) last(x) else x
}

guess_names <- function(x) {
    if (is.null(names(x))) names(x) = seq_along(x)
    x
}

which.notna <- function(x) which(!is.na(x))

match2 <- function (x, y) {
    if (is.null(x) || is.null(y)) return(NULL)

    I <- match(x, y)
    I_x <- which.notna(I)
    I_y <- I[I_x]
    d <- data.table(x = x[I_x], y = y[I_y], I_x, I_y, grp = cumsum(c(TRUE, 
        diff(I_y) != 1)))
    d
}

contain <- function(d, pattern = "NDVI|EVI") {
    names(d) %>% .[grep(pattern, .)]
}

listk <- function(...) {
    # get variable names from input expressions
    vars <- as.list(substitute(list(...)))[-1] # the first element is `list`
    names <- names(vars)

    Id_noname <- if (is.null(names)) seq_along(vars) else which(names == "")
    if (length(Id_noname) > 0)
        names[Id_noname] <- sapply(vars[Id_noname], deparse)
    setNames(list(...), names)
}

export2glob <- function(...) {
    list2env(listk(...), envir = .GlobalEnv)
    invisible()
}

#' melt_list
#' @keywords internal
#' @export
melt_list <- function (list, var.name = "variable", na.rm = TRUE, ...)
{
    if (is.null(names(list)))
        names(list) <- seq_along(list)
    list <- rm_empty(list)
    if (is.null(list) || length(list) == 0) {
        return(NULL)
    }
    first <- list[[1]]
    if (is.data.frame(first)) {
        names <- names(list)
        for (i in seq_along(list)) {
            x <- list[[i]]
            eval(parse(text = sprintf("x$%s <- names[i]",
                var.name)))
            list[[i]] <- x
        }
        res <- do.call(rbind, list) %>% data.table()
    }
    reorder_name(res, var.name)
}

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

# d <- data.table(a = 1:10, f1 = 1:10, f2 = 1:10, f3 = 1:10)
# select_vars(d, "f*")
# select_vars(d, "f")
select_vars <- function(x, pattern) {
    names <- names(x)
    # if (is.data.frame(x)) {
    #     names <- colnames(x)
    # }
    ind <- grep(pattern, names)
    # vars = names[ind]
    if (is.data.table(x)) {
        x[, .SD, .SDcols = ind]
    } else if (is.data.frame(x)) {
        x[, ind]
    }
}

#' @keywords internal
#' @export
zeallot::`%<-%`

#' @export
magrittr::`%>%`

#' @export
magrittr::`%<>%`

#' @importFrom lubridate make_date
#' @export
lubridate::make_date

clamp <- function(x, lims = c(0, 1), fill.na = FALSE) {
    if (fill.na) {
        x[x < lims[1]] <- NA_real_
        x[x > lims[2]] <- NA_real_
    } else {
        x[x < lims[1]] <- lims[1]
        x[x > lims[2]] <- lims[2]
    }
    x
}
