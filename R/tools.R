fprintf <- function(fmt, ...) cat(sprintf(fmt, ...))

runningId <- function(i, step = 1, N, prefix = "") {
    perc <- ifelse(missing(N), "", sprintf(", %.1f%% ", i/N*100))
    if (mod(i, step) == 0) cat(sprintf("[%s] running%s %d ...\n", prefix, perc, i))    
}

rm_empty <- function(x){
    if (is.list(x)){
        x[sapply(x, length) > 0]
    }else {
        x[!is.na(x)]
    }
}

contain <- function(d, pattern = "NDVI|EVI") {
    names(d) %>% .[grep(pattern, .)]
}

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

#' @importFrom data.table melt
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

first <- function(x, order_by = NULL, default = NA_real_) {
    nth(x, 1L, order_by = order_by, default = default)
}

last <- function(x, order_by = NULL, default = NA_real_) {
    nth(x, -1L, order_by = order_by, default = default)
}

nth <- function (x, n, order_by = NULL, default = NA_real_) {
    stopifnot(length(n) == 1, is.numeric(n))
    n <- trunc(n)
    if (n == 0 || n > length(x) || n < -length(x)) {
        return(default)
    }
    if (n < 0) {
        n <- length(x) + n + 1
    }
    if (is.null(order_by)) {
        x[[n]]
    }
    else {
        x[[order(order_by)[[n]]]]
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

mutate <- function (.data, ...) 
{
    stopifnot(is.data.frame(.data) || is.list(.data) || is.environment(.data))
    cols <- as.list(substitute(list(...))[-1])
    cols <- cols[names(cols) != ""]
    for (col in names(cols)) {
        .data[[col]] <- eval(cols[[col]], .data, parent.frame())
    }
    .data
}

check_function <- function(fun) {
    if (is.character(fun)) fun = get(fun)
    return(fun)
}

merge_pdf <- function(outfile = "RPlot.pdf", indir = "Figure",
         pattern = "*.pdf", del = FALSE) {
    files <- dir(indir, pattern, full.names = TRUE)
    order <- str_extract(basename(files), "(?<=\\[)\\d*(?=.*\\])") %>%
        as.numeric() %>% order()
    if (all(is.finite(order))) files <- files[order]
    # print(basename(files))
    pdftools::pdf_combine(files, outfile)
    if (del) file.remove(files)
}
