# , .progress = "none", .inform = FALSE, .parallel = FALSE, .paropts = NULL
#' @importFrom purrr as_mapper
llply <- function(.data, .f = NULL, ...) {
    if (is_empty(.data)) return(.data)
    .f <- as_mapper(.f, ...)
    lapply(.data, .f, ...)
}

ldply <- function(.data, .f = NULL, ...) {
    llply(.data, .f, ...) %>% as.data.table()
}
