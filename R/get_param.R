#' Get parameters from curve fitting result
#'
#' @param x `fFITs` object returned by [curvefit()], or list of `fFITs` objects
#'
#' @return 
#' A list of `tibble` with the length being equal to the number of methods.
#' Each line of `tibble` cotains the corresponding parameters of each growing 
#' season.
#' 
#' @example inst/examples/ex-get_fitting_param_GOF.R
#' @export
get_param <- function(x) UseMethod("get_param", x)

#' @rdname get_param
#' @export
get_param.list <- function(x) {
    lapply(x, get_param.fFITs) %>%
        purrr::transpose() %>%
        map(~ map_df(., ~.x, .id = "flag")) # could improve
}

#' @rdname get_param
#' @export
get_param.fFITs <- function(x){
    map(x$model, ~.x[["par"]] %>% as_tibble())
}

#' @rdname get_param
#' @export
get_param.fFIT <- function(x) {
    x[["par"]] %>% as_tibble()
}
