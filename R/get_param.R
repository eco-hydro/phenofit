#' Get parameters from curve fitting result
#'
#' @param fits Multiple methods curve fitting results by `curvefits` result.
#' @inheritParams get_GOF
#'
#' @example inst/examples/ex-get_fitting_param_GOF.R
#' @export
get_param <- function(fits){
    lapply(fits, get_param.fFITs) %>%
        purrr::transpose() %>%
        map(~map_df(., ~ .x, .id = "flag")) # could improve
}

#' @rdname get_param
#' @export
get_param.fFITs <- function(fFITs){
    map(fFITs$model, ~.x[["par"]] %>% as_tibble())
}
