#' Get parameters from curve fitting result
#'
#' @param fits Multiple methods curve fitting results by `curvefits` result.
#' @inheritParams get_GOF
#'
#' @example inst/examples/ex-get_fitting_param_GOF.R
#' @export
get_param <- function(fits){
    llply(fits, get_param.fFITs) %>%
        purrr::transpose() %>%
        map(~ldply(., function(x) x, .id = "flag") %>% as_tibble())
}

#' @rdname get_param
#' @export
get_param.fFITs <- function(fFITs){
    map(fFITs$fFIT, "par")
}
