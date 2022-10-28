env <- list2env(list(.julia = FALSE))

#' julia_init
#' 
#' @keywords internal
#' @import JuliaCall
#' @export
julia_setup <- function() {
    JuliaCall::julia_setup()
    # print(infile)
    JuliaCall::julia_library("phenofit")
    # JuliaCall::julia_library("nlminb")
    # infile <- system.file("julia/wBisquare.jl", package = "phenofit")
    # JuliaCall::julia_source(infile)
}

#' @importFrom JuliaCall julia_setup julia_source julia_call
julia_init <- function() {
    if (!env$.julia) {
        env$.julia = TRUE
        julia_setup()
    }
}


wBisquare_julia <- function(y, yfit, w, ..., wmin = 0.2, 
    trs_high = 0.6,
    trs_low  = trs_high,
    trs_bg = 0.2, 
    .toUpper = TRUE)
{
    julia_init()
    if (missing(w)) w  <- rep(1, length(y))
    wnew = JuliaCall::julia_call("phenofit.wBisquare", y, yfit, w,
        wmin = wmin,
        trs_high = trs_high,
        trs_low  = trs_low,
        trs_bg = trs_bg)
    return(wnew)
}
