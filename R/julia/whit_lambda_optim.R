#' @rdname lambda_vcurve
#' @export
lambda_vcurve_jl <- function(y, w, lg_lambda_min = 0.1, lg_lambda_max = 3) {
  JuliaCall::julia_call("phenofit.lambda_vcurve", y, w,
    lg_lambda_min = lg_lambda_min,
    lg_lambda_max = lg_lambda_max
  )
}

#' @rdname lambda_vcurve
#' @export
lambda_cv_jl <- function(y, w, d = 2, lg_lambda_min = 0.1, lg_lambda_max = 3) {
  JuliaCall::julia_call("phenofit.lambda_cv", y, w,
    lg_lambda_min = lg_lambda_min,
    lg_lambda_max = lg_lambda_max
  )
}
