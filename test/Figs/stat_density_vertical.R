stat_densityV <- function(mapping = NULL, data = NULL,
                         geom = "area", position = "stack",
                         ...,
                         bw = "nrd0",
                         adjust = 1,
                         kernel = "gaussian",
                         n = 512,
                         trim = FALSE,
                         na.rm = FALSE,
                         show.legend = NA,
                         inherit.aes = TRUE) {

  layer(
    data = data,
    mapping = mapping,
    stat = StatDensityV,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      bw = bw,
      adjust = adjust,
      kernel = kernel,
      n = n,
      trim = trim,
      na.rm = na.rm,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
StatDensityV <- ggproto("StatDensityV", Stat,
  required_aes = "y",
  default_aes = aes(x = calc(density), fill = NA),

  compute_group = function(data, scales, bw = "nrd0", adjust = 1, kernel = "gaussian",
                           n = 512, trim = FALSE, na.rm = FALSE) {
    if (trim) {
      range <- range(data$y, na.rm = TRUE)
    } else {
      range <- scales$y$dimension()
    }

    compute_density(data$y, data$weight, from = range[1], to = range[2],
      bw = bw, adjust = adjust, kernel = kernel, n = n)
  }

)


compute_density <- function(x, w, from, to, bw = "nrd0", adjust = 1,
                            kernel = "gaussian", n = 512) {
  nx <- length(x)
  if (is.null(w)) {
    w <- rep(1 / nx, nx)
  }

  # if less than 2 points return data frame of NAs and a warning
  if (nx < 2) {
    warning("Groups with fewer than two data points have been dropped.", call. = FALSE)
    return(data.frame(
      x = NA_real_,
      density = NA_real_,
      scaled = NA_real_,
      count = NA_real_,
      n = NA_integer_
    ))
  }

  dens <- stats::density(x, weights = w, bw = bw, adjust = adjust,
    kernel = kernel, n = n, from = from, to = to)

  res <- data.frame(
    x = dens$x,
    density = dens$y,
    scaled =  dens$y / max(dens$y, na.rm = TRUE),
    count =   dens$y * nx,
    n = nx
  )
  res <- data.frame(
      x = dens$y,
      y = dens$x,
      density = dens$y,
      scaled =  dens$y / max(dens$y, na.rm = TRUE),
      count =   dens$y * nx,
      n = nx
  )
  res

}

# compute_density <- function(x, w, from, to, bw = "nrd0", adjust = 1,
#                             kernel = "gaussian", n = 512) {
#   nx <- length(x)
#   if (is.null(w)) {
#     w <- rep(1 / nx, nx)
#   }
#
#   # if less than 2 points return data frame of NAs and a warning
#   if (nx < 2) {
#     warning("Groups with fewer than two data points have been dropped.", call. = FALSE)
#     return(data.frame(
#       x = NA_real_,
#       density = NA_real_,
#       scaled = NA_real_,
#       count = NA_real_,
#       n = NA_integer_
#     ))
#   }
#
#   dens <- stats::density(x, weights = w, bw = bw, adjust = adjust,
#     kernel = kernel, n = n, from = from, to = to)
#
#   data.frame(
#     y = dens$x,
#     x = dens$y,
#     density = dens$x,
#     scaled =  dens$y / max(dens$y, na.rm = TRUE),
#     count =   dens$y * nx,
#     n = nx
#   )
# }


# ggplot(d, aes(y = GPP_avg)) +
#     geom_density(alpha = 0.6)
#     stat_densityV()
