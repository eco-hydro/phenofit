# stronger statDensity to show density prob
StatProb2d <- ggproto("StatProb2d", Stat,
    default_aes = aes(colour = "#3366FF", size = 0.5),
    required_aes = c("x", "y"),

    compute_group = function(data, scales, na.rm = FALSE, h = NULL,
                           contour = TRUE, n = 100, bins = NULL,
                           binwidth = NULL, breaks, ...) {
        if (is.null(h)) {
          h <- c(MASS::bandwidth.nrd(data$x), MASS::bandwidth.nrd(data$y))
        }

        dens <- MASS::kde2d(
          data$x, data$y, h = h, n = n,
          lims = c(scales$x$dimension(), scales$y$dimension())
        )

        dx <- diff(dens$x[1:2])  # lifted from emdbook::HPDregionplot()
        dy <- diff(dens$y[1:2])
        sz <- sort(dens$z)
        c1 <- cumsum(sz) * dx * dy

        prob <- approx(sz, 1-c1, dens$z)$y

        df <- data.frame(expand.grid(x = dens$x, y = dens$y), z = as.vector(dens$z))
        df$group <- data$group[1]
        df$z <- prob#replaced density with probability

        # print(str(data))
        # print(str(df))
        # print(scales)
        if (contour) {
            StatContour$compute_panel(df, scales, bins, binwidth, breaks = breaks)
            # temp$level %<>% as.factor()
            # print(str(temp))
            # print(unique(temp$level))
            # temp
        } else {
            names(df) <- c("x", "y", "density", "group")
            df$level <- 1
            df$piece <- 1
            df
        }
    }
)

stat_prob_2d <- function(mapping = NULL, data = NULL,
                            geom = "density_2d", position = "identity",
                            ...,
                            contour = TRUE,
                            n = 100,
                            h = NULL,
                            na.rm = FALSE,
                            show.legend = NA,
                            inherit.aes = TRUE) {
  layer(
    data        = data,
    mapping     = mapping,
    stat        = StatProb2d, #StatDensity2d, 
    geom        = geom,
    position    = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params      = list(na.rm = na.rm, contour = contour, n = n, h = h, ...)
  )
}




# compute_panel <- function (self, data, scales, ...) {
#     if (empty(data)) 
#         return(data.frame())
#     groups <- split(data, data$group)
#     stats <- lapply(groups, function(group) {
#         self$compute_group(data = group, scales = scales, ...)
#     })
#     stats <- mapply(function(new, old) {
#         if (empty(new)) 
#             return(data.frame())
#         unique <- uniquecols(old)
#         missing <- !(names(unique) %in% names(new))
#         cbind(new, unique[rep(1, nrow(new)), missing, drop = FALSE])
#     }, stats, groups, SIMPLIFY = FALSE)
#     do.call(plyr::rbind.fill, stats)
# }

# StatContour <- ggproto("StatContour", Stat,
#   required_aes = c("x", "y", "z"),
#   default_aes = aes(order = calc(level)),

#   compute_group = function(data, scales, bins = NULL, binwidth = NULL,
#                            breaks = NULL, complete = FALSE, na.rm = FALSE) {
#     # If no parameters set, use pretty bins
#     if (is.null(bins) && is.null(binwidth) && is.null(breaks)) {
#       breaks <- pretty(range(data$z), 10)
#     }
#     # If provided, use bins to calculate binwidth
#     if (!is.null(bins)) {
#       binwidth <- diff(range(data$z)) / bins
#     }
#     # If necessary, compute breaks from binwidth
#     if (is.null(breaks)) {
#       breaks <- fullseq(range(data$z), binwidth)
#     }

#     contour_lines(data, breaks, complete = complete)
#   }
# )

# geom_prob2d <- function(mapping = NULL, data = NULL,
#                             stat = "stat_prob2d", position = "identity",
#                             ...,
#                             lineend = "butt",
#                             linejoin = "round",
#                             linemitre = 10,
#                             na.rm = FALSE,
#                             show.legend = NA,
#                             inherit.aes = TRUE) 
# {
#   layer(
#     data = data,
#     mapping = mapping,
#     stat = stat,
#     geom = GeomDensity2d,
#     position = position,
#     show.legend = show.legend,
#     inherit.aes = inherit.aes,
#     params = list(
#       lineend = lineend,
#       linejoin = linejoin,
#       linemitre = linemitre,
#       na.rm = na.rm,
#       ...
#     )
#   )
# }

# GeomDensity2d <- ggproto("GeomDensity2d", GeomPath,
#   default_aes = aes(colour = "#3366FF", size = 0.5, linetype = 1, alpha = NA)
# )