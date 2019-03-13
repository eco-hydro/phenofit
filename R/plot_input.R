#' Plot INPUT returned by check_input
#'
#' @inheritParams season
#' @param ... other parameter will be ignored.
#' 
#' @examples
#' library(phenofit)
#' data("MOD13A1")
#' 
#' dt <- tidy_MOD13.gee(MOD13A1$dt)
#' st <- MOD13A1$st
#' 
#' sitename <- dt$site[1]
#' d     <- dt[site == sitename, ] # get the first site data
#' sp    <- st[site == sitename, ] # station point
#' # global parameter
#' IsPlot = TRUE
#' print  = FALSE
#' nptperyear = 23
#' ypeak_min  = 0.05
#' 
#' dnew     <- add_HeadTail(d, nptperyear = nptperyear) # add one year in head and tail
#' INPUT    <- check_input(dnew$t, dnew$y, dnew$w, d$QC_flag, nptperyear, 
#'                         maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
#' plot_input(INPUT)
#' @export
plot_input <- function(INPUT, wmin = 0.2, ...){
    y <- INPUT$y
    t <- INPUT$t
    w <- INPUT$w
    n <- length(y)

    if (is.null(t)) {
        t <- 1:n
        years <- NULL
    } else {
        years  <- year(t)
    }
    if (is.null(w)) w <- rep(1, n)

    QC_flag <- INPUT$QC_flag
    if (is.null(QC_flag)){
        # divide into three categories: good, marginal and bad
        wf <- 4 - findInterval(w, c(-Inf, wmin, 0.5, 1), left.open = TRUE)
        QC_flag <- factor(wf, 1:3, c("good", "marginal", "cloud")) %>% as.character() %>% 
            factor(qc_levels)
    }

    # nptperyear <- INPUT$nptperyear
    npt <- length(y)
    # show grid lines
    par(mgp = c(1.5, 0.5, 0)) #oma = c(1, 2, 3, 1)
    # at <- t[seq(1, npt, nptperyear)]
    # fmt <- ifelse(yday(at[1]) == 1, "%Y", "%Y/%m/%d")
    # axis(side=1, at = at, labels = format(at, fmt))

    # colors <- c("grey60", "#00BFC4", "#C77CFF", "#F8766D", "blue", "red", "black")
    pch    <- c(19, 15, 17) # 4

    main <- 'Vegetation Index'
    if (!is.null(years) && length(unique(years)) < 3){
        plot(t, y, type = "l", xaxt="n", ann = FALSE, main = main, ...)
        axis.Date(1, at=seq(min(t), max(t), by="month"), format="%Y-%m")
    } else {
        plot(t, y, type = "l", ann = FALSE, main = main, ...)
    }

    # Ids <- unique(wf)
    for (i in 1:length(qc_levels)){
        I = QC_flag == qc_levels[i]
        add <- ifelse(i == 1, F, TRUE)
        points(t[I], y[I], pch = qc_shapes[i], bg = qc_colors[i], col = qc_colors[i], cex = 0.8)
    }
    
    # ylab = expression(paste('GPP ( gC ', m^-2, d^-1, ')'))
    if (!is.null(years)){
        date_beg <- ymd( min(years) *1e4 + 0101 )
        date_end <- ymd( max(years) *1e4 + 0101 )

        t_grids  <- seq.Date(date_beg, date_end, by = "year")
        abline(v = t_grids, col = "grey60", lty = 3)
    }
    grid(nx = NA, NULL)
    ylu <- INPUT$ylu
    if (!is.null(ylu)) abline(h=ylu, col="red", lty=2) # show ylims
}

# still under debug, 20180920
# plot_input.plotly <- function(){
#     w <- INPUT$w
#     if (is.null(w)) w <- rep(1, length(t))
#     # nptperyear <- INPUT$nptperyear

#     npt <- length(y)
#     # show grid lines
#     par(mgp = c(1.5, 0.5, 0)) #oma = c(1, 2, 3, 1)
#     # at <- t[seq(1, npt, nptperyear)]
#     # fmt <- ifelse(yday(at[1]) == 1, "%Y", "%Y/%m/%d")
#     # axis(side=1, at = at, labels = format(at, fmt))

#     # divide into three categories
#     wf <- 4 - findInterval(w, c(-Inf, wmin, 0.5, 1), left.open = TRUE)

#     colors <- c("grey60", "#00BFC4", "#F8766D", "#C77CFF", "blue", "red", "black")
#     pch    <- c(19, 15, 4)

#     # library(plotly)
#     scales <- ((0:2)*463) %>% paste0("m")
#     colors    <- c("red", "blue", "green") %>% set_names(scales)
#     colors <- c("grey60", "#00BFC4", "#F8766D")
#     shapes <- c(19, 15, 4)

#     d <- lst_sm$MOD13A1[site == sitename & scale == "0m"]
#     d[is.na(w), w := 0.1]
#     d$wf <- 4 - findInterval(d$w, c(-Inf, wmin, 0.5, 1), left.open = TRUE)
#     d$wf %<>% as.factor()

#     plot_ly(data = d, x = ~t, y = ~NDVI) %>%
#         #
#         add_trace(color = ~wf, symbol = ~wf,
#             colors = colors, symbols  = shapes,
#             type = 'scatter', mode = "markers") %>%
#         add_trace(name = "y", type = 'scatter', mode = "lines",
#                   line = list(width = 1, color = "black")) %>%
#         add_segments(x = -Inf, xend = Inf, y = 1, yend = 1)

#         #%>% #, shape = NULL, inherit = F) %>%
#     p <- ggplot(d, aes(t, NDVI)) +
#         geom_line() +
#         geom_point(aes(shape = wf, color = wf)) +
#         scale_shape_manual(values = shapes) +
#         scale_color_manual(values = colors)

# }