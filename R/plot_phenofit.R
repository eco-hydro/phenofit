#' plot_phenofit
#'
#' @param obj `phenofit` object, list(INPUT, fit, seasons)
#' @param type one of c("season", "fitting", "pheno", "all")
#' @inheritParams curvefits
#' @inheritParams plot_curvefits
#' @param IsPlot boolean. If false, a ggplot object will be returned.
#' @param show.legend If now show legend, ggplot object will be returned, else
#' grid object will be returned.
#' @param newpage boolean, whether draw figure in a new page?
#' @param ... other parameters to [plot_curvefits()]
#' @keywords internel
#' @export
plot_phenofit <- function(obj, type = "all",
                          methods,
                          title = NULL, ylab = "Vegetation Index",
                          IsPlot = TRUE, show.legend = TRUE, newpage = TRUE, ...) {
    if (missing(methods) || is.null(methods)) {
        methods <- names(obj$fit[[1]]$fFIT)
    }
    g <- NULL

    plot_fitting <- function() {
        df_fit <- get_fitting(obj$fit)
        df_fit <- df_fit[meth %in% methods]

        g <- plot_curvefits(df_fit, obj$brks, title = title, ylab = ylab, show.legend = show.legend, ...)

        if (IsPlot) {
            if (newpage) grid::grid.newpage()
            grid::grid.draw(g)
        }
        return(g)
    }
    
    if (type == "fitting") {
        g <- plot_fitting()
    } else if (type == "season") {
        if (!is.null(obj$INPUT)) plot_season(obj$INPUT, obj$seasons)
    } else if (type == "pheno") {
        l_pheno <- get_pheno(obj$fit, methods, IsPlot = T)
    } else if (type == "all") {
        # fitting
        g <- plot_fitting()
        # season
        if (!is.null(obj$INPUT)) plot_season(obj$INPUT, obj$seasons)
        # pheno
        l_pheno <- get_pheno(obj$fit, methods, IsPlot = T)
    } else {
        stop(sprintf("[e] wrong type: %sn", type))
    }
    return(g)
}
