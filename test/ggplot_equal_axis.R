# x, y, facet all quote variable
#'
#' @param d data.table or data.frame
#' @param x expression
#' @param y expression
#' @param facet quote variable, e.g. `quote(.(gof, index))`. facet can't put in
#' other variables, and then pass to this function.
#'
get_facetlims <- function(d, x, y, facet){
    x     <- substitute(x)
    y     <- substitute(y)
    facet <- substitute(facet)
    # facet <- substitute(facet, list(facet = facet))

    if (!is.data.table(d)) d <- data.table(d)

    # substitute(quote(list("xname" = xname)), list(xname = quote(GPP_avg)))
    # substitute(list(xname = xname), list2env(list(xname = quote(GPP_avg)), envir = .GlobalEnv))

    comd <- substitute(
        rbind(d[, .( xname = max(xname, yname), yname = max(xname, yname)), by],
              d[, .( xname = min(xname, yname), yname = min(xname, yname)), by])
        , listk(xname, yname, by = facet) )

    facetlims <- eval(comd)

    ncol <- ncol(facetlims)
    colnames(facetlims)[(ncol-1):ncol] <- c(deparse(xname), deparse(yname))
    facetlims
}

# need also update major_source and minor_source
equal_range <- function(p){
    b <- ggplot_build(p)

    b$layout$panel_ranges %<>% map(function(panel_range){
        x.range <- panel_range$x.range
        y.range <- panel_range$y.range

        min     <- pmin(x.range[1], y.range[1])
        max     <- pmax(x.range[1], y.range[1])
        range   <- c(min, max)
        panel_range$x.range <- range
        panel_range$y.range <- range

        return(panel_range)#fixed range
    })
    b$plot
}

