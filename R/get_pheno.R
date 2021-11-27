# ' PhenoPlot
# '
# ' @inheritParams check_input
# ' @param main figure title
# ' @param ... ignored parameters
# '
# ' @export
PhenoPlot <- function(t, y, main = "", ...){
    plot(t, y, main = main, ...,
             type= "l", cex = 2, col = "black", lwd = linewidth) #pch = 20,
    grid(nx = NA)
    # grid(ny = 4, nx = NA)
}

#' get_pheno
#'
#' Get yearly vegetation phenological metrics of a curve fitting method
#'
#' @inheritParams get_GOF
#' @inheritParams D
#' @param x One of:
#' - `rfit` (rought fitting object), returned by [brks2rfit()].
#' - `fFITs` (fine fitting object), return by multiple curve fitting methods by [curvefit()] for 
#'    a growing season.
#' - list of [fFITs()] object, for multiple growing seasons.
#' @param method Which fine curve fitting method to be extracted?
#' @param TRS Threshold for `PhenoTrs`.
#' @param IsPlot Boolean. Whether to plot figure?
#' @param show.title Whether to show the name of fine curve fitting method
#' in top title?
#' @param ... ignored.
#'
#' @return List of every year phenology metrics
#'
#' @example inst/examples/ex-get_fitting_param_GOF.R
#' @export
get_pheno <- function(x, ...) UseMethod("get_pheno", x)

#' @inheritParams PhenoTrs
#' @rdname get_pheno
#' @export
get_pheno.rfit <- function(x, TRS = c(0.2, 0.5), asymmetric = TRUE, ...) {
    dt = x$dt
    tout = x$tout
    doys = as.numeric(difftime(tout, date.origin, units = "day"))
    TRS %<>% set_names(., paste0("TRS", .*10))
    pheno = dt %>% group_by(flag) %>% group_map(function(d, .y) {
        ind = which(tout >= d$beg & tout <= d$end)
        t = doys[ind]
        values = last(x$zs)[ind]
        der1 = diff(values) %>% c(NA, .)

        p_TRS = map(TRS, ~ PhenoTrs.default(values, t, trs = ., IsPlot = FALSE, asymmetric = asymmetric))
        p_DER = PhenoDeriv.default(values, t, der1, IsPlot = FALSE)
        p_GU  = PhenoGu.default(values, t, der1, IsPlot = FALSE)
        c(p_TRS, list(DER = p_DER, p_GU)) %>% unlist()
    }) %>% set_names(dt$flag)
    tidy_pheno(pheno)
}

#' @rdname get_pheno
#' @export
get_pheno.list <- function(x, method,
    TRS = c(0.2, 0.5, 0.6),
    analytical = TRUE, smoothed.spline = FALSE,
    IsPlot = FALSE, show.title = TRUE, ...)
{
    if (!is.list(x)) stop("Unsupported input type!")
    if (class(x) == 'fFITs')  x <- list(x)

    names   <- names(x) # methods
    methods <- if (missing(method)) names(x[[1]]$model) else method

    # pheno_list
    res <- lapply(set_names(seq_along(methods), methods), function(k){
        method <- methods[k]
        if (IsPlot){
            op <- par(mfrow = c(length(x), 5),
                mgp = c(3, 0.6, 0), mar = rep(0, 4), yaxt = "n", xaxt = "n")
            if (isTRUE(all.equal(par("oma"), c(0, 0, 0, 0)))) {
                margin_l = 5.5
                oma <- if (show.title) c(1, margin_l, 4, 1) else c(1, margin_l, 2, 1)
                par(oma = oma)
            }
        }

        pheno_list <- lapply(seq_along(names) %>% set_names(names), function(i){
            fFITs <- x[[i]]
            .params = listk(
                x = fFITs, method, TRS,
                analytical, smoothed.spline, IsPlot,
                show.PhenoName = ifelse(i == 1, TRUE, FALSE),
                title.left = names[i]
            )
            do.call(get_pheno.fFITs, .params)
        })

        if (IsPlot){
            par(new = TRUE, mfrow = c(1, 1), oma = rep(0, 4), mar = rep(0, 4))
            plot(0, axes = F, type = "n", xaxt = "n", yaxt = "n") #
            if (show.title) {
                text(1, 1, method, font = 2, cex = 1.3)
            }
        }
        pheno_list
    })
    lapply(res, tidy_pheno) %>% purrr::transpose()
}

#' @param fFITs `fFITs` object returned by [curvefits()]
#' @param title.left String of growing season flag.
#' @param show.PhenoName Whether to show phenological methods names in the top panel?
#'
#' @rdname get_pheno
#' @export
get_pheno.fFITs <- function(x, method,
    TRS = c(0.2, 0.5),
    analytical = TRUE, smoothed.spline = FALSE,
    IsPlot = FALSE,
    title.left = "", show.PhenoName = TRUE, ...)
{
    meths <- names(x$model)
    if (missing(method)) {
        method <- meths[1]
        warning(sprintf("method is missing and set to %s!", method))
    }
    if (!(method %in% meths)){
        warning(sprintf("%s not in methods and set to %s!", method, meths[1]))
        method <- meths[1]
    }
    model <- x$model[[method]]

    # get_pheno methods
    methods  <- c(paste0("TRS", TRS*10),"DER","GU", "ZHANG")
    TRS_last <- last(TRS) # only display last threshold figure

    ypred  <- last2(model$zs)
    all_na <- all(is.na(ypred))

    show.legend = FALSE
    if (IsPlot && !all_na){
        ti <- x$data$t
        yi <- x$data$y # may have NA values.
        # constrain plot ylims
        ylim0    <- c( pmin(min(yi, na.rm = TRUE), min(ypred)),
                       pmax(max(yi, na.rm = TRUE), max(ypred)))
        A        <- diff(ylim0)
        ylim     <- ylim0 + c(-1, 0.2) * 0.05 *A
        ylim_trs <- (ylim - ylim0) / A # TRS:0-1

        PhenoPlot(x$tout, ypred, ylim = ylim, yaxt = "s")
        lines(ti, yi, lwd = 1, col = "grey60")

        QC_flag <- x$data$QC_flag
        # Different quality points with different color and shape.
        if (is.null(QC_flag)){
            points(ti, yi)
        } else {
            for (j in seq_along(qc_levels)){
                ind = which(QC_flag == qc_levels[j])
                if (!is_empty(ind)) {
                    points(ti[ind], yi[ind], pch = qc_shapes[j],
                           col = qc_colors[j], bg = qc_colors[j])
                }
            }
        }

        # show legend in first subplot
        if (show.PhenoName){
            show.legend <- TRUE
            legend('topright', c('y', "f(t)"), lty = c(1, 1), pch =c(1, NA), bty='n')
        }
        stat <- get_GOF.fFITs(x)
        stat <- subset(stat, meth == method) %>% select(-meth) %>% unlist() %>% round(3)

        exprs = list(
            eval(substitute(bquote(R2 == .R2), list(.R2 = stat['R2']))),
            eval(substitute(bquote(NSE == .NSE), list(.NSE = stat['NSE']))),
            eval(substitute(bquote(RMSE == .RMSE), list(.RMSE = stat['RMSE'])))
        )

        legend('topleft', do.call(expression, exprs), adj = c(0.2, 0.2), bty='n', text.col = "red")
        # mtext(title.left, side = 2, line = 0.2)
        mtext(title.left, side = 2, line = 1.8)
    }
    if (show.PhenoName && IsPlot) mtext("Fine fitting", line = 0.2)

    p_TRS <- lapply(TRS, function(trs) {
        PhenoTrs(model, t = x$tout, approach = "White", trs = trs, IsPlot = FALSE)
    })

    if (IsPlot && !all_na) {
        p_TRS_last <- PhenoTrs(model, approach="White", trs=TRS_last, IsPlot = IsPlot, ylim = ylim)
        if (show.PhenoName) mtext(sprintf('TRS%d', TRS_last*10))
    }

    param_common  <- list(model, t = x$tout,
        analytical = analytical, smoothed.spline = smoothed.spline,
        IsPlot, ylim = ylim)
    param_common2 <- c(param_common, list(show.legend = show.legend))

    der <- do.call(PhenoDeriv, param_common2)
    if (show.PhenoName && IsPlot) mtext("DER", line = 0.2)

    zhang <- do.call(PhenoKl, param_common2)
    if (show.PhenoName && IsPlot) mtext("Inflexion ", line = 0.2)

    gu <- do.call(PhenoGu, param_common)[1:4]
    if (show.PhenoName && IsPlot) mtext("Gu", line = 0.2)

    c(p_TRS, list(der, gu, zhang)) %>% set_names(methods)
}

#' tidy_pheno
#'
#' Tidy for every method with multiple years phenology data
#'
#' @param pheno Phenology metrics extracted from `get_pheno`
#'
#' @keywords internal
#' @examples
#' library(phenofit)
#' # simulate vegetation time-series
#' fFUN <- doubleLog.Beck
#' par <- c(mn = 0.1, mx = 0.7, sos = 50, rsp = 0.1, eos = 250, rau = 0.1)
#'
#' t <- seq(1, 365, 8)
#' tout <- seq(1, 365, 1)
#' y <- fFUN(par, t)
#'
#' methods <- c("AG", "Beck", "Elmore", "Gu", "Zhang") # "Klos" too slow
#' fFITs <- curvefit(y, t, tout, methods)
#'
#' # multiple years
#' fits <- list(`2001` = fFITs, `2002` = fFITs)
#' pheno <- get_pheno(fits, "AG", IsPlot = FALSE)
#' @export
tidy_pheno <- function(pheno) {
    doy2date <- function(datenum) as.Date(unlist(datenum), origin = date.origin)
    names <- unlist(pheno[[1]]) %>% names()

    p_date <- map_df(pheno, function(x) {
        doy2date(x) %>% set_names(names) %>% as.list() %>% as.data.table()
    }, .id = "flag") %>%
        mutate(origin = ymd(paste0(substr(flag, 1, 4), "-01-01"))) %>%
        reorder_name(c("flag", "origin")) %>%
        set_colnames(c("flag", "origin", names)) %>%
        data.table()

    colnames(p_date) %<>% gsub("GU\\.|ZHANG\\.", "", .)
    phenonames <- setdiff(colnames(p_date), c("flag", "origin", "meth"))

    p_doy <- p_date %>% date2doy()
    vars <- c("flag", "origin", phenonames)
    list(doy = p_doy[, ..vars], date = p_date[, ..vars])
}

# phenonames <- c('TRS2.sos', 'TRS2.eos', 'TRS5.sos', 'TRS5.eos', 'TRS6.sos', 'TRS6.eos',
#                 'DER.sos', 'DER.pop', 'DER.eos',
#                 'GU.UD', 'GU.SD', 'GU.DD', 'GU.RD',
#                 'ZHANG.Greenup', 'ZHANG.Maturity', 'ZHANG.Senescence', 'ZHANG.Dormancy')
# fix the error of \code{pheno} without name
# names <- names(pheno)
# if (is.null(names)){
#     years <- seq(year(origin), by = 1, length.out = length(pheno))
#     names(pheno) <- years
# }

#' @rdname tidy_pheno
#' @export
date2doy <- function(p_date){
    p_date %>% mutate(across(3:ncol(.), ~ as.numeric(.x - origin + 1)))
}

#' @rdname tidy_pheno
#' @export
doy2date <- function(p_doy){
    p_doy %>% mutate(across(3:ncol(.), ~ (origin + .x - 1)))
}
