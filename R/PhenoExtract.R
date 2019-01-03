colors    <- c("blue", "green3", "orange", "red")
linewidth <- 1.2

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
    # grid(nx = NA)
    grid(ny = 4, nx = NA)
}

#' PhenoExtract
#'
#' Get yearly vegetation phenological metrics of a curve fitting method
#'
#' @inheritParams GOF_fFITs
#' @inheritParams D
#' @param fits A list of \code{\link{fFITs}} object, for a single curve fitting
#' method.
#' @param method Which fine curve fitting method to be extracted?
#' @param TRS Threshold for \code{PhenoTrs}.
#' @param IsPlot Boolean. Whether to plot figure?
#' @param title_left String of growing season flag.
#' @param show.title_top Whether to show top title? Generally, only show top
#' title in the first row.
#' @param ... ignored.
#'
#' @note
#' Please note that only a single fine curve fitting method allowed here!
#'
#' @return List of every year phenology metrics
#'
#' @examples
#' library(phenofit)
#' # simulate vegetation time-series
#' fFUN = doubleLog.Beck
#' par  = c(
#'     mn  = 0.1,
#'     mx  = 0.7,
#'     sos = 50,
#'     rsp = 0.1,
#'     eos = 250,
#'     rau = 0.1)
#' t    <- seq(1, 365, 8)
#' tout <- seq(1, 365, 1)
#' y <- fFUN(par, t)
#'
#' methods <- c("AG", "Beck", "Elmore", "Gu", "Zhang") # "Klos" too slow
#' fFITs <- curvefit(y, t, tout, methods)
#'
#' par(mfrow = c(1, 5),
#'     oma = c(1, 2, 3, 1), mar = rep(0, 4), yaxt = "n", xaxt = "n")
#' pheno <- PhenoExtract.fFITs(fFITs, "AG", IsPlot = TRUE)
#'
#' # multiple years
#' fits <- list(`2001` = fFITs, `2002` = fFITs)
#' pheno <- PhenoExtract(fits, "AG", IsPlot=TRUE)
#' @export
PhenoExtract <- function(fits, method,
    TRS = c(0.2, 0.5),
    analytical = TRUE, smoothed.spline = FALSE,
    IsPlot = FALSE, ...)
{
    if (!is.list(fits)) {
        stop("Unsupported input type!")
    }

    if (class(fits) == 'fFITs') {
        fits <- list(fits)
    }

    names   <- names(fits)
    # nseason <- length(fits)

    if (missing(method)) {
        methods <- names(fits[[1]]$fFIT)
    } else {
        methods <- method
    }

    # pheno_list
    res <- llply(set_names(seq_along(methods), methods), function(k){
        method <- methods[k]

        if (IsPlot){
            op <- par(mfrow = c(length(fits), 5),
                oma = c(1, 2, 4, 1), mar = rep(0, 4), yaxt = "n", xaxt = "n")
        }

        # fFITs
        pheno_list <- llply(seq_along(names) %>% set_names(names), function(i){
            fFITs <- fits[[i]]
            title_left     <- names[i]
            show.title_top <- ifelse(i == 1, TRUE, FALSE)

            # browser()
            PhenoExtract.fFITs(fFITs, method,
                TRS = TRS,
                analytical = analytical, smoothed.spline = smoothed.spline,
                IsPlot = IsPlot,
                title_left = title_left, show.title_top = show.title_top)
        })

        if (IsPlot){
            par(new = TRUE, mfrow = c(1, 1), oma = rep(0, 4), mar = rep(0, 4))
            plot(0, axes = F, type = "n", xaxt = "n", yaxt = "n") #
            text(1, 1, method, font = 2, cex = 1.3)
        }

        pheno_list
    })

    # return
    map(res, tidyFitPheno) %>% purrr::transpose() 
}

#' @rdname PhenoExtract
#' @export
PhenoExtract.fFITs <- function(fFITs, method,
    TRS = c(0.2, 0.5),
    analytical = TRUE, smoothed.spline = FALSE,
    IsPlot = FALSE,
    title_left = "", show.title_top = TRUE)
{
    meths <- names(fFITs$fFIT)
    if (missing(method)) {
        method <- meths[1]
        warning(sprintf("method is missing and set to %s!", method))
    }

    # browser()
    if (!(method %in% meths)){
        warning(sprintf("%s not in methods and set to %s!", method, meths[1]))
        method <- meths[1]
    }

    # PhenoExtract methods
    methods  <- c(paste0("TRS", TRS*10),"DER","GU", "ZHANG")
    TRS_last <- last(TRS) # only display last threshold figure

    fFIT   <- fFITs$fFIT[[method]]
    ypred  <- last(fFIT$zs)
    all_na <- all(is.na(ypred))

    show.lgd = FALSE

    if (IsPlot && !all_na){
        ti <- fFITs$data$t
        yi <- fFITs$data$y # may have NA values.
        # constrain plot ylims
        ylim0    <- c( pmin(min(yi, na.rm = TRUE), min(ypred)),
                       pmax(max(yi, na.rm = TRUE), max(ypred)))
        A        <- diff(ylim0)
        ylim     <- ylim0 + c(-1, 0.2) * 0.05 *A
        ylim_trs <- (ylim - ylim0) / A # TRS:0-1

        PhenoPlot(fFITs$tout, ypred, ylim = ylim)
        lines(ti, yi, lwd = 1, col = "grey60")

        QC_flag <- fFITs$data$QC_flag
        # Different quality points with different color and shape.
        if (is.null(QC_flag)){
            points(ti, yi)
        } else {
            # browser()
            for (j in seq_along(qc_levels)){
                ind = which(QC_flag == qc_levels[j])
                if (!is_empty(ind)) {
                    points(ti[ind], yi[ind], pch = qc_shapes[j],
                           col = qc_colors[j], bg = qc_colors[j])
                }
            }
        }

        # show legend in first subplot
        if (show.title_top){
            show.lgd <- TRUE
            legend('topright', c('y', "f(t)"), lty = c(1, 1), pch =c(1, NA), bty='n')
        }
        stat <- GOF_fFITs(fFITs)
        stat <- subset(stat, meth == method)

        stat_txt <- sprintf("  R=%.2f, p=%.3f\n RMSE=%.3f\nNSE=%.2f\n",
            stat[['R']], stat[['pvalue']], stat[['RMSE']], stat[['NSE']])

        legend('topleft', stat_txt, adj = c(0.2, 0.2), bty='n', text.col = "red")
        mtext(title_left, side = 2)
    }
    if (show.title_top && IsPlot) mtext("fitting")
    # browser()

    p_TRS <- map(TRS, function(trs) {
        PhenoTrs(fFIT, approach = "White", trs = trs, IsPlot = FALSE)
    })

    if (IsPlot && !all_na) {
        p_TRS_last <- PhenoTrs(fFIT, approach="White", trs=TRS_last, IsPlot = IsPlot, ylim = ylim)
        if (show.title_top) mtext(sprintf('TRS%d', TRS_last*10))
    }

    param_common  <- list(fFIT,
        analytical = analytical, smoothed.spline = smoothed.spline,
        IsPlot, ylim = ylim)
    param_common2 <- c(param_common, list(show.lgd = show.lgd))

    der   <- do.call(PhenoDeriv, param_common2);  if (show.title_top && IsPlot) mtext("DER")
    gu    <- do.call(PhenoGu, param_common)[1:4]; if (show.title_top && IsPlot) mtext("GU")
    zhang <- do.call(PhenoKl, param_common2);     if (show.title_top && IsPlot) mtext("ZHANG")

    c(p_TRS, list(der, gu, zhang)) %>% set_names(methods)
}

