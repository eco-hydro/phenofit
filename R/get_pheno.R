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

#' get_pheno
#'
#' Get yearly vegetation phenological metrics of a curve fitting method
#'
#' @inheritParams get_GOF
#' @inheritParams D
#' @param fits A list of [fFITs()] object, for a single curve fitting
#' method.
#' @param method Which fine curve fitting method to be extracted?
#' @param TRS Threshold for `PhenoTrs`.
#' @param IsPlot Boolean. Whether to plot figure?
#' @param show_title Whether to show the name of fine curve fitting method
#' in top title?
#' @param ... ignored.
#'
#' @note
#' Please note that only a single fine curve fitting method allowed here!
#'
#' @return List of every year phenology metrics
#'
#' @example inst/examples/ex-get_fitting_param_GOF.R
#' @export
get_pheno <- function(fits, method,
    TRS = c(0.2, 0.5, 0.6),
    analytical = TRUE, smoothed.spline = FALSE,
    IsPlot = FALSE, show_title = TRUE, ...)
{
    if (!is.list(fits)) stop("Unsupported input type!")
    if (class(fits) == 'fFITs')  fits <- list(fits)

    names   <- names(fits) # methods
    methods <- if (missing(method)) names(fits[[1]]$model) else method

    # pheno_list
    res <- lapply(set_names(seq_along(methods), methods), function(k){
        method <- methods[k]
        if (IsPlot){
            oma <- if (show_title) c(1, 2, 4, 1) else c(1, 2, 2, 1)
            op <- par(mfrow = c(length(fits), 5), oma = oma,
                mar = rep(0, 4), yaxt = "n", xaxt = "n")
        }

        # fFITs
        pheno_list <- lapply(seq_along(names) %>% set_names(names), function(i){
            fFITs <- fits[[i]]
            .params = listk(
                fFITs, method, TRS,
                analytical, smoothed.spline, IsPlot,
                showName_pheno = ifelse(i == 1, TRUE, FALSE),
                title_left = names[i]
            )
            do.call(get_pheno.fFITs, .params)
        })

        if (IsPlot){
            par(new = TRUE, mfrow = c(1, 1), oma = rep(0, 4), mar = rep(0, 4))
            plot(0, axes = F, type = "n", xaxt = "n", yaxt = "n") #
            if (show_title) {
                text(1, 1, method, font = 2, cex = 1.3)
            }
        }
        pheno_list
    })
    lapply(res, tidy_pheno) %>% purrr::transpose()
}

#' @param title_left String of growing season flag.
#' @param showName_pheno Whether to show phenological methods names in the top panel?
#'
#' @rdname get_pheno
#' @export
get_pheno.fFITs <- function(fFITs, method,
    TRS = c(0.2, 0.5),
    analytical = TRUE, smoothed.spline = FALSE,
    IsPlot = FALSE,
    title_left = "", showName_pheno = TRUE)
{
    meths <- names(fFITs$model)
    if (missing(method)) {
        method <- meths[1]
        warning(sprintf("method is missing and set to %s!", method))
    }
    if (!(method %in% meths)){
        warning(sprintf("%s not in methods and set to %s!", method, meths[1]))
        method <- meths[1]
    }
    fFIT <- fFITs$model[[method]]

    # get_pheno methods
    methods  <- c(paste0("TRS", TRS*10),"DER","GU", "ZHANG")
    TRS_last <- last(TRS) # only display last threshold figure

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
            for (j in seq_along(qc_levels)){
                ind = which(QC_flag == qc_levels[j])
                if (!is_empty(ind)) {
                    points(ti[ind], yi[ind], pch = qc_shapes[j],
                           col = qc_colors[j], bg = qc_colors[j])
                }
            }
        }

        # show legend in first subplot
        if (showName_pheno){
            show.lgd <- TRUE
            legend('topright', c('y', "f(t)"), lty = c(1, 1), pch =c(1, NA), bty='n')
        }
        stat <- get_GOF.fFITs(fFITs)
        stat <- subset(stat, meth == method) %>% select(-meth) %>% unlist() %>% round(3)

        exprs = list(
            eval(substitute(bquote(R2 == .R2), list(.R2 = stat['R2']))),
            eval(substitute(bquote(NSE == .NSE), list(.NSE = stat['NSE']))),
            eval(substitute(bquote(RMSE == .RMSE), list(.RMSE = stat['RMSE'])))
        )

        legend('topleft', do.call(expression, exprs), adj = c(0.2, 0.2), bty='n', text.col = "red")
        mtext(title_left, side = 2, line = 0.2)
    }
    if (showName_pheno && IsPlot) mtext("Fitting", line = 0.2)

    p_TRS <- lapply(TRS, function(trs) {
        PhenoTrs(fFIT, approach = "White", trs = trs, IsPlot = FALSE)
    })

    if (IsPlot && !all_na) {
        p_TRS_last <- PhenoTrs(fFIT, approach="White", trs=TRS_last, IsPlot = IsPlot, ylim = ylim)
        if (showName_pheno) mtext(sprintf('TRS%d', TRS_last*10))
    }

    param_common  <- list(fFIT,
        analytical = analytical, smoothed.spline = smoothed.spline,
        IsPlot, ylim = ylim)
    param_common2 <- c(param_common, list(show.lgd = show.lgd))

    der   <- do.call(PhenoDeriv, param_common2);  if (showName_pheno && IsPlot) mtext("DER", line = 0.2)
    gu    <- do.call(PhenoGu, param_common)[1:4]; if (showName_pheno && IsPlot) mtext("GU", line = 0.2)
    zhang <- do.call(PhenoKl, param_common2);     if (showName_pheno && IsPlot) mtext("ZHANG", line = 0.2)

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

date2doy <- function(p_date){
    p_date %>% mutate(across(3:ncol(.), ~ as.numeric(.x - origin + 1)))
}
