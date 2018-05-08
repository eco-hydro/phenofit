#' getFits_pheno
#'
#' Get yearly vegetation phenological metrics of a curve fitting method
#'
#' @param outdir save calculated phenology in oudir txts
#' @export
getFits_pheno <- function(fits, TRS = c(0.1, 0.2, 0.5, 0.6), IsPlot = FALSE){
    names <- names(fits)
    pheno_list <- list()
    methods    <- c(paste0("TRS", TRS*10),"DER","GU", "ZHANG")
    if (IsPlot)
        op <- par(mfrow = c(length(fits), 5),
            oma = c(1, 2, 3, 1), mar = rep(0, 4), yaxt = "n", xaxt = "n")
    for (i in seq_along(fits)){
        fit    <- fits[[i]]
        # 1. show curve fitting RMSE
        # need to fix here, about status variable. 31 Jan, 2018
        if (IsPlot){
            PhenoPlot(index(fit$pred), fit$pred)
            lines(fit$data$t, fit$data$x, type = "b")
            legend(min(fit$data$t) - 40, max(fit$pred),
                sprintf("rmse = %.2e\nnash = %.2f\nR=%.2f, p = %.3f",
                    status$rmse, status$nash, status$r, status$pvalue),
                adj = c(0, 0), bty='n', text.col = "red")
            mtext(names[i], side = 2)
        }
        if (i == 1 && IsPlot) mtext("Fitting")

        p_TRS <- map(TRS, function(trs) {
            PhenoTrs(fit$pred, approach = "White", trs = trs, IsPlot = FALSE)
        })

        if (i == 1 && IsPlot) {
            trs  <- last(TRS)
            trs6 <- PhenoTrs(last(fit$fits), approach="White", trs=trs, IsPlot = IsPlot)
            mtext(sprintf('TRS%d', trs*10))
        }

        der   <- PhenoDeriv(fit, IsPlot);    if (i == 1 && IsPlot) mtext("DER")
        gu    <- PhenoGu(fit, IsPlot)[1:4];  if (i == 1 && IsPlot) mtext("GU")
        zhang <- PhenoKl(fit, IsPlot);       if (i == 1 && IsPlot) mtext("ZHANG")
        pheno_list[[i]] <- c(p_TRS, list(der, gu, zhang)) %>% set_names(methods)
    }
    pheno_list %<>% set_names(names)
    return(pheno_list)
}

#' getYears_pheno
#'
#' Deprecated function.
#'
#' Get multi-years vegetation phenology metrics
#' @param outdir save calculated phenology in oudir txts
#' @export
getYears_pheno <- function(lst, tt, nptperyear, nyear, years, outdir = "./", IsPlot = FALSE, IsSave = FALSE){
    xx <- lst$xx
    i  <- lst$i
    # check outfile whether exist
    outfile <- sprintf("%spheno_%05d.rda", outdir, i)
    if (file.exists(outfile)){
        warning(sprintf('[%s] already exist!\n', outfile))
        return()
    }

    if (all(is.na(xx))){
        warning(sprintf("\t[%d]:all values are na!", i))
    }else{
        pheno <- list()
        for (j in 1:nyear){
            # runningId(j)
            cat(sprintf('running i = %d, j = %d\n', i, j))
            I <- ((j - 1)*nptperyear + 1) : (j*nptperyear)
            t <- tt[I]
            x <- xx[I]
            tout <- t[1]:t[length(t)]

            if (IsPlot) CairoPNG(filename = sprintf("%sTP_avhrr_pheno%05d.png", outdir, i), width = 960*2, height = 960*2*3/4, res = 200)
            pheno[[j]] <- tryCatch(
                {
                    fits <- curvefit(x, t, tout) #x didn't need to be a zoo object
                    getFits_pheno(fits, IsPlot = IsPlot)#return pheno
                },
                # warning = function(w) {
                #     cat(sprintf("[%d]WARNING %s", j, w))
                #     # suppressWarnings(PhenoFUN(xs, t, tout, ...))
                # },
                error = function(e){
                    # NULL will lead to inconsistent list length
                    sprintf("[%d]ERROR %s", j, e); NA
                }
            )
            if (IsPlot){
                title(main = list(sprintf("i = %05d, j = %02d", i, j), cex = 1.5, col = "red"), outer = T, line = 1.5)
                dev.off()
            }
        }
        pheno %<>% set_names(years)
        # should year by year
        if (IsSave){
            save(pheno, file = outfile)
        }else{
            return(pheno)
        }
    }
}
