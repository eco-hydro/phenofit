ExtractPheno <- function(fits, TRS = c(0.1, 0.2, 0.5, 0.6), IsPlot = FALSE){
    names <- names(fits)
    pheno_list <- list()
    methods    <- c(paste0("TRS", TRS*10),"DER","GU", "ZHANG")
    TRS_last   <- last(TRS) # only display last threshold figure

    if (IsPlot)
        op <- par(mfrow = c(length(fits), 5),
                  oma = c(1, 2, 3, 1), mar = rep(0, 4), yaxt = "n", xaxt = "n")
    ylim <- NULL
    for (i in seq_along(fits)){
        fit    <- fits[[i]]
        ypred  <- last(fit$fits)
        all_na <- all(is.na(ypred))
        # 1. show curve fitting RMSE
        # need to fix here, about status variable. 31 Jan, 2018
        show.lgd = FALSE
        if (IsPlot && !all_na){
            ti <- fit$data$t
            yi <- fit$data$y

            ylim0    <- range(yi, na.rm = T); A = diff(ylim0);
            ylim     <- ylim0 + c(-1, 1) * 0.05 *A
            ylim_trs <- (ylim - ylim0[1]) / A # TRS:0-1

            PhenoPlot(fit$tout, ypred, ylim = ylim)
            lines(ti, yi, lwd = 1, col = "grey60")
            # pch = 19, col = "grey60"
            wi <- as.numeric(fit$data$w)

            ## Just designed for MOD13A1
            # Levels:  good  margin  snow&ice  cloud
            labels <- c(" good", " margin", " snow&ice", " cloud")
            colors <- c("grey60", "#00BFC4", "#F8766D", "#C77CFF")
            pch <- c(19, 15, 4, 17)
            for (j in 1:4){
                ind = which(wi == j)
                if (!is_empty(ind))
                    points(ti[ind], yi[ind], pch = pch[j], col = colors[j])
            }

            if (i == 1){
                show.lgd = TRUE
                legend('topright', c('y', "f(t)"), lty = c(1, 1), pch =c(1, NA), bty='n')
            }
            stat     <- statistic.phenofit(fit)
            stat_txt <- sprintf("  R=%.2f, p=%.3f\n RMSE=%.3f\nNSE=%.2f\n",
                                stat[['rmse']], stat[['nash']], stat[['R']], stat[['pvalue']])
            legend('topleft', stat_txt, adj = c(0.2, 0.2), bty='n', text.col = "red")
            mtext(names[i], side = 2)
        }
        if (i == 1 && IsPlot) mtext("fitting")

        p_TRS <- map(TRS, function(trs) {
            PhenoTrs(fit, approach = "White", trs = trs, IsPlot = FALSE)
        })

        if (IsPlot && !all_na) {
            trs6 <- PhenoTrs(fit, approach="White", trs=TRS_last, IsPlot = IsPlot, ylim = ylim_trs)
            if (i == 1) mtext(sprintf('TRS%d', TRS_last*10))
        }

        param_common  <- list(fit, IsPlot, ylim = ylim)
        param_common2 <- list(fit, IsPlot, ylim = ylim, show.lgd = show.lgd)

        der   <- do.call(PhenoDeriv, param_common2);   if (i == 1 && IsPlot) mtext("DER")
        gu    <- do.call(PhenoGu, param_common)[1:4];  if (i == 1 && IsPlot) mtext("GU")
        zhang <- do.call(PhenoKl, param_common2);      if (i == 1 && IsPlot) mtext("ZHANG")
        pheno_list[[i]] <- c(p_TRS, list(der, gu, zhang)) %>% set_names(methods)
    }
    pheno_list %<>% set_names(names)
    return(pheno_list)
}
