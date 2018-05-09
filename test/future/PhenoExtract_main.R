
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
