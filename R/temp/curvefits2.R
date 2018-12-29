##  curvefits2
## 
##  @param t A date vector
##  @param y A numberic vector, same length as t
##  @param methods A character vector, c('spline', 'beck', 'elmore', 'klos', 'AG', 'zhang')
##  @param ... other parameters pass to season or curvefit_site
##  @export
curvefits2 <- function(t, y, w, nptperyear = 46,
                          wFUN = wTSM, iters = 2,
                          lambda, south = FALSE,
                          IsPlot = FALSE,
                          Aymin_less = 0.6, ymax_min = 0.1,
                          methods = c('AG', 'zhang', 'beck', 'elmore', 'Gu'),
                          debug = FALSE, ...)
{
    # 1. Check input data and initial parameters for phenofit
    INPUT <- check_input(t, y, w, maxgap = nptperyear / 4, alpha = 0.02)

    # 2. The detailed information of those parameters can be seen in `season`.
    brks  <- season(INPUT, lambda, nptperyear, iters = 3, wFUN = wFUN,
                    IsPlot = IsPlot, south = south,
                    Aymin_less = Aymin_less, ymax_min = ymax_min,
                    max_MaxPeaksperyear=2.5, max_MinPeaksperyear=3.5, ...)
    ## 3.1. curve fitting
    fit <- curvefits(INPUT, brks, lambda =lambda,
                         methods = c("AG", "zhang", "beck", "elmore", 'Gu'), #,"klos"
                         nptperyear = nptperyear, debug = F, wFUN = wTSM, ...)
    ## 3.2 Get GOF information
    stat  <- ldply(fit$fits, function(fits_meth){
        ldply(fits_meth, statistic.phenofit, .id = "flag")
    }, .id = "meth")

    # 4. extract phenology
    # pheno: list(p_date, p_doy)
    p <- lapply(fit$fits, PhenoExtract)
    pheno  <- map(p, tidyFitPheno, origin = t[1]) %>% purrr::transpose()

    fit$INPUT   <- INPUT
    fit$seasons <- brks
    fit$stat    <- stat
    fit$pheno   <- pheno
    return(fit)
}
