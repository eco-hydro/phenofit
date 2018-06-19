library(phenofit)
library(data.table)
library(plyr)
library(data.table)

source("test/GEE/pkg_seasonality.R")

get_phenofit <- function(d){
    tryCatch({
        #d     <- df[Id == i, ] # get the first site data
        titlestr <- with(d[1, ], sprintf('[%03d]%s, lat = %.2f, lon = %.2f', 
                                         Id, IGBPname, lat, lon))
        dnew  <- add_HeadTail(d)
        # 1. Check input data and initial parameters for phenofit
        INPUT <- check_input(dnew$t, dnew$y, dnew$w, maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
        INPUT$y0 <- dnew$y
        IsPlot   <- FALSE # for brks
        # 2. The detailed information of those parameters can be seen in `season`.
        lambda <- init_lambda(INPUT$y)#*2
        # brks   <- season(INPUT, nptperyear,
        #                FUN = whitsmw2, wFUN = wFUN, iters = 2, 
        #                lambda = lambda,
        #                IsPlot = IsPlot, plotdat = d, 
        #                south = d$lat[1] < 0,
        #                rymin_less = 0.6, ymax_min = ymax_min,
        #                max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5) #, ...
        # get growing season breaks in a 3-year moving window
        brks2 <- season_3y(INPUT, nptperyear, FUN = whitsmw2, 
                           IsPlot = IsPlot, print = print, south = d$lat[1] < 0, partial = F)

        # 3. curve fitting
        fit  <- curvefits(INPUT, brks2, lambda =lambda,
                          methods = c("AG", "zhang", "beck", "elmore"), #,"klos",, 'Gu'
                          nptperyear = nptperyear, debug = F, wFUN = wTSM,
                          extend_month = 2, nextent = 1,
                          qc = as.numeric(dnew$SummaryQA), minPercValid = 0.2, 
                          print = print)
        fit$INPUT   <- INPUT
        fit$seasons <- brks2

        # svg("Figure1_phenofit_curve_fitting.svg", 11, 7)
        file_pdf = sprintf('Figures/phenofit_v3_%03d.pdf', d$Id[1])
        Cairo::CairoPDF(file_pdf, 11, 6)
        # grid::grid.newpage()
        plot_phenofit(fit, d, titlestr) %>% grid::grid.draw()# plot to check the curve fitting
        dev.off()
        
        # temp <- ExtractPheno(fit$fits$ELMORE[1:5], IsPlot = T) # check extracted phenology
        ## 3.2 Get GOF information
        stat  <- ldply(fit$fits, function(fits_meth){
            ldply(fits_meth, statistic.phenofit, .id = "flag")
        }, .id = "meth")
        fit$stat <- stat

        # 4. extract phenology, check extracted phenology extraction for one method.
        # ratio = 1.15
        # file <- "Figure5_Phenology_Extraction_temp.pdf"
        # cairo_pdf(file, 8*ratio, 6*ratio)
        # temp <- ExtractPheno(fit$fits$ELMORE[2:6], IsPlot = T)
        # dev.off()
        # file.show(file)

        # pheno: list(p_date, p_doy)
        p <- lapply(fit$fits, ExtractPheno)
        pheno  <- map(p, tidyFitPheno, origin = INPUT$t[1]) %>% purrr::transpose()

        fit$pheno  <- pheno
        return(fit)
    }, error = function(e){
        message(sprintf('[e]: %s, %s', titlestr, e$message))
    })    
}

################################################################################
# lambda   <- 5    # Whittaker parameter
ymax_min   <- 0.1  # the maximum ymax shoud be greater than `ymax_min` 
rymin_less <- 0.8  # trough < ymin + A*rymin_less
nptperyear <- 23   # How many points for a single year
wFUN       <- wTSM # Weights updating function, could be one of `wTSM`, 'wBisquare', `wChen` and `wSELF`. 

levels <- c(" good", " margin", " snow&ice", " cloud")
df <- fread("vignettes/phenofit_MOD13A1_st_20.csv", strip.white = F)
df[, `:=`( t = ymd(t), 
    SummaryQA = factor(SummaryQA, levels))]

Ids        <- unique(df$Id)
Id         <- df$Id[1]
d          <- df[Id == Ids[1], ] # get the first site data

par_sbatch(Ids, get_phenofit, df = df, save = T, outdir = "result_v3")
