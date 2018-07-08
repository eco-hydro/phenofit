
<!-- README.md is generated from README.Rmd. Please edit that file -->
phenofit
--------

[![Travis Build Status](https://travis-ci.org/kongdd/phenofit.svg?branch=master)](https://travis-ci.org/kongdd/phenofit) [![codecov](https://codecov.io/gh/kongdd/phenofit/branch/master/graph/badge.svg)](https://codecov.io/gh/kongdd/phenofit)

A state-of-the-art **remote sensing vegetation phenology** extraction package: `phenofit`

-   `phenofit` combine merits of TIMESAT and phenopix
-   A simple and stable growing season dividing methods was proposed
-   Provide a practical snow elimination method, based on Whittaker
-   7 curve fitting methods and 4 phenology extraction methods
-   We add parameters boundary for every curve fitting methods according to their ecological meaning.
-   `optimx` is used to select best optimization method for different curve fitting methods.

Installation
------------

You can install phenofit from github with:

``` r
# install.packages("devtools")
devtools::install_github("kongdd/phenofit")
```

Example
-------

Here, we illustrate how to use `phenofit` to extract vegetation phenology from MOD13A1 in the sampled points. Regional analysis also can be conducted in the similar way.

1. Preparing input data
=======================

1.1 Download MOD13A1 data
-------------------------

Upload point shapefile into GEE, clip MOD13A1 and download vegetation index data. [Here](https://code.earthengine.google.com/ee3ec39cf3061374dab435c358d008a3) is the corresponding GEE script.

1.2 Initial weights for input data
----------------------------------

Load packages.

``` r
library(phenofit)
library(data.table)
library(magrittr)
library(lubridate)
library(purrr)
library(plyr)
```

Read the point shapefile to get points coordinate information. Read Enhanced Vegetation Index (EVI) exported by `GEE`.

-   Add date according to composite day of the year (DayOfYear), other than image date.
-   Add weights according to `SummaryQA`.

For MOD13A1, Weights can by initialed by `SummaryQA` band (also suit for MOD13A2 and MOD13Q1). We have written a qc function for `SummaryQA`, `qc_summary`.

| SummaryQA       | Pixel reliability summary QA                        | weight |
|-----------------|-----------------------------------------------------|--------|
| -1 Fill/No data | Not processed                                       | `wmin` |
| 0 Good data     | Use with confidence                                 | 1      |
| 1 Marginal data | Useful but look at detailed QA for more information | 0.5    |
| 2 Snow/ice      | Pixel covered with snow/ice                         | `wmin` |
| 3 Cloudy        | Pixel is cloudy                                     | `wmin` |

``` r
data('MOD13A1')
df <- MOD13A1$dt
st <- MOD13A1$st

df[, `:=`(date = ymd(date), year = year(date), doy = as.integer(yday(date)))]
df[is.na(DayOfYear), DayOfYear := doy] # If DayOfYear is missing
    
# In case of last scene of a year, doy of last scene could in the next year
df[abs(DayOfYear - doy) >= 300, t := as.Date(sprintf("%d-%03d", year+1, DayOfYear), "%Y-%j")] # last scene
df[abs(DayOfYear - doy) <  300, t := as.Date(sprintf("%d-%03d", year  , DayOfYear), "%Y-%j")]

# MCD12Q1.006 land cover 1-17, IGBP scheme
IGBPnames_006 <- c("ENF", "EBF", "DNF", "DBF", "MF" , "CSH", 
              "OSH", "WSA", "SAV", "GRA", "WET", "CRO", 
              "URB", "CNV", "SNOW", "BSV", "water", "UNC")
# Initial weights
df[, w := qc_summary(SummaryQA)]
# Remap SummaryQA factor level, plot_phenofit use this variable. For other 
# remote sensing data without `SummaryQA`, need to modify `plot_phenofit`
if ('SummaryQA' %in% colnames(df)){
    values <- c("0", "1", "2", "3")
    levels <- c("good", "margin", "snow&ice", "cloud")
    df$SummaryQA %<>% factor() %>% mapvalues(values, levels)
}

df <- df[, .(site, y = EVI/1e4, t, w, date, SummaryQA)]
```

Add one year in head and tail, for growing season dividing. For example, the input data period is 20000218 ~ 20171219. After adding one year in head and tail, it becomes 19990101 ~ 20181219.

2. Running `phenofit`
=====================

2.1 Set global parameters for `phenofit`
----------------------------------------

``` r
# lambda   <- 5    # non-parameter Whittaker, only suit for 16-day. Other time-scale
# should assign a lambda.
ymax_min   <- 0.1  # the maximum ymax shoud be greater than `ymax_min` 
rymin_less <- 0.8  # trough < ymin + A*rymin_less
nptperyear <- 23   # How many points for a single year
wFUN       <- wTSM # Weights updating function, could be one of `wTSM`, 'wBisquare', `wChen` and `wSELF`. 

sites        <- unique(df$site)
sitename     <- df$site[1]
d            <- df[site == sitename] # get the first site data
sp           <- st[site == sitename]

print      <- TRUE
IsPlot     <- TRUE # for brks

prefix_fig <- "phenofit"
titlestr   <- with(sp, sprintf('[%03d,%s] %s, lat = %5.2f, lon = %6.2f',
                                     ID, site, IGBPname, lat, lon))
file_pdf   <- sprintf('Figure/%s_[%03d]_%s.pdf', prefix_fig, sp$ID[1], sp$site[1])
```

If need night temperature (Tn) to constrain ungrowing season backgroud value, NA values in Tn should be filled.

``` r
d$Tn %<>% zoo::na.approx(maxgap = 4)
plot(d$Tn, type = "l"); abline(a = 5, b = 0, col = "red")
```

2.1 Check input data
--------------------

``` r
dnew  <- add_HeadTail(d) # add additional one year in head and tail
INPUT <- check_input(dnew$t, dnew$y, dnew$w, maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
# y0 is used for plot. Original y value has been interpolated and changed.
INPUT$y0 <- dnew$y 
```

2.2 Divide growing seasons
--------------------------

Simply treating calendar year as a complete growing season will induce a considerable error for phenology extraction. A simple growing season dividing method was proposed in `phenofit`.

The growing season dividing method rely on heavily in Whittaker smoother.

Procedures of initial weight, growing season dividing and curve fitting are separated. Phenology extraction and curve fitting are combined together.

``` r
par(mar = c(3, 2, 2, 1), mgp = c(3, 0.6, 0))
lambda <- init_lambda(INPUT$y)
# The detailed information of those parameters can be seen in `season`.
# brks   <- season(INPUT, nptperyear,
#                FUN = whitsmw2, wFUN = wFUN, iters = 2,
#                lambda = lambda,
#                IsPlot = IsPlot, plotdat = d,
#                south = d$lat[1] < 0,
#                rymin_less = 0.6, ymax_min = ymax_min,
#                max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5) #, ...
# get growing season breaks in a 3-year moving window
brks2 <- season_3y(INPUT, nptperyear, south = sp$lat[1] < 0, FUN = whitsmw2,
                   IsPlot = IsPlot, print = print, partial = F)
#   [season_3y]  running 1 ...
#   [season_3y]  running 2 ...
#   [season_3y]  running 3 ...
#   [season_3y]  running 4 ...
#   [season_3y]  running 5 ...
#   [season_3y]  running 6 ...
#   [season_3y]  running 7 ...
#   [season_3y]  running 8 ...
#   [season_3y]  running 9 ...
#   [season_3y]  running 10 ...
#   [season_3y]  running 11 ...
#   [season_3y]  running 12 ...
#   [season_3y]  running 13 ...
#   [season_3y]  running 14 ...
#   [season_3y]  running 15 ...
#   [season_3y]  running 16 ...
#   [season_3y]  running 17 ...
#   [season_3y]  running 18 ...
```

<img src="man/Figure/divide growing season-1.png" style="display: block; margin: auto;" />

2.3 Curve fitting
-----------------

``` r
fit  <- curvefits(INPUT, brks2, lambda =lambda,
                  methods = c("AG", "zhang", "beck", "elmore"), #,"klos",, 'Gu'
                  nptperyear = nptperyear, debug = F, wFUN = wTSM,
                  nextent = 2, maxExtendMonth = 3, minExtendMonth = 1,
                  qc = as.numeric(dnew$SummaryQA), minPercValid = 0.2,
                  print = print)
#   [curvefits]  running 1 ...
#   [curvefits]  running 2 ...
#   [curvefits]  running 3 ...
#   [curvefits]  running 4 ...
#   [curvefits]  running 5 ...
#   [curvefits]  running 6 ...
#   [curvefits]  running 7 ...
#   [curvefits]  running 8 ...
#   [curvefits]  running 9 ...
#   [curvefits]  running 10 ...
#   [curvefits]  running 11 ...
#   [curvefits]  running 12 ...
#   [curvefits]  running 13 ...
#   [curvefits]  running 14 ...
#   [curvefits]  running 15 ...
# Warning in optim_pheno(prior, FUN, y, t, tout, optimFUN, method, w, lower =
# lower, : Not convergent!

# Warning in optim_pheno(prior, FUN, y, t, tout, optimFUN, method, w, lower =
# lower, : Not convergent!
#   [curvefits]  running 16 ...
# Warning in optim_pheno(prior, FUN, y, t, tout, optimFUN, method, w, lower =
# lower, : Not convergent!

# Warning in optim_pheno(prior, FUN, y, t, tout, optimFUN, method, w, lower =
# lower, : Not convergent!
#   [curvefits]  running 17 ...
#   [curvefits]  running 18 ...
fit$INPUT   <- INPUT
fit$seasons <- brks2

## check the curve fitting parameters
params <- getparam(fit)
print(str(params, 1))
# List of 4
#  $ AG    :Classes 'tbl_df', 'tbl' and 'data.frame':   18 obs. of  8 variables:
#  $ BECK  :Classes 'tbl_df', 'tbl' and 'data.frame':   18 obs. of  7 variables:
#  $ ELMORE:Classes 'tbl_df', 'tbl' and 'data.frame':   18 obs. of  8 variables:
#  $ ZHANG :Classes 'tbl_df', 'tbl' and 'data.frame':   18 obs. of  8 variables:
# NULL
print(params$AG)
# # A tibble: 18 x 8
#    flag      t0    mn    mx     rsp    a3     rau    a5
#    <fct>  <dbl> <dbl> <dbl>   <dbl> <dbl>   <dbl> <dbl>
#  1 2000_1  604. 0.418 0.621 0.00630  6    0.0233   6   
#  2 2001_1  913. 0.404 0.647 0.0198   2    0.0123   2   
#  3 2002_1 1247. 0.419 0.624 0.0228   2    0.00930  2   
#  4 2003_1 1601. 0.431 0.592 0.0406   2    0.00870  2   
#  5 2004_1 2088. 0.429 0.564 0.00655  2    0.0305   2   
#  6 2005_1 2430. 0.425 0.616 0.0102   2    0.0440   2   
#  7 2006_1 2790. 0.428 0.660 0.0418   3.50 0.0309   2   
#  8 2007_1 3145. 0.386 0.562 0.00965  2.18 0.0128   2   
#  9 2008_1 3520. 0.391 0.570 0.00870  2.00 0.0214   2   
# 10 2009_1 3861. 0.398 0.615 0.00767  4.45 0.0121   6   
# 11 2010_1 4247. 0.420 0.628 0.00984  2    0.0155   2   
# 12 2011_1 4526. 0.402 0.623 0.0289   2    0.00819  2   
# 13 2012_1 4905. 0.402 0.621 0.0158   2    0.00645  6   
# 14 2013_1 5237. 0.458 0.599 0.0394   2    0.00649  6   
# 15 2014_1 5724. 0.426 0.608 0.00828  2    0.0346   2   
# 16 2015_1 6122. 0.425 0.599 0.00522  6    0.0427   6   
# 17 2016_1 6334. 0.427 0.581 0.0186   2    0.00670  2.55
# 18 2017_1 6695. 0.424 0.561 0.0260   2    0.00573  6

## Get GOF information
stat  <- ldply(fit$fits, function(fits_meth){
    ldply(fits_meth, statistic.phenofit, .id = "flag")
}, .id = "meth")
fit$stat <- stat
print(head(stat))
#   meth   flag      RMSE         NSE         R       pvalue  n
# 1   AG 2000_1 0.2632999 -0.13463590 0.7813034 2.462271e-06 26
# 2   AG 2001_1 0.1624429  0.27917479 0.8328435 2.798427e-06 21
# 3   AG 2002_1 0.1810090  0.15679409 0.8284686 1.051100e-06 23
# 4   AG 2003_1 0.1824374 -0.01015221 0.7843137 5.738810e-06 24
# 5   AG 2004_1 0.2396933 -0.22206966 0.8652529 4.848286e-08 24
# 6   AG 2005_1 0.1931176  0.12270772 0.7317145 2.456480e-04 20

## visualization
# svg("Figure1_phenofit_curve_fitting.svg", 11, 7)
# Cairo::CairoPDF(file_pdf, 11, 6) #
# dev.off()
g <- plot_phenofit(fit, d, titlestr)
# Warning: Removed 316 rows containing missing values (geom_point).
grid::grid.newpage(); grid::grid.draw(g)# plot to check the curve fitting
```

<img src="man/Figure/curve fitting-1.png" style="display: block; margin: auto;" />

2.4 Extract phenology.
----------------------

``` r
# pheno: list(p_date, p_doy)
p <- lapply(fit$fits, ExtractPheno)
pheno  <- map(p, tidyFitPheno, origin = INPUT$t[1]) %>% purrr::transpose()
fit$pheno  <- pheno

# ratio = 1.15
# file <- "Figure5_Phenology_Extraction_temp.pdf"
# cairo_pdf(file, 8*ratio, 6*ratio)
# temp <- ExtractPheno(fit$fits$ELMORE[2:6], IsPlot = T)
# dev.off()
# file.show(file)

## check the extracted phenology
temp <- ExtractPheno(fit$fits$ELMORE[1:6], IsPlot = T, TRS = 0.5)
```

<img src="man/Figure/Extract phenology-1.png" style="display: block; margin: auto;" />

``` r

print(str(pheno, 2))
# List of 2
#  $ date:List of 4
#   ..$ AG    :Classes 'tbl_df', 'tbl' and 'data.frame':    18 obs. of  21 variables:
#   ..$ BECK  :Classes 'tbl_df', 'tbl' and 'data.frame':    18 obs. of  21 variables:
#   ..$ ELMORE:Classes 'tbl_df', 'tbl' and 'data.frame':    18 obs. of  21 variables:
#   ..$ ZHANG :Classes 'tbl_df', 'tbl' and 'data.frame':    18 obs. of  21 variables:
#  $ doy :List of 4
#   ..$ AG    :Classes 'tbl_df', 'tbl' and 'data.frame':    18 obs. of  21 variables:
#   ..$ BECK  :Classes 'tbl_df', 'tbl' and 'data.frame':    18 obs. of  21 variables:
#   ..$ ELMORE:Classes 'tbl_df', 'tbl' and 'data.frame':    18 obs. of  21 variables:
#   ..$ ZHANG :Classes 'tbl_df', 'tbl' and 'data.frame':    18 obs. of  21 variables:
# NULL
head(pheno$doy$AG)
# # A tibble: 6 x 21
#   flag   origin     TRS1.sos TRS1.eos TRS2.sos TRS2.eos TRS5.sos TRS5.eos
#   <fct>  <date>        <dbl>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl>
# 1 2000_1 2000-01-01       57      291       68      287       91      281
# 2 2001_1 2001-01-01      105      310      119      287      141      251
# 3 2002_1 2002-01-01       86      320       98      290      117      243
# 4 2003_1 2003-01-01      103      321      110      289      121      237
# 5 2004_1 2004-01-01       42      316       69      306      137      292
# 6 2005_1 2005-01-01       87      274      114      268      158      258
# # ... with 13 more variables: TRS6.sos <dbl>, TRS6.eos <dbl>,
# #   DER.sos <dbl>, DER.pop <dbl>, DER.eos <dbl>, GU.UD <dbl>, GU.SD <dbl>,
# #   GU.DD <dbl>, GU.RD <dbl>, ZHANG.Greenup <dbl>, ZHANG.Maturity <dbl>,
# #   ZHANG.Senescence <dbl>, ZHANG.Dormancy <dbl>
```
