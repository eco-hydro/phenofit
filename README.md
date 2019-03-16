
phenofit
========

[![Travis Build Status](https://travis-ci.org/kongdd/phenofit.svg?branch=master)](https://travis-ci.org/kongdd/phenofit)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/kongdd/phenofit?branch=master&svg=true)](https://ci.appveyor.com/project/kongdd/phenofit)
[![codecov](https://codecov.io/gh/kongdd/phenofit/branch/master/graph/badge.svg)](https://codecov.io/gh/kongdd/phenofit)
[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
[![CRAN](http://www.r-pkg.org/badges/version/phenofit)](https://cran.r-project.org/package=phenofit)

A state-of-the-art **remote sensing vegetation phenology** extraction package: `phenofit`

-   `phenofit` combine merits of TIMESAT and phenopix
-   A simple and stable growing season dividing methods was proposed
-   Provide a practical snow elimination method, based on Whittaker
-   7 curve fitting methods and 4 phenology extraction methods
-   We add parameters boundary for every curve fitting methods according to their ecological meaning.
-   `optimx` is used to select best optimization method for different curve fitting methods.


***Task lists***

- [ ] Improve computational efficiency of fine fitting;
- [ ] Complete script automatic generating module of shinyapp;
- [ ] Uncertainty analysis of curve fitting and phenological metrics;
- [ ] Support spatial analysis;
- [ ] Support annual season in curve fitting;
- [ ] flexible fine fitting input ( original time-series or smoothed time-series by rough fitting).


![title](man/Figure/Figure1_phenofit_flowchart.svg)   
*<u>Figure 1. The flowchart of phenology extraction in `phenofit`.</u>*

Installation
============

You can install phenofit from github with:

``` r
# install.packages("devtools")
devtools::install_github("kongdd/phenofit")
```

Or run with `shiny`:

``` r
shiny::runGitHub("phenofit", "kongdd", subdir = "inst/shiny/phenofit")
# Or run local installed
shiny::runApp(system.file("shiny/phenofit", package = "phenofit"))
```

R code Example
==============

Here, we illustrate how to use `phenofit` to extract vegetation phenology from MOD13A1 in the sampled points. Regional analysis also can be conducted in the similar way.

1.1 Download MOD13A1 data
-------------------------

Upload point shapefile into GEE, clip MOD13A1 and download vegetation index data. [Here](https://code.earthengine.google.com/ee3ec39cf3061374dab435c358d008a3) is the corresponding GEE script.

1.2 Initial weights for input data
----------------------------------

Load packages.

``` r
suppressMessages({
    library(data.table)
    library(magrittr)
    library(lubridate)
    library(purrr)
    library(plyr)
    
    library(phenofit)
})
```

Set global parameters for `phenofit`

``` r
# lambda   <- 5    # non-parameter Whittaker, only suit for 16-day. Other time-scale
# should assign a lambda.
ymax_min   <- 0.1  # the maximum ymax shoud be greater than `ymax_min` 
rymin_less <- 0.8  # trough < ymin + A*rymin_less
nptperyear <- 23   # How many points for a single year
wFUN       <- wBisquare #wTSM #wBisquare # Weights updating function, could be one of `wTSM`, 'wBisquare', `wChen` and `wSELF`. 
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

df <- df[!duplicated(df[, .(site, t)]), ]

# MCD12Q1.006 land cover 1-17, IGBP scheme
IGBPnames_006 <- c("ENF", "EBF", "DNF", "DBF", "MF" , "CSH", 
              "OSH", "WSA", "SAV", "GRA", "WET", "CRO", 
              "URB", "CNV", "SNOW", "BSV", "water", "UNC")
# Initial weights
df[, c("QC_flag", "w") := qc_summary(SummaryQA)]
# Remap SummaryQA factor level, plot_phenofit use this variable. For other 
# remote sensing data without `SummaryQA`, need to modify `plot_phenofit`

df <- df[, .(site, y = EVI/1e4, t, date, w, QC_flag)]
```

Add one year in head and tail, for growing season dividing. For example, the input data period is 20000218 ~ 20171219. After adding one year in head and tail, it becomes 19990101 ~ 20181219.

2.1 load data
-------------

``` r
sites        <- unique(df$site)
sitename     <- sites[3]
d            <- df[site == sitename] # get the first site data
sp           <- st[site == sitename]

south      <- sp$lat < 0
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

2.2 Check input data
--------------------

``` r
dnew  <- add_HeadTail(d, south, nptperyear = 23) # add additional one year in head and tail
INPUT <- check_input(dnew$t, dnew$y, dnew$w, dnew$QC_flag,
                     nptperyear, south, 
                     maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
```

2.3 Divide growing seasons
--------------------------

Simply treating calendar year as a complete growing season will induce a considerable error for phenology extraction. A simple growing season dividing method was proposed in `phenofit`.

The growing season dividing method rely on heavily in Whittaker smoother.

Procedures of initial weight, growing season dividing and curve fitting are separated. Phenology extraction and curve fitting are combined together.

``` r
par(mar = c(3, 2, 2, 1), mgp = c(3, 0.6, 0))
lambda <- init_lambda(INPUT$y)
# The detailed information of those parameters can be seen in `season`.
# brks   <- season(INPUT, nptperyear,
#                FUN = wWHIT, wFUN = wFUN, iters = 2,
#                lambda = lambda,
#                IsPlot = IsPlot, plotdat = d,
#                south = d$lat[1] < 0,
#                rymin_less = 0.6, ymax_min = ymax_min,
#                max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5) #, ...
# get growing season breaks in a 3-year moving window
brks2 <- season_mov(INPUT, 
                   FUN = wWHIT, wFUN = wFUN,
                   maxExtendMonth = 6, r_min = 0.1,
                   IsPlot = IsPlot, IsPlot.OnlyBad = FALSE, print = print)
#   [season_mov]  running 1 ...
#   [season_mov]  running 2 ...
#   [season_mov]  running 3 ...
#   [season_mov]  running 4 ...
#   [season_mov]  running 5 ...
#   [season_mov]  running 6 ...
#   [season_mov]  running 7 ...
#   [season_mov]  running 8 ...
#   [season_mov]  running 9 ...
#   [season_mov]  running 10 ...
#   [season_mov]  running 11 ...
#   [season_mov]  running 12 ...
#   [season_mov]  running 13 ...
#   [season_mov]  running 14 ...
#   [season_mov]  running 15 ...
#   [season_mov]  running 16 ...
#   [season_mov]  running 17 ...
#   [season_mov]  running 18 ...
```

<img src="man/Figure/divide growing season-1.svg" style="display: block; margin: auto;" />

2.4 Curve fitting
-----------------

``` r
fit  <- curvefits(INPUT, brks2,
                  methods = c("AG", "Zhang", "Beck", "Elmore"), #,"klos",, 'Gu'
                  wFUN = wFUN,
                  nextent = 2, maxExtendMonth = 3, minExtendMonth = 1, minPercValid = 0.2,
                  print = print, verbose = FALSE)
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
#   [curvefits]  running 16 ...
#   [curvefits]  running 17 ...
#   [curvefits]  running 18 ...

## check the curve fitting parameters
l_param <- get_param(fit)
print(str(l_param, 1))
# List of 4
#  $ AG    :Classes 'tbl_df', 'tbl' and 'data.frame':   18 obs. of  8 variables:
#  $ BECK  :Classes 'tbl_df', 'tbl' and 'data.frame':   18 obs. of  7 variables:
#  $ ELMORE:Classes 'tbl_df', 'tbl' and 'data.frame':   18 obs. of  8 variables:
#  $ ZHANG :Classes 'tbl_df', 'tbl' and 'data.frame':   18 obs. of  8 variables:
# NULL
print(l_param$AG)
# # A tibble: 18 x 8
#    flag      t0    mn    mx    rsp    a3    rau    a5
#    <fct>  <dbl> <dbl> <dbl>  <dbl> <dbl>  <dbl> <dbl>
#  1 2000_1  201. 0.168 0.407 0.0348  3.38 0.0148  6   
#  2 2001_1  564. 0.174 0.405 0.0215  4.48 0.0168  6   
#  3 2002_1  938. 0.189 0.494 0.0310  2.25 0.0178  3.57
#  4 2003_1 1276. 0.168 0.432 0.0263  2    0.0124  3.71
#  5 2004_1 1660. 0.181 0.449 0.0395  2    0.0181  4.79
#  6 2005_1 2054. 0.180 0.462 0.0132  6    0.0332  2   
#  7 2006_1 2383. 0.175 0.436 0.0195  2.12 0.0136  3.74
#  8 2007_1 2753. 0.166 0.483 0.0207  2    0.0151  2.89
#  9 2008_1 3118. 0.175 0.501 0.0265  2    0.0152  3.90
# 10 2009_1 3526. 0.174 0.485 0.0135  5.25 0.0319  2   
# 11 2010_1 3841. 0.190 0.492 0.0248  2.01 0.0150  2.09
# 12 2011_1 4208. 0.184 0.465 0.0308  2    0.0127  4.77
# 13 2012_1 4564. 0.170 0.511 0.0391  2.57 0.0115  6   
# 14 2013_1 4956. 0.171 0.489 0.0165  6    0.0153  2.20
# 15 2014_1 5305. 0.201 0.500 0.0342  2    0.0137  3.57
# 16 2015_1 5678. 0.217 0.504 0.0244  4.14 0.0201  6   
# 17 2016_1 6024. 0.195 0.485 0.0407  2    0.0126  4.76
# 18 2017_1 6405. 0.175 0.445 0.0210  6    0.0126  3.70

d_fit <- get_fitting(fit)
## Get GOF information
d_gof <- get_GOF(fit)
# fit$stat <- stat
print(head(d_gof))
#      flag   meth       RMSE       NSE         R       pvalue  n
# 1: 2000_1     AG 0.10118189 0.5059248 0.8993144 5.417640e-09 23
# 2: 2000_1   BECK 0.10123009 0.5054540 0.9041806 3.292129e-09 23
# 3: 2000_1 ELMORE 0.10123255 0.5054300 0.9048602 3.064459e-09 23
# 4: 2000_1  ZHANG 0.10122847 0.5054698 0.9042313 3.274661e-09 23
# 5: 2001_1     AG 0.08849003 0.5841794 0.9040939 3.322259e-09 23
# 6: 2001_1   BECK 0.08833724 0.5856141 0.9069664 2.445684e-09 23

# print(fit$fits$AG$`2002_1`$ws)
print(fit$`2002_1`$fFIT$AG$ws)
# $iter1
#  [1] 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000
#  [8] 0.2000000 0.2000000 0.5000000 1.0000000 0.9507082 0.9452160 0.1000000
# [15] 0.5000000 0.6913968 1.0000000 1.0000000 0.1000000 0.2000000 0.1000000
# [22] 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000
# [29] 0.2000000 0.2000000 0.2000000 1.0000000 0.8705899
# 
# $iter2
#  [1] 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000
#  [8] 0.2000000 0.2000000 0.5000000 1.0000000 0.9507082 0.9149648 0.2000000
# [15] 0.5000000 0.6790401 1.0000000 0.9928448 0.2000000 0.2000000 0.2000000
# [22] 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000
# [29] 0.2000000 0.2000000 0.2000000 1.0000000 0.8705899
## visualization
# svg("Figure1_phenofit_curve_fitting.svg", 11, 7)
# Cairo::CairoPDF(file_pdf, 11, 6) #
# dev.off()
g <- plot_phenofit(d_fit, brks2, titlestr)
grid::grid.newpage(); grid::grid.draw(g)# plot to check the curve fitting
```

<img src="man/Figure/curve fitting-1.svg" style="display: block; margin: auto;" />

2.5 Extract phenology
---------------------

``` r
# pheno: list(p_date, p_doy)
l_pheno <- get_pheno(fit, IsPlot = F) #%>% map(~melt_list(., "meth"))

# ratio = 1.15
# file <- "Figure5_Phenology_Extraction_temp.pdf"
# cairo_pdf(file, 8*ratio, 6*ratio)
# temp <- get_pheno(fit$fits$ELMORE[2:6], IsPlot = T)
# dev.off()
# file.show(file)

## check the extracted phenology
pheno <- get_pheno(fit[1:6], "ELMORE", IsPlot = T)
```

<img src="man/Figure/Extract phenology-1.svg" style="display: block; margin: auto;" />

``` r

# print(str(pheno, 1))
head(l_pheno$doy$AG)
#      flag     origin TRS2.sos TRS2.eos TRS5.sos TRS5.eos DER.sos DER.pop
# 1: 2000_1 2000-01-01      169      275      176      265     176     202
# 2: 2001_1 2001-01-01      147      263      156      255     155     199
# 3: 2002_1 2002-01-01      168      272      180      258     183     208
# 4: 2003_1 2003-01-01      134      271      150      253     154     182
# 5: 2004_1 2004-01-01      168      262      179      252     183     201
# 6: 2005_1 2005-01-01      146      266      157      254     155     228
#    DER.eos  UD  SD  DD  RD Greenup Maturity Senescence Dormancy
# 1:     267 165 187 250 281     138      175        269      309
# 2:     257 142 170 241 269     114      154        258      295
# 3:     259 163 198 239 279      72      186        260      326
# 4:     255 128 172 228 281      NA      163        257      354
# 5:     253 165 194 237 268      76      186        254      297
# 6:     250 140 174 235 272     109      153        245       NA
```
