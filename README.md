
# phenofit

[![Travis Build
Status](https://travis-ci.org/kongdd/phenofit.svg?branch=master)](https://travis-ci.org/kongdd/phenofit)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/kongdd/phenofit?branch=master&svg=true)](https://ci.appveyor.com/project/kongdd/phenofit)
[![codecov](https://codecov.io/gh/kongdd/phenofit/branch/master/graph/badge.svg)](https://codecov.io/gh/kongdd/phenofit)
[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
[![CRAN](http://www.r-pkg.org/badges/version/phenofit)](https://cran.r-project.org/package=phenofit)
[![total](http://cranlogs.r-pkg.org/badges/grand-total/phenofit)](https://www.rpackages.io/package/phenofit)
[![monthly](http://cranlogs.r-pkg.org/badges/phenofit)](https://www.rpackages.io/package/phenofit)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3605560.svg)](https://doi.org/10.5281/zenodo.3605560)

A state-of-the-art **remote sensing vegetation phenology** extraction
package: `phenofit`

-   `phenofit` combine merits of TIMESAT and phenopix
-   A simple and stable growing season dividing methods was proposed
-   Provide a practical snow elimination method, based on Whittaker
-   7 curve fitting methods and 4 phenology extraction methods
-   We add parameters boundary for every curve fitting methods according
    to their ecological meaning.
-   `optimx` is used to select best optimization method for different
    curve fitting methods.

***Task lists***

-   [ ] Test the performance of `phenofit` in multiple growing season
    regions (e.g. the North China Plain);
-   [ ] Uncertainty analysis of curve fitting and phenological metrics;
-   [x] shiny app has been moved to
    [phenofit.shiny](https://github.com/kongdd/phenofit.shiny);
-   [x] Complete script automatic generating module in shinyapp;
-   [x] `Rcpp` improve double logistics optimization efficiency by 60%;
-   [x] Support spatial analysis;
-   [x] Support annual season in curve fitting;
-   [x] flexible fine fitting input ( original time-series or smoothed
    time-series by rough fitting).
-   [x] Asymmetric of Threshold method

![title](man/Figure/Figure1_phenofit_flowchart.svg)

*<u>Figure 1. The flowchart of phenology extraction in `phenofit`.</u>*

# Installation

You can install phenofit from github with:

``` r
# install.packages("devtools")
devtools::install_github("kongdd/phenofit")
```

# Example

Here, we illustrate how to use `phenofit` to extract vegetation
phenology from MOD13A1 in the sampled points. Regional analysis also can
be conducted in the similar way.

<!-- ## 1.1 Download MOD13A1 data

Upload point shapefile into GEE, clip MOD13A1 and download vegetation index
data. [Here](https://code.earthengine.google.com/ee3ec39cf3061374dab435c358d008a3) is the corresponding GEE script. 
 -->

## 1.1 Initial weights for input data

Load packages.

``` r
suppressMessages({
    library(data.table)
    library(magrittr)
    library(lubridate)
    library(purrr)
    library(plyr)
    library(ggplot2)
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

-   Add date according to composite day of the year (DayOfYear), other
    than image date.
-   Add weights according to `SummaryQA`.

For MOD13A1, Weights can by initialed by `SummaryQA` band (also suit for
MOD13A2 and MOD13Q1). There is already a `QC` function for `SummaryQA`,
i.e. `qc_summary`.

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
df <- df[, .(site, y = EVI/1e4, t, date, w, QC_flag)]
```

-   `add_HeadTail` is used to deal with such situation, e.g. MOD13A2
    begins from 2000-02-08. We need to construct a series with complete
    year, which begins from 01-01 for NH, and 07-01 for SH. For example,
    the input data period is 20000218 \~ 20171219. After adding one year
    in head and tail, it becomes 19990101 \~ 20181219.

## 2.1 load site data

``` r
sites        <- unique(df$site)
sitename     <- sites[3]
d            <- df[site == sitename] # get the first site data
sp           <- st[site == sitename]

south      <- sp$lat < 0
print      <- FALSE # whether print progress
IsPlot     <- TRUE  # for brks

prefix_fig <- "phenofit"
titlestr   <- with(sp, sprintf('[%03d,%s] %s, lat = %5.2f, lon = %6.2f',
                                     ID, site, IGBPname, lat, lon))
file_pdf   <- sprintf('Figure/%s_[%03d]_%s.pdf', prefix_fig, sp$ID[1], sp$site[1])
```

If need night temperature (Tn) to constrain ungrowing season backgroud
value, NA values in Tn should be filled.

``` r
d$Tn %<>% zoo::na.approx(maxgap = 4)
plot(d$Tn, type = "l"); abline(a = 5, b = 0, col = "red")
```

## 2.2 Check input data

``` r
dnew  <- add_HeadTail(d, south, nptperyear = 23) # add additional one year in head and tail
INPUT <- check_input(dnew$t, dnew$y, dnew$w, dnew$QC_flag,
                     nptperyear, south, 
                     maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
```

## 2.3 Divide growing seasons

Simply treating calendar year as a complete growing season will induce a
considerable error for phenology extraction. A simple growing season
dividing method was proposed in `phenofit`.

The growing season dividing method rely on heavily in Whittaker
smoother.

Procedures of initial weight, growing season dividing, curve fitting,
and phenology extraction are conducted separately.

``` r
par(mar = c(3, 2, 2, 1), mgp = c(3, 0.6, 0))
lambda <- init_lambda(INPUT$y)
# The detailed information of those parameters can be seen in `season`.
# brks   <- season(INPUT, nptperyear,
#                FUN = smooth_wWHIT, wFUN = wFUN, iters = 2,
#                lambda = lambda,
#                IsPlot = IsPlot, plotdat = d,
#                south = d$lat[1] < 0,
#                rymin_less = 0.6, ymax_min = ymax_min,
#                max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5) #, ...
# get growing season breaks in a 3-year moving window
brks2 <- season_mov(INPUT, 
                   FUN = smooth_wWHIT, wFUN = wFUN,
                   maxExtendMonth = 6, r_min = 0.1,
                   IsPlot = IsPlot, IsPlot.OnlyBad = FALSE, print = print)
```

<img src="man/Figure/divide growing season-1.svg" style="display: block; margin: auto;" />

## 2.4 Curve fitting

``` r
fit  <- curvefits(INPUT, brks2,
                  methods = c("AG", "Zhang", "Beck", "Elmore"), #,"klos",, 'Gu'
                  wFUN = wFUN,
                  nextend = 2, maxExtendMonth = 3, minExtendMonth = 1, minPercValid = 0.2,
                  print = print, verbose = FALSE)

## check the curve fitting parameters
l_param <- get_param(fit)
print(str(l_param, 1))
# List of 4
#  $ AG    : tibble [18 x 8] (S3: tbl_df/tbl/data.frame)
#  $ Beck  : tibble [18 x 7] (S3: tbl_df/tbl/data.frame)
#  $ Elmore: tibble [18 x 8] (S3: tbl_df/tbl/data.frame)
#  $ Zhang : tibble [18 x 8] (S3: tbl_df/tbl/data.frame)
# NULL
print(l_param$AG)
# # A tibble: 18 x 8
#    flag      t0    mn    mx    rsp    a3     rau    a5
#    <fct>  <dbl> <dbl> <dbl>  <dbl> <dbl>   <dbl> <dbl>
#  1 2000_1  199. 0.167 0.407 0.0351  3.19 0.0146   6   
#  2 2001_1  579. 0.169 0.397 0.0157  6    0.0227   4.57
#  3 2002_1  934. 0.167 0.528 0.0329  2    0.0161   2   
#  4 2003_1 1275. 0.167 0.432 0.0268  2    0.0123   3.73
#  5 2004_1 1689. 0.169 0.445 0.0166  5.08 0.0372   2   
#  6 2005_1 2052. 0.172 0.470 0.0141  6    0.0301   2   
#  7 2006_1 2370. 0.166 0.420 0.0280  2.02 0.0112   3.12
#  8 2007_1 2751. 0.165 0.479 0.0212  2    0.0142   2.98
#  9 2008_1 3123. 0.169 0.481 0.0215  2.70 0.0164   6   
# 10 2009_1 3522. 0.168 0.479 0.0142  6    0.0278   2   
# 11 2010_1 3840. 0.167 0.489 0.0229  2    0.0131   2   
# 12 2011_1 4208. 0.175 0.468 0.0288  2    0.0133   5.00
# 13 2012_1 4558. 0.166 0.505 0.0478  2    0.0109   4.16
# 14 2013_1 4968. 0.166 0.484 0.0137  6    0.0204   2.26
# 15 2014_1 5328. 0.168 0.501 0.0166  3.90 0.0164   2.35
# 16 2015_1 5701. 0.189 0.484 0.0146  6    0.0287   2.03
# 17 2016_1 6027. 0.172 0.472 0.0299  2    0.0130   5.57
# 18 2017_1 6381. 0.168 0.441 0.0403  2.95 0.00979  5.98

d_fit <- get_fitting(fit)
## Get GOF information
d_gof <- get_GOF(fit)
# fit$stat <- stat
print(head(d_gof))
#      flag   meth       RMSE       NSE         R       pvalue  n
# 1: 2000_1     AG 0.09691142 0.4963370 0.9104446 4.481677e-11 27
# 2: 2000_1   Beck 0.09688513 0.4966102 0.9127208 3.289646e-11 27
# 3: 2000_1 Elmore 0.09690657 0.4963874 0.9109116 4.209014e-11 27
# 4: 2000_1  Zhang 0.09688479 0.4966137 0.9127908 3.258098e-11 27
# 5: 2001_1     AG 0.08197474 0.6537794 0.9011935 5.927392e-08 20
# 6: 2001_1   Beck 0.08139990 0.6586180 0.9047362 4.324325e-08 20

# print(fit$fits$AG$`2002_1`$ws)
print(fit$`2002_1`$fFIT$AG$ws)
# $iter1
#  [1] 0.8 0.2 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.5 1.0 1.0 1.0 1.0 0.5 1.0
# [20] 1.0 1.0 1.0 0.2 0.2 0.2 0.2 0.2 0.8 0.8 0.8 0.8 0.8
# 
# $iter2
#  [1] 0.8000000 0.2000000 0.8000000 0.8000000 0.8000000 0.8000000 0.8000000
#  [8] 0.8000000 0.8000000 0.8000000 0.8000000 0.8000000 0.5000000 1.0000000
# [15] 1.0000000 0.2000000 0.2000000 0.2873749 0.2000000 0.2000000 1.0000000
# [22] 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.8000000
# [29] 0.8000000 0.8000000 0.8000000 0.8000000
## visualization
g <- plot_phenofit(d_fit, brks2, NULL, title.ylab = "NDVI", "Time",
                   theme = coord_cartesian(xlim = c(ymd("2000-04-01"), ymd("2017-07-31"))))
# Coordinate system already present. Adding new coordinate system, which will replace the existing one.
grid::grid.newpage(); grid::grid.draw(g)# plot to check the curve fitting
```

<img src="man/Figure/curve fitting-1.svg" style="display: block; margin: auto;" />

``` r
# write_fig(g, "Figure1_phenofit_curve_fitting.pdf", 10, 6)
```

## 2.5 Extract phenology

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
pheno <- get_pheno(fit[1:6], "Elmore", IsPlot = T)
```

<img src="man/Figure/Extract phenology-1.svg" style="display: block; margin: auto;" />

``` r
# print(str(pheno, 1))
head(l_pheno$doy$AG)
#      flag     origin TRS2.sos TRS2.eos TRS5.sos TRS5.eos TRS6.sos TRS6.eos
# 1: 2000_1 2000-01-01      167      275      175      264      177      261
# 2: 2001_1 2001-01-01      145      263      154      255      157      252
# 3: 2002_1 2002-01-01      164      284      178      256      182      248
# 4: 2003_1 2003-01-01      132      273      149      253      153      247
# 5: 2004_1 2004-01-01      162      264      173      252      176      248
# 6: 2005_1 2005-01-01      148      268      159      254      162      250
#    DER.sos DER.pop DER.eos  UD  SD  DD  RD Greenup Maturity Senescence Dormancy
# 1:     175     200     266 163 186 249 281     157      194        243      287
# 2:     152     214     256 140 168 241 268     133      174        235      274
# 3:     182     204     248 161 196 222 292     153      203        219      307
# 4:     153     180     254 127 170 228 282     118      179        197      294
# 5:     171     229     248 157 189 235 268     150      196        274       NA
# 6:     157     225     249 143 175 233 273     135      181        281       NA
```

# **References**

> \[1\] Kong, D., Zhang, Y., Wang, D., Chen, J., & Gu, X. (2020).
> Photoperiod Explains the Asynchronization Between Vegetation Carbon
> Phenology and Vegetation Greenness Phenology. *Journal of Geophysical
> Research: Biogeosciences*, 125(8), e2020JG005636.
> <https://doi.org/10.1029/2020JG005636>
>
> \[2\] Kong, D., Zhang, Y., Gu, X., & Wang, D. (2019). A robust method
> for reconstructing global MODIS EVI time series on the Google Earth
> Engine. *ISPRS Journal of Photogrammetry and Remote Sensing*, 155,
> 13–24.
>
> \[3\] Kong, D., (2020). R package: A state-of-the-art Vegetation
> Phenology extraction package, `phenofit` version 0.2.6,
> <https://doi.org/10.5281/zenodo.3605560>
>
> \[4\] Zhang, Q., Kong, D., Shi, P., Singh, V.P., Sun, P., 2018.
> Vegetation phenology on the Qinghai-Tibetan Plateau and its response
> to climate change (1982–2013). Agric. For. Meteorol. 248, 408–417.
> <https://doi.org/10.1016/j.agrformet.2017.10.026>

# Acknowledgements

Keep in mind that this repository is released under a GPL2 license,
which permits commercial use but requires that the source code (of
derivatives) is always open even if hosted as a web service.
