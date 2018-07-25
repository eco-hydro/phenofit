
<!-- README.md is generated from README.Rmd. Please edit that file -->

## phenofit

[![Travis Build
Status](https://travis-ci.org/kongdd/phenofit.svg?branch=master)](https://travis-ci.org/kongdd/phenofit)
[![codecov](https://codecov.io/gh/kongdd/phenofit/branch/master/graph/badge.svg)](https://codecov.io/gh/kongdd/phenofit)

A state-of-the-art **remote sensing vegetation phenology** extraction
package: `phenofit`

  - `phenofit` combine merits of TIMESAT and phenopix
  - A simple and stable growing season dividing methods was proposed
  - Provide a practical snow elimination method, based on Whittaker
  - 7 curve fitting methods and 4 phenology extraction methods
  - We add parameters boundary for every curve fitting methods according
    to their ecological meaning.
  - `optimx` is used to select best optimization method for different
    curve fitting methods.

## Installation

You can install phenofit from github with:

``` r
# install.packages("devtools")
devtools::install_github("kongdd/phenofit")
```

## Example

Here, we illustrate how to use `phenofit` to extract vegetation
phenology from MOD13A1 in the sampled points. Regional analysis also can
be conducted in the similar way.

## 1.1 Download MOD13A1 data

Upload point shapefile into GEE, clip MOD13A1 and download vegetation
index data.
[Here](https://code.earthengine.google.com/ee3ec39cf3061374dab435c358d008a3)
is the corresponding GEE script.

## 1.2 Initial weights for input data

Load packages.

``` r
library(phenofit)
library(data.table)
library(magrittr)
library(lubridate)
library(purrr)
library(plyr)
```

Set global parameters for
`phenofit`

``` r
# lambda   <- 5    # non-parameter Whittaker, only suit for 16-day. Other time-scale
# should assign a lambda.
ymax_min   <- 0.1  # the maximum ymax shoud be greater than `ymax_min` 
rymin_less <- 0.8  # trough < ymin + A*rymin_less
nptperyear <- 23   # How many points for a single year
wFUN       <- wBisquare #wTSM #wBisquare # Weights updating function, could be one of `wTSM`, 'wBisquare', `wChen` and `wSELF`. 
```

Read the point shapefile to get points coordinate information. Read
Enhanced Vegetation Index (EVI) exported by `GEE`.

  - Add date according to composite day of the year (DayOfYear), other
    than image date.
  - Add weights according to `SummaryQA`.

For MOD13A1, Weights can by initialed by `SummaryQA` band (also suit for
MOD13A2 and MOD13Q1). We have written a qc function for `SummaryQA`,
`qc_summary`.

| SummaryQA        | Pixel reliability summary QA                        | weight |
| ---------------- | --------------------------------------------------- | ------ |
| \-1 Fill/No data | Not processed                                       | `wmin` |
| 0 Good data      | Use with confidence                                 | 1      |
| 1 Marginal data  | Useful but look at detailed QA for more information | 0.5    |
| 2 Snow/ice       | Pixel covered with snow/ice                         | `wmin` |
| 3 Cloudy         | Pixel is cloudy                                     | `wmin` |

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

Add one year in head and tail, for growing season dividing. For example,
the input data period is 20000218 ~ 20171219. After adding one year in
head and tail, it becomes 19990101 ~ 20181219.

## 2.1 load data

``` r
sites        <- unique(df$site)
sitename     <- sites[3]
d            <- df[site == sitename] # get the first site data
sp           <- st[site == sitename]

print      <- TRUE
IsPlot     <- TRUE # for brks

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

## 2.1 Check input data

``` r
dnew  <- add_HeadTail(d) # add additional one year in head and tail
INPUT <- check_input(dnew$t, dnew$y, dnew$w, maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
# y0 is used for plot. Original y value has been interpolated and changed.
INPUT$y0 <- dnew$y 
```

## 2.2 Divide growing seasons

Simply treating calendar year as a complete growing season will induce a
considerable error for phenology extraction. A simple growing season
dividing method was proposed in `phenofit`.

The growing season dividing method rely on heavily in Whittaker
smoother.

Procedures of initial weight, growing season dividing and curve fitting
are separated. Phenology extraction and curve fitting are combined
together.

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
brks2 <- season_3y(INPUT, nptperyear, south = sp$lat[1] < 0, 
                   FUN = whitsmw2, wFUN = wFUN,
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

<img src="man/Figure/divide growing season-1.svg" style="display: block; margin: auto;" />

## 2.3 Curve fitting

``` r
fit  <- curvefits(INPUT, brks2, lambda =lambda,
                  methods = c("AG", "zhang", "beck", "elmore"), #,"klos",, 'Gu'
                  nptperyear = nptperyear, debug = F, 
                  wFUN = wFUN,
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
#   [curvefits]  running 16 ...
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
#    flag      t0    mn    mx    rsp    a3    rau    a5
#    <fct>  <dbl> <dbl> <dbl>  <dbl> <dbl>  <dbl> <dbl>
#  1 2000_1  558. 0.176 0.408 0.0367  3.68 0.0150  6   
#  2 2001_1  925. 0.176 0.401 0.0202  5.05 0.0178  6   
#  3 2002_1 1291. 0.176 0.541 0.0368  2    0.0157  2   
#  4 2003_1 1636. 0.178 0.448 0.0270  2    0.0123  2.77
#  5 2004_1 2021. 0.182 0.464 0.0361  2    0.0202  5.69
#  6 2005_1 2413. 0.183 0.462 0.0131  6    0.0337  2   
#  7 2006_1 2739. 0.184 0.438 0.0213  2    0.0134  3.52
#  8 2007_1 3111. 0.186 0.483 0.0216  2    0.0149  2.72
#  9 2008_1 3484. 0.188 0.505 0.0224  2.08 0.0173  6   
# 10 2009_1 3885. 0.196 0.479 0.0140  6    0.0365  2   
# 11 2010_1 4197. 0.193 0.488 0.0266  2    0.0132  2   
# 12 2011_1 4566. 0.198 0.466 0.0340  2    0.0138  6   
# 13 2012_1 4926. 0.200 0.510 0.0388  4.24 0.0124  6   
# 14 2013_1 5287. 0.205 0.482 0.0402  6    0.0110  3.68
# 15 2014_1 5713. 0.208 0.494 0.0122  6    0.0455  6   
# 16 2015_1 6025. 0.219 0.512 0.0346  2.63 0.0160  6   
# 17 2016_1 6383. 0.214 0.485 0.0428  2    0.0130  4.99
# 18 2017_1 6753. 0.211 0.441 0.0284  6    0.0120  5.02

## Get GOF information
stat  <- ldply(fit$fits, function(fits_meth){
    ldply(fits_meth, statistic.phenofit, .id = "flag")
}, .id = "meth")
fit$stat <- stat
print(head(stat))
#   meth   flag       RMSE       NSE         R       pvalue  n
# 1   AG 2000_1 0.10229316 0.4913260 0.8964904 7.151335e-09 23
# 2   AG 2001_1 0.08966131 0.5730986 0.9040977 3.320912e-09 23
# 3   AG 2002_1 0.08394268 0.7169239 0.9143854 1.056728e-09 23
# 4   AG 2003_1 0.06395071 0.7639834 0.9488097 5.565019e-12 23
# 5   AG 2004_1 0.06179483 0.7424845 0.9292973 4.198802e-10 22
# 6   AG 2005_1 0.08001220 0.7216373 0.9301156 3.750332e-10 22

print(fit$fits$AG$`2002_1`$ws)
# $iter1
#  [1] 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000
#  [8] 0.2000000 0.9434168 0.6945310 0.4900482 0.9295393 1.0000000 1.0000000
# [15] 0.5000000 1.0000000 1.0000000 1.0000000 0.3578768 0.6543742 0.9554592
# [22] 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000
# 
# $iter2
#  [1] 0.9452294 0.9452294 0.9452294 0.9452294 0.9452294 0.9452294 0.9452277
#  [8] 0.9973005 0.9434168 0.9012773 0.9794146 0.9295393 0.2459014 1.0000000
# [15] 0.8185856 0.9768257 1.0000000 0.9966571 0.2000000 0.6543742 0.9554592
# [22] 0.9648546 0.9318852 0.9424703 0.9443088 0.9450985 0.9452044 0.9452287
## visualization
# svg("Figure1_phenofit_curve_fitting.svg", 11, 7)
# Cairo::CairoPDF(file_pdf, 11, 6) #
# dev.off()
g <- plot_phenofit(fit, d, titlestr)
# Warning: Removed 712 rows containing missing values (geom_point).
grid::grid.newpage(); grid::grid.draw(g)# plot to check the curve fitting
```

<img src="man/Figure/curve fitting-1.svg" style="display: block; margin: auto;" />

## 2.4 Extract phenology.

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

<img src="man/Figure/Extract phenology-1.svg" style="display: block; margin: auto;" />

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
#   flag  origin     TRS1.sos TRS1.eos TRS2.sos TRS2.eos TRS5.sos TRS5.eos
#   <fct> <date>        <dbl>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl>
# 1 2000… 2000-01-01      166      278      169      272      176      263
# 2 2001… 2001-01-01      142      266      146      262      155      254
# 3 2002… 2002-01-01      159      301      167      283      179      255
# 4 2003… 2003-01-01      125      294      135      279      152      254
# 5 2004… 2004-01-01      159      259      166      255      179      248
# 6 2005… 2005-01-01      140      274      145      265      156      253
# # ... with 13 more variables: TRS6.sos <dbl>, TRS6.eos <dbl>,
# #   DER.sos <dbl>, DER.pop <dbl>, DER.eos <dbl>, GU.UD <dbl>, GU.SD <dbl>,
# #   GU.DD <dbl>, GU.RD <dbl>, ZHANG.Greenup <dbl>, ZHANG.Maturity <dbl>,
# #   ZHANG.Senescence <dbl>, ZHANG.Dormancy <dbl>
```
