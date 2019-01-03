
<!-- README.md is generated from README.Rmd. Please edit that file -->

# phenofit

[![Travis Build Status](https://travis-ci.org/kongdd/phenofit.svg?branch=master)](https://travis-ci.org/kongdd/phenofit)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/kongdd/phenofit?branch=master&svg=true)](https://ci.appveyor.com/project/kongdd/phenofit)
[![codecov](https://codecov.io/gh/kongdd/phenofit/branch/master/graph/badge.svg)](https://codecov.io/gh/kongdd/phenofit)
[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
[![CRAN](http://www.r-pkg.org/badges/version/phenofit)](https://cran.r-project.org/package=phenofit)

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

![title](man/Figure/Figure1_phenofit_flowchart.svg)  
*<u>Figure 1. The flowchart of phenology extraction in `phenofit`.</u>*

# Installation

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

![](man/Figure/phenofit_shiny.png)  
*<u>Figure 2. Shiny interface of `phenofit`.</u>*

# R code Example

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
suppressMessages({
    library(data.table)
    library(magrittr)
    library(lubridate)
    library(purrr)
    library(plyr)
    
    library(phenofit)
})
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

Add one year in head and tail, for growing season dividing. For example,
the input data period is 20000218 ~ 20171219. After adding one year in
head and tail, it becomes 19990101 ~ 20181219.

## 2.1 load data

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

If need night temperature (Tn) to constrain ungrowing season backgroud
value, NA values in Tn should be filled.

``` r
d$Tn %<>% zoo::na.approx(maxgap = 4)
plot(d$Tn, type = "l"); abline(a = 5, b = 0, col = "red")
```

## 2.2 Check input data

``` r
dnew  <- add_HeadTail(d, south, nptperyear = 23) # add additional one year in head and tail
INPUT <- check_input(dnew$t, dnew$y, dnew$w, 
                     nptperyear, south, 
                     maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
```

## 2.3 Divide growing seasons

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
#                FUN = wWHIT, wFUN = wFUN, iters = 2,
#                lambda = lambda,
#                IsPlot = IsPlot, plotdat = d,
#                south = d$lat[1] < 0,
#                rymin_less = 0.6, ymax_min = ymax_min,
#                max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5) #, ...
# get growing season breaks in a 3-year moving window
brks2 <- season_3y(INPUT, 
                   FUN = wWHIT, wFUN = wFUN,
                   maxExtendMonth = 6, threshold_min = 0.1,
                   IsPlot = IsPlot, IsPlot.OnlyBad = FALSE, print = print)
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

## 2.4 Curve fitting

``` r
fit  <- curvefits(INPUT, brks2,
                  methods = c("AG", "Zhang", "Beck", "Elmore"), #,"klos",, 'Gu'
                  verbose = FALSE, 
                  wFUN = wFUN,
                  nextent = 2, maxExtendMonth = 3, minExtendMonth = 1,
                  QC_flag = dnew$QC_flag, minPercValid = 0.2,
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

## check the curve fitting parameters
params <- get_param(fit)
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
#  1 2000_1  223. 0.158 0.410 0.0189  5.50 0.0220  3.06
#  2 2001_1  565. 0.171 0.402 0.0208  4.76 0.0171  6   
#  3 2002_1  939. 0.186 0.501 0.0296  2.27 0.0181  2.85
#  4 2003_1 1278. 0.162 0.447 0.0251  2    0.0128  3.46
#  5 2004_1 1661. 0.176 0.449 0.0392  2    0.0178  4.55
#  6 2005_1 2054. 0.175 0.461 0.0131  6    0.0319  2   
#  7 2006_1 2383. 0.177 0.435 0.0201  2.29 0.0138  3.88
#  8 2007_1 2755. 0.156 0.480 0.0189  2    0.0146  2.72
#  9 2008_1 3117. 0.169 0.502 0.0265  2    0.0148  3.57
# 10 2009_1 3527. 0.181 0.475 0.0133  5.92 0.0345  2   
# 11 2010_1 3840. 0.174 0.482 0.0238  2    0.0121  2   
# 12 2011_1 4208. 0.180 0.465 0.0301  2    0.0131  4.97
# 13 2012_1 4578. 0.160 0.509 0.0242  3.91 0.0136  4.98
# 14 2013_1 4929. 0.188 0.482 0.0386  6    0.0112  4.04
# 15 2014_1 5306. 0.201 0.501 0.0341  2    0.0137  3.51
# 16 2015_1 5692. 0.211 0.503 0.0178  6.00 0.0289  4.44
# 17 2016_1 6024. 0.192 0.484 0.0398  2    0.0126  4.78
# 18 2017_1 6404. 0.170 0.445 0.0213  5.88 0.0123  3.71

df_fit <- get_fitting(fit)
## Get GOF information
stat <- ldply(fit, GOF_fFITs, .id = "flag") %>% data.table()
# fit$stat <- stat
print(head(stat))
#      flag   meth       RMSE       NSE         R       pvalue  n
# 1: 2000_1     AG 0.10118884 0.5058570 0.9009524 4.594461e-09 23
# 2: 2000_1   BECK 0.10118971 0.5058484 0.9029269 3.752268e-09 23
# 3: 2000_1 ELMORE 0.10122565 0.5054974 0.9008057 4.663310e-09 23
# 4: 2000_1  ZHANG 0.10118802 0.5058650 0.9029436 3.745764e-09 23
# 5: 2001_1     AG 0.09080665 0.5421617 0.9087420 1.372862e-10 26
# 6: 2001_1   BECK 0.09110927 0.5391050 0.9142823 6.666951e-11 26

# print(fit$fits$AG$`2002_1`$ws)
print(fit$`2002_1`$fFIT$AG$ws)
# $iter1
#  [1] 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.5000000
#  [8] 1.0000000 0.9192712 0.9979039 0.1000000 0.5000000 0.8826912 1.0000000
# [15] 1.0000000 0.1000000 0.2000000 0.1000000 0.2000000 0.2000000 0.2000000
# [22] 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000
# [29] 1.0000000 0.7849687
# 
# $iter2
#  [1] 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.5000000
#  [8] 1.0000000 0.9192712 0.9739109 0.2000000 0.5000000 0.8736643 1.0000000
# [15] 0.9935016 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000
# [22] 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000
# [29] 1.0000000 0.7849687
## visualization
# svg("Figure1_phenofit_curve_fitting.svg", 11, 7)
# Cairo::CairoPDF(file_pdf, 11, 6) #
# dev.off()
g <- plot_phenofit(df_fit, brks2, titlestr)
grid::grid.newpage(); grid::grid.draw(g)# plot to check the curve fitting
```

<img src="man/Figure/curve fitting-1.svg" style="display: block; margin: auto;" />

## 2.5 Extract phenology

``` r
# pheno: list(p_date, p_doy)
p <- PhenoExtract(fit, IsPlot = F)
pheno <- map(p, tidyFitPheno) %>% purrr::transpose() %>% 
    map(~melt_list(., "meth"))
fit$pheno <- pheno

# ratio = 1.15
# file <- "Figure5_Phenology_Extraction_temp.pdf"
# cairo_pdf(file, 8*ratio, 6*ratio)
# temp <- PhenoExtract(fit$fits$ELMORE[2:6], IsPlot = T)
# dev.off()
# file.show(file)

## check the extracted phenology
p <- PhenoExtract(fit[1:6], "ELMORE", IsPlot = T)
```

<img src="man/Figure/Extract phenology-1.svg" style="display: block; margin: auto;" />

``` r

print(str(pheno, 2))
# List of 2
#  $ date:Classes 'data.table' and 'data.frame':    72 obs. of  18 variables:
#   ..$ flag      : Factor w/ 18 levels "2000_1","2001_1",..: 1 2 3 4 5 6 7 8 9 10 ...
#   ..$ origin    : Date[1:72], format: "2000-01-01" ...
#   ..$ TRS2.sos  : Date[1:72], format: "2000-06-15" ...
#   ..$ TRS2.eos  : Date[1:72], format: "2000-10-01" ...
#   ..$ TRS5.sos  : Date[1:72], format: "2000-06-23" ...
#   ..$ TRS5.eos  : Date[1:72], format: "2000-09-19" ...
#   ..$ DER.sos   : Date[1:72], format: "2000-06-21" ...
#   ..$ DER.pop   : Date[1:72], format: "2000-08-10" ...
#   ..$ DER.eos   : Date[1:72], format: "2000-09-20" ...
#   ..$ UD        : Date[1:72], format: "2000-06-11" ...
#   ..$ SD        : Date[1:72], format: "2000-07-05" ...
#   ..$ DD        : Date[1:72], format: "2000-09-01" ...
#   ..$ RD        : Date[1:72], format: "2000-10-07" ...
#   ..$ Greenup   : Date[1:72], format: "2000-05-14" ...
#   ..$ Maturity  : Date[1:72], format: "2000-06-20" ...
#   ..$ Senescence: Date[1:72], format: "2000-09-20" ...
#   ..$ Dormancy  : Date[1:72], format: "2000-11-26" ...
#   ..$ meth      : chr [1:72] "AG" "AG" "AG" "AG" ...
#   ..- attr(*, ".internal.selfref")=<externalptr> 
#  $ doy :Classes 'data.table' and 'data.frame':    72 obs. of  18 variables:
#   ..$ flag      : Factor w/ 18 levels "2000_1","2001_1",..: 1 2 3 4 5 6 7 8 9 10 ...
#   ..$ origin    : Date[1:72], format: "2000-01-01" ...
#   ..$ TRS2.sos  : num [1:72] 167 147 168 135 168 146 131 138 151 159 ...
#   ..$ TRS2.eos  : num [1:72] 275 263 275 271 262 268 274 275 271 277 ...
#   ..$ TRS5.sos  : num [1:72] 175 156 180 151 179 156 150 158 166 169 ...
#   ..$ TRS5.eos  : num [1:72] 263 255 258 253 252 254 258 256 256 264 ...
#   ..$ DER.sos   : num [1:72] 173 154 183 155 182 154 153 162 170 167 ...
#   ..$ DER.pop   : num [1:72] 223 201 210 184 201 228 193 200 197 240 ...
#   ..$ DER.eos   : num [1:72] 264 257 257 255 254 250 259 257 258 261 ...
#   ..$ UD        : num [1:72] 163 142 163 129 165 140 124 130 145 153 ...
#   ..$ SD        : num [1:72] 187 169 199 174 194 173 176 186 187 186 ...
#   ..$ DD        : num [1:72] 245 242 235 228 237 235 235 230 235 247 ...
#   ..$ RD        : num [1:72] 281 269 283 281 269 273 283 285 280 282 ...
#   ..$ Greenup   : num [1:72] 135 114 NA NA 75 108 NA NA NA 121 ...
#   ..$ Maturity  : num [1:72] 172 154 186 164 186 152 157 174 178 165 ...
#   ..$ Senescence: num [1:72] 264 258 256 256 255 245 261 256 259 256 ...
#   ..$ Dormancy  : num [1:72] 331 295 387 391 301 384 338 NA 343 390 ...
#   ..$ meth      : chr [1:72] "AG" "AG" "AG" "AG" ...
#   ..- attr(*, ".internal.selfref")=<externalptr> 
# NULL
head(pheno$doy$AG)
# NULL
```
