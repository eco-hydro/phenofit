
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

***Task lists***

  - [ ] Test the performance of `phenofit` in multiple growing season regions (e.g. the North China Plain);
  - [ ] separate shiny application into a independent repository;
  - [ ] Improve computational efficiency of fine fitting;
  - [x] Complete script automatic generating module in shinyapp;
  - [ ] Uncertainty analysis of curve fitting and phenological metrics;
  - [x] Support spatial analysis;
  - [x] Support annual season in curve fitting;
  - [x] flexible fine fitting input ( original time-series or smoothed
    time-series by rough fitting).
  - [x] Asymmetric of Threshold method

![title](man/Figure/Figure1_phenofit_flowchart.svg)

*<u>Figure 1. The flowchart of phenology extraction in `phenofit`.</u>*

# Installation

You can install phenofit from github with:

``` r
# install.packages("devtools")
devtools::install_github("kongdd/phenofit")
```

Run `shinyapp`:

``` r
shiny::runGitHub("phenofit", "kongdd", subdir = "inst/shiny/phenofit")
# Or run locally
shiny::runApp(system.file("shiny/phenofit", package = "phenofit"))
```

Running in docker: 

```bash
docker run -it -p 8787:8787 `
    -v E:/ubuntu:/ubuntu `
    -v E:/ubuntu/site-library:/usr/local/lib/R/site-library `
    -v E:/ubuntu/etc:/usr/local/lib/R/etc `
    -v E:/github:/github 
    --name R-3.5.3 kongdd/phenofit bash 

# powershell docker rm @(docker ps -aq)
# powershell docker rmi @(docker images -f "dangling=true" -q)
```

<!-- ![](man/Figure/phenofit_shiny.png)
*<u>Figure 2. Shiny interface of `phenofit`.</u>*
 -->

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

  - Add date according to composite day of the year (DayOfYear), other
    than image date.
  - Add weights according to `SummaryQA`.

For MOD13A1, Weights can by initialed by `SummaryQA` band (also suit for
MOD13A2 and MOD13Q1). There is already a `QC` function for `SummaryQA`,
i.e. `qc_summary`.

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
df <- df[, .(site, y = EVI/1e4, t, date, w, QC_flag)]
```

  - `add_HeadTail` is used to deal with such situation, e.g. MOD13A2
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
#  $ AG    :Classes 'tbl_df', 'tbl' and 'data.frame':   18 obs. of  8 variables:
#  $ Beck  :Classes 'tbl_df', 'tbl' and 'data.frame':   18 obs. of  7 variables:
#  $ Elmore:Classes 'tbl_df', 'tbl' and 'data.frame':   18 obs. of  8 variables:
#  $ Zhang :Classes 'tbl_df', 'tbl' and 'data.frame':   18 obs. of  8 variables:
# NULL
print(l_param$AG)
# # A tibble: 18 x 8
#    flag      t0    mn    mx    rsp    a3    rau    a5
#    <fct>  <dbl> <dbl> <dbl>  <dbl> <dbl>  <dbl> <dbl>
#  1 2000_1  201. 0.168 0.407 0.0323  3.53 0.0151  5.83
#  2 2001_1  561. 0.173 0.407 0.0225  4.06 0.0162  6   
#  3 2002_1  931. 0.188 0.511 0.0383  2    0.0169  4.21
#  4 2003_1 1275. 0.167 0.434 0.0268  2    0.0122  3.37
#  5 2004_1 1660. 0.175 0.451 0.0363  2    0.0179  4.46
#  6 2005_1 2052. 0.180 0.466 0.0141  6    0.0314  2   
#  7 2006_1 2379. 0.174 0.436 0.0226  2.51 0.0131  3.26
#  8 2007_1 2752. 0.165 0.483 0.0208  2    0.0150  2.88
#  9 2008_1 3134. 0.177 0.492 0.0180  3.50 0.0199  6   
# 10 2009_1 3525. 0.172 0.480 0.0133  5.40 0.0313  2   
# 11 2010_1 3838. 0.194 0.487 0.0269  2    0.0146  2.39
# 12 2011_1 4206. 0.189 0.464 0.0327  2    0.0135  6   
# 13 2012_1 4558. 0.166 0.512 0.0473  2    0.0109  3.56
# 14 2013_1 4966. 0.168 0.484 0.0140  6    0.0190  2   
# 15 2014_1 5350. 0.204 0.479 0.0127  6    0.0376  6   
# 16 2015_1 5690. 0.215 0.494 0.0182  6    0.0275  4.81
# 17 2016_1 6023. 0.193 0.484 0.0428  2    0.0126  4.93
# 18 2017_1 6405. 0.171 0.447 0.0204  6    0.0129  3.65

d_fit <- get_fitting(fit)
## Get GOF information
d_gof <- get_GOF(fit)
# fit$stat <- stat
print(head(d_gof))
#      flag   meth       RMSE       NSE         R       pvalue  n
# 1: 2000_1     AG 0.10058989 0.4952465 0.9016506 1.809700e-09 24
# 2: 2000_1   Beck 0.10059953 0.4951497 0.9036266 1.461349e-09 24
# 3: 2000_1 Elmore 0.10060466 0.4950983 0.9021757 1.710497e-09 24
# 4: 2000_1  Zhang 0.10059866 0.4951585 0.9036627 1.455586e-09 24
# 5: 2001_1     AG 0.08808864 0.5879431 0.9037624 3.439665e-09 23
# 6: 2001_1   Beck 0.08817911 0.5870963 0.9047721 3.093163e-09 23

# print(fit$fits$AG$`2002_1`$ws)
print(fit$`2002_1`$fFIT$AG$ws)
# $iter1
#  [1] 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.5 1.0 1.0 1.0 1.0 0.5 1.0 1.0
# [18] 1.0 1.0 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 1.0 1.0
# 
# $iter2
#  [1] 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000
#  [8] 0.2000000 0.2000000 0.5000000 1.0000000 1.0000000 0.7893075 1.0000000
# [15] 0.4984446 0.8873205 1.0000000 1.0000000 1.0000000 0.2000000 0.2000000
# [22] 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000 0.2000000
# [29] 0.2000000 0.2000000 0.2000000 1.0000000 1.0000000
## visualization
# svg("Figure1_phenofit_curve_fitting.svg", 11, 7)
# Cairo::CairoPDF(file_pdf, 11, 6) #
# dev.off()
g <- plot_phenofit(d_fit, brks2, titlestr)
grid::grid.newpage(); grid::grid.draw(g)# plot to check the curve fitting
```

<img src="man/Figure/curve fitting-1.svg" style="display: block; margin: auto;" />

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
#      flag     origin TRS2.sos TRS2.eos TRS5.sos TRS5.eos DER.sos DER.pop
# 1: 2000_1 2000-01-01      167      273      175      264     174     203
# 2: 2001_1 2001-01-01      145      263      155      254     154     196
# 3: 2002_1 2002-01-01      168      268      179      256     183     202
# 4: 2003_1 2003-01-01      133      273      149      253     153     180
# 5: 2004_1 2004-01-01      166      262      178      251     181     201
# 6: 2005_1 2005-01-01      150      266      159      252     157     225
#    DER.eos  UD  SD  DD  RD Greenup Maturity Senescence Dormancy
# 1:     266 164 186 249 280     157      193        243      287
# 2:     256 141 170 240 268     133      177        234      275
# 3:     257 164 195 238 274     157      201        223      282
# 4:     254 128 170 225 283     118      179        196      298
# 5:     253 162 193 236 268     154      200        226      276
# 6:     248 143 175 233 272     136      181        279       NA
```

# **References**

> \[1\] Dongdong Kong, R package: A state-of-the-art Vegetation Phenology extraction package, `phenofit` version 0.2.2, <https://github.com/kongdd/phenofit>
> 
> \[2\] Zhang, Q., Kong, D., Shi, P., Singh, V.P., Sun, P., 2018. Vegetation phenology on the Qinghai-Tibetan Plateau and its response to climate change (1982–2013). Agric. For. Meteorol. 248, 408–417. <https://doi.org/10.1016/j.agrformet.2017.10.026>

# Acknowledgements

Keep in mind that this repository is released under a GPL2 license,
which permits commercial use but requires that the source code (of
derivatives) is always open even if hosted as a web service.