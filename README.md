
# phenofit

[![R-CMD-check](https://github.com/eco-hydro/phenofit/workflows/R-CMD-check/badge.svg)](https://github.com/eco-hydro/phenofit/actions)
[![codecov](https://codecov.io/gh/eco-hydro/phenofit/branch/master/graph/badge.svg)](https://app.codecov.io/gh/eco-hydro/phenofit)
[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
[![CRAN](http://www.r-pkg.org/badges/version/phenofit)](https://cran.r-project.org/package=phenofit)
[![total](http://cranlogs.r-pkg.org/badges/grand-total/phenofit)](https://www.rpackages.io/package/phenofit)
[![monthly](http://cranlogs.r-pkg.org/badges/phenofit)](https://www.rpackages.io/package/phenofit)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6320537.svg)](https://doi.org/10.5281/zenodo.6320537)

A state-of-the-art **remote sensing vegetation phenology** extraction
package: `phenofit`

-   `phenofit` combine merits of TIMESAT and phenopix
-   A simple and stable growing season dividing method was proposed
-   Provide a practical snow elimination method based on Whittaker
-   7 curve fitting methods and 4 phenology extraction methods
-   We add parameters boundary for every curve fitting method according
    to their ecological meaning.
-   `optimx` is used to select the best optimization method for different
    curve fitting methods.

***Task lists***

-   [x] Test the performance of `phenofit` in multiple growing seasons
    regions (e.g., the North China Plain);
-   [ ] Uncertainty analysis of curve fitting and phenological metrics;
-   [x] shiny app has been moved to
    [phenofit.shiny](https://github.com/eco-hydro/phenofit.shiny);
-   [x] Complete script automatic generating module in shinyapp;
-   [x] `Rcpp` improve double logistics optimization efficiency by 60%;
-   [x] Support spatial analysis;
-   [x] Support annual season in curve fitting;
-   [x] flexible fine fitting input ( original time-series or smoothed
    time-series by rough fitting).
-   [x] Asymmetric Threshold method

<!-- ![title](man/Figure/Figure1_phenofit_flowchart.svg)

*<u>Figure 1. The flowchart of phenology extraction in `phenofit`.</u>* -->

# Installation

You can install phenofit from github with:

``` r
# install.packages("remotes")
remotes::install_github("eco-hydro/phenofit")
```

# Note

Users can through the following options to improve the performance of phenofit in multiple growing 
season regions:

- Users can decrease those three parameters `nextend`, `minExtendMonth` and
  `maxExtendMonth` to a relative low value, by setting option 
  `set_options(fitting = list(nextend = 1, minExtendMonth = 0, maxExtendMonth = 0.5))`.

- Use `wHANTS` as the rough fitting function. Due to the nature of Fourier
  functions, `wHANTS` is more stable for multiple growing seasons, but it is
  less flexible than `wWHIT.` `wHANTS` is suitable for regions with the static
  growing season pattern across multiple years, `wWHIT` is more suitable for
  regions with the dynamic growing season pattern. Dynamic growing season
  pattern is the most challenging task, which also means that a large
  uncertainty might exist.

  When using `wHANTS` as the rough fitting function, `r_min` is suggested to be
  set as zero.

- Use only one iteration in the fine fitting procedure.


# **References**

> [1]. Dongdong Kong, Tim R. McVicar, Mingzhong Xiao\*, Yongqiang Zhang, Jorge L. Peña-Arancibia, Gianluca Filippa, Yuxuan Xie, and Xihui Gu\*. (2022), 
> phenofit: A R package for extracting vegetation phenology from time series remote sensing. __*Methods in Ecology and Evolution*__. (Accepted April 08, 2022). <https://doi.org/10.1111/2041-210X.13870>
> 
> [2] Kong, D., Zhang, Y.\*, Wang, D., Chen, J., & Gu, X\*. (2020).
> Photoperiod Explains the Asynchronization Between Vegetation Carbon
> Phenology and Vegetation Greenness Phenology. __*Journal of Geophysical
> Research: Biogeosciences*__, 125(8), e2020JG005636.
> <https://doi.org/10.1029/2020JG005636>
>
> [3] Kong, D., Zhang, Y.\*, Gu, X., & Wang, D. (2019). A robust method
> for reconstructing global MODIS EVI time series on the Google Earth
> Engine. __*ISPRS Journal of Photogrammetry and Remote Sensing*__, 155,
> 13–24.
>
> [4] Kong, D., (2020). R package: A state-of-the-art Vegetation
> Phenology extraction package, `phenofit` version 0.3.5,
> <https://doi.org/10.5281/zenodo.6320537>
>
> [5] Zhang, Q.\*, Kong, D.\*, Shi, P., Singh, V.P., Sun, P., 2018.
> Vegetation phenology on the Qinghai-Tibetan Plateau and its response
> to climate change (1982–2013). __*Agricultural and Forest Meteorology*__. 248, 408–417.
> <https://doi.org/10.1016/j.agrformet.2017.10.026>


# Acknowledgements

Keep in mind that this repository is released under a GPL2 license,
which permits commercial use but requires that the source code (of
derivatives) is always open even if hosted as a web service.
