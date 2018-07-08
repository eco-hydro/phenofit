## phenofit

[![Travis Build
Status](https://travis-ci.org/kongdd/phenofit.svg?branch=master)](https://travis-ci.org/kongdd/phenofit)
[![codecov](https://codecov.io/gh/kongdd/phenofit/branch/master/graph/badge.svg)](https://codecov.io/gh/kongdd/phenofit)

A state-of-the-art **remote sensing vegetation phenology** extraction package: `phenofit`


 - `phenofit` combine merits of TIMESAT and phenopix
 - A simple and stable growing season dividing methods was proposed
 - Provide a practical snow elimination method, based on Whittaker
 - 7 curve fitting methods and 4 phenology extraction methods
 - We add parameters boundary for every curve fitting methods according to their ecological meaning.
 - `optimx` is used to select best optimization method for different curve fitting methods.


## Installation
```r
devtools::install_github('kongdd/phenofit')
```
