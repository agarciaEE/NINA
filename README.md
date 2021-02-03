
<!-- README.md is generated from README.Rmd. Please edit that file -->

# NINA

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/agarciaEE/NINA.svg?branch=main)](https://travis-ci.com/agarciaEE/NINA)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/agarciaEE/NINA?branch=main&svg=true)](https://ci.appveyor.com/project/agarciaEE/NINA)
[![Codecov test
coverage](https://codecov.io/gh/agarciaEE/NINA/branch/main/graph/badge.svg)](https://codecov.io/gh/agarciaEE/NINA?branch=main)
<!-- badges: end -->

The goal of NINA is the the analysis and implementation of biotic
interactions into environmental niche models using kernel density
estimations.

## Installation

You can install the released version of NINA from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("NINA")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(NINA)
#> Registered S3 methods overwritten by 'adehabitatMA':
#>   method                       from
#>   print.SpatialPixelsDataFrame sp  
#>   print.SpatialPixels          sp
#> Registered S3 method overwritten by 'spatstat':
#>   method     from
#>   print.boxx cli
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/master/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
