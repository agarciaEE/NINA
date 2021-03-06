
<!-- README.md is generated from README.Rmd. Please edit that file -->

# NINA package

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/agarciaEE/NINA.svg?branch=main)](https://travis-ci.com/agarciaEE/NINA)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/agarciaEE/NINA?branch=main&svg=true)](https://ci.appveyor.com/project/agarciaEE/NINA)
[![Codecov test
coverage](https://codecov.io/gh/agarciaEE/NINA/branch/main/graph/badge.svg)](https://codecov.io/gh/agarciaEE/NINA?branch=main)
[![CRAN
status](https://www.r-pkg.org/badges/version/NINA)](https://CRAN.R-project.org/package=NINA)
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

This is a basic example which shows you how to use NINA’s functions:

``` r
library(NINA)
#> Registered S3 methods overwritten by 'adehabitatMA':
#>   method                       from
#>   print.SpatialPixelsDataFrame sp  
#>   print.SpatialPixels          sp
## basic example code
```

Estimate Environmental Niche Models:

``` r
g1_EN = EN_model(env_data, occ_data1)
#> Carrying out unique EN model...
#> Occurrence dataset ... OK
#> Environmental dataset ... OK
#>  - Conforming environmental space...
#>  - Carrying out species EN models...
#>      - Modelling spA Environmental Niche...
#>      - Modelling spB Environmental Niche...
#>      - Modelling spC Environmental Niche...
#>      - Modelling spD Environmental Niche...
#>      - Modelling spE Environmental Niche...
#> Species EN models succesfully completed!
g2_EN = EN_model(env_data, occ_data2)
#> Carrying out unique EN model...
#> Occurrence dataset ... OK
#> Environmental dataset ... OK
#>  - Conforming environmental space...
#>  - Carrying out species EN models...
#>      - Modelling sp1 Environmental Niche...
#>      - Modelling sp2 Environmental Niche...
#>      - Modelling sp3 Environmental Niche...
#>      - Modelling sp4 Environmental Niche...
#>      - Modelling sp5 Environmental Niche...
#> Species EN models succesfully completed!
```

Correct EN models by biotic interactions:

``` r
g2_BC <- BC_model(g2_EN, g1_EN, A.matrix = int_matrix,  type = "global")
#>  Adding biotic constrains to sp1...  ...Success!
#>  Adding biotic constrains to sp2...  ...Success!
#>  Adding biotic constrains to sp3...  ...Success!
#>  Adding biotic constrains to sp4...  ...Success!
#>  Adding biotic constrains to sp5...  ...Success!
#> Models successfully corrected!
```

Transform environmental niche space into ecological niche space
(prioritizing the effect of biotic interactions).

``` r
g2_EC <- EC_model(g2_BC, type = "global")
#> Estimating ecological niche of sp1...    ...Success!
#> Estimating ecological niche of sp2...    ...Success!
#> Estimating ecological niche of sp3...    ...Success!
#> Estimating ecological niche of sp4...    ...Success!
#> Estimating ecological niche of sp5...    ...Success!
#> Models successfully transformed!
```

You can summarize the output by using `summary` function or `print`:

``` r
summary(g1_EN)
#> Niche model type:
#> Environmental-only
#> Predictors:
#> bio1 bio2 bio3 bio4 bio5 bio6
#> 2 axis-components used
#> Class: pca dudi
#> Call: dudi.pca(df = stats::na.exclude(env[, env.var]), center = T, 
#>     scale = T, scannf = F, nf = 2)
#> 
#> Total inertia: 6
#> 
#> Eigenvalues:
#>      Ax1      Ax2      Ax3      Ax4      Ax5 
#> 3.418906 1.779810 0.773672 0.018853 0.007016 
#> 
#> Projected inertia (%):
#>     Ax1     Ax2     Ax3     Ax4     Ax5 
#> 56.9818 29.6635 12.8945  0.3142  0.1169 
#> 
#> Cumulative projected inertia (%):
#>     Ax1   Ax1:2   Ax1:3   Ax1:4   Ax1:5 
#>   56.98   86.65   99.54   99.85   99.97 
#> 
#> (Only 5 dimensions (out of 6) are shown)
#> 
#> Spatially constrained:
#> FALSE
#> Geographical extents:
#> 1
#> 
#> spA spB spC spD spE 
#>   1   1   1   1   1 
#> Ensemble of regional models:
#> FALSE
#> Failures or warnings:
#> FALSE
#> Models evaluation:
#> FALSE
print(g2_BC)
#> Object class: NINA
#> 
#> Niche model type:
#> Environmental-constrained
#> Predictors:
#> bio1 bio2 bio3 bio4 bio5 bio6
#> Spatially constrained:
#> FALSE
#> Geographical extents:
#> 1
#> Ensemble of regional models:
#> FALSE
#> Number of species:
#> 5
#> Failures or warnings:
#> FALSE
#> Models evaluation:
#> FALSE
```

You can also plot a summary of the output model using `plot` function:

    #> Ploting only the first four species maps of a total of 5

<img src="man/figures/README-plot-1.png" width="100%" />

Models can be evaluated using NINA’s function `models_evaluation` and
visualize it by plotting the output.

``` r
eval <- models_evaluation(g2_BC, plot = F)
#> Observations ... OK
#> Environmental predictors ... OK
#> Model predictions ... OK
#> Sampling pseudo_absences...
#> Performing models evaluation...
#>  ...Evaluating sp1 niche model...    ...Evaluating sp2 niche model...    ...Evaluating sp3 niche model...    ...Evaluating sp4 niche model...    ...Evaluating sp5 niche model...Models evaluation performed.
plot(eval)
```

<img src="man/figures/README-eval-1.png" width="100%" />
