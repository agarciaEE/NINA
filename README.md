
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

You can install NINA from github repository with:

``` r
devtools::install_github("agarciaEE/NINA")
```

## Example

This is a basic example which shows you how to use NINA’s functions:

Estimate Environmental Niche Models:

Correct EN models by biotic interactions:

``` r
g2_BC <- BC_model(g2_EN, g1_EN, A.matrix = int_matrix, type = "region")
#> Estimating biotic constrains of species in region A...
#>  Adding biotic constrains to sp1...  ...Success!
#>  Adding biotic constrains to sp4...  ...Success!
#> Estimating biotic constrains of species in region B...
#>  Adding biotic constrains to sp1...  ...Success!
#>  Adding biotic constrains to sp2...  ...Success!
#>  Adding biotic constrains to sp3...  ...Success!
#>  Adding biotic constrains to sp5...  ...Success!
#> Estimating biotic constrains of species in region C...
#>  Adding biotic constrains to sp1...  ...Success!
#>  Adding biotic constrains to sp2...  ...Success!
#>  Adding biotic constrains to sp3...  ...Success!
#>  Adding biotic constrains to sp5...  ...Success!
#> Estimating biotic constrains of species in region D...
#>  Adding biotic constrains to sp1...  ...Success!
#>  Adding biotic constrains to sp2...  ...Success!
#>  Adding biotic constrains to sp3...  ...Success!
#>  Adding biotic constrains to sp5...  ...Success!
#> Estimating biotic constrains of species in region E...
#>  Adding biotic constrains to sp1...  ...Success!
#>  Adding biotic constrains to sp4...  ...Success!
#>  Adding biotic constrains to sp2...  ...Success!
#>  Adding biotic constrains to sp3...  ...Success!
#>  - Carrying out models evaluations...
#> Observations ... OK
#> Environmental predictors ... OK
#> Model predictions ... OK
#> Sampling pseudo_absences... 
#> Performing models evaluation...
#>  ...Evaluating sp1 niche model...
#>  ...Evaluating sp4 niche model...
#>  ...Evaluating sp2 niche model...
#>  ...Evaluating sp3 niche model...
#> Warning in cor(x, y): the standard deviation is zero
#>  ...Evaluating sp5 niche model...
#> Models evaluation performed.
#> Models successfully corrected!
```

Transform environmental niche space into ecological niche space
(prioritizing the effect of biotic interactions).

``` r
g2_EC <- EC_model(g2_BC, type = "region")
#> Transforming niche space of species in region A...
#>  Estimating ecological niche of sp1...   ...Success!
#>  Estimating ecological niche of sp4...   ...Success!
#> Transforming niche space of species in region B...
#>  Estimating ecological niche of sp1...   ...Success!
#>  Estimating ecological niche of sp2...   ...Success!
#>  Estimating ecological niche of sp3...   ...Success!
#>  Estimating ecological niche of sp5...   ...Success!
#> Transforming niche space of species in region C...
#>  Estimating ecological niche of sp1...   ...Success!
#>  Estimating ecological niche of sp2...   ...Success!
#>  Estimating ecological niche of sp5...   ...Success!
#> Transforming niche space of species in region D...
#>  Estimating ecological niche of sp1...   ...Success!
#>  Estimating ecological niche of sp2...   ...Success!
#>  Estimating ecological niche of sp3...   ...Success!
#>  Estimating ecological niche of sp5...   ...Success!
#> Transforming niche space of species in region E...
#>  Estimating ecological niche of sp1...   ...Success!
#>  Estimating ecological niche of sp4...   ...Success!
#>  Estimating ecological niche of sp2...   ...Success!
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
#> TRUE
#> Geographical extents:
#> 5
#>   spA spB spC spD spE
#> A   0   1   1   1   0
#> B   1   1   1   1   1
#> C   0   1   1   1   0
#> D   1   1   1   1   1
#> E   0   1   1   1   1
#> Ensemble of regional models:
#> FALSE
#> Failures or warnings:
#> TRUE
#>   region species
#> 1      A     spE
#> 2      C     spA
#> 3      E     spA
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
#> TRUE
#> Geographical extents:
#> 5
#> Ensemble of regional models:
#> FALSE
#> Number of species:
#> 5
#> Failures or warnings:
#> TRUE
#> Models evaluation:
#> TRUE
```

You can also plot a summary of the output model using `plot` function:

``` r
plot(g2_EN)
#> Ploting only the first four species maps of a total of 5
```

Models can be evaluated using NINA’s function `models_evaluation` and
visualize it by plotting the output.
