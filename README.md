# PCEV
[![Build Status](https://travis-ci.org/GreenwoodLab/pcev.svg?branch=master)](https://travis-ci.org/GreenwoodLab/pcev)


R package which implements Principal components of explained variance (PCEV).

PCEV is a statistical tool for the analysis of a mutivariate response vector. It is a dimension-reduction technique, similar to Principal Components Analysis (PCA), which seeks to maximize the proportion of variance (in the response vector) being explained by a set of covariates.

## Installation

You can install from GitHub using the [devtools](http://cran.r-project.org/web/packages/devtools/index.html) package:

``` r
library(devtools)
devtools::install_github('GreenwoodLab/pcev', build_vignettes = TRUE)
```

The main function is ```computePCEV```, and indeed most users will only need this one function. See the documentation for more information about its parameters and for some examples.
