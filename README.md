# PCEV
[![Build Status](https://travis-ci.org/GreenwoodLab/pcev.svg?branch=master)](https://travis-ci.org/GreenwoodLab/pcev) [![CRAN](http://www.r-pkg.org/badges/version/pcev?color=blue)](http://cran.rstudio.com/package=pcev) [![Downloads](http://cranlogs.r-pkg.org/badges/pcev?color=green)](http://www.r-pkg.org/pkg/pcev)


R package which implements Principal components of explained variance (PCEV).

PCEV is a statistical tool for the analysis of a mutivariate response vector. It is a dimension-reduction technique, similar to Principal Components Analysis (PCA), which seeks to maximize the proportion of variance (in the response vector) being explained by a set of covariates.

For more information you can look at the [vignette](https://cran.rstudio.com/web/packages/pcev/vignettes/pcev.pdf). Alternatively, if you have already installed the package along with the vignette, you can access the vignette from within ```R``` by using the following command:

``` r
vignette("pcev")
```

## Installation

This package is available on [CRAN](https://cran.r-project.org/web/packages/pcev/). Alternatively, you can install from GitHub using the [devtools](http://cran.r-project.org/web/packages/devtools/index.html) package:

``` r
library(devtools)
devtools::install_github('GreenwoodLab/pcev', build_vignettes = TRUE)
```

The main function is ```computePCEV```, and indeed most users will only need this one function. See the documentation for more information about its parameters and for some examples.
