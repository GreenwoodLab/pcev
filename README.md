# PCEV
[![Build Status](https://travis-ci.org/GreenwoodLab/pcev.svg?branch=master)](https://travis-ci.org/GreenwoodLab/pcev) [![CRAN](http://www.r-pkg.org/badges/version/pcev?color=blue)](http://cran.rstudio.com/package=pcev) [![Downloads](http://cranlogs.r-pkg.org/badges/pcev?color=green)](http://www.r-pkg.org/pkg/pcev)


R package which implements Principal components of explained variance (PCEV).

PCEV is a statistical tool for the analysis of a mutivariate response vector. It is a dimension-reduction technique, similar to Principal Components Analysis (PCA), that seeks to maximize the proportion of variance (in the response vector) being explained by a set of covariates. It implements three versions:

 - the classic version, when p < n;
 - the singular version, when p > n;
 - the block version, our extension of the algorithm for the case of a high number of data points (p>>n).
 
For all three versions, we provide hypothesis testing based on Roy's largest root.

For more information you can look at the [vignette](https://cran.rstudio.com/web/packages/pcev/vignettes/pcev.pdf). Alternatively, if you have already installed the package along with the vignette, you can access the vignette from within ```R``` by using the following command:

``` r
vignette("pcev")
```

## Installation

This package is available on [CRAN](https://cran.r-project.org/package=pcev). Alternatively, you can install from GitHub using the [devtools](http://cran.r-project.org/package=devtools) package:

``` r
library(devtools)
devtools::install_github('GreenwoodLab/pcev', build_vignettes = TRUE)
```

The main function is ```computePCEV```, and indeed most users will only need this one function. See the documentation for more information about its parameters and for some examples.

### Singular case

Note that hypothesis testing in the singular case can be performed using an approximation of Roy's Largest Root's distribution. This approximation uses the distribution of the largest eigenvalue of a Wishart matrix, and it is implemented using the package [rootWishart](https://cran.r-project.org/package=rootWishart). This package is only suggested, and therefore you will need to install it separately if you want to take advantage of this approximation. Otherwise, the p-value is computed using a combination of permutations, method of moments, and the Tracy-Widom distribution.