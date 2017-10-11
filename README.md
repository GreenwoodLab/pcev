# PCEV
[![Build Status](https://travis-ci.org/GreenwoodLab/pcev.svg?branch=master)](https://travis-ci.org/GreenwoodLab/pcev) [![CRAN](http://www.r-pkg.org/badges/version/pcev?color=blue)](https://cran.r-project.org/package=pcev) [![Downloads](http://cranlogs.r-pkg.org/badges/pcev?color=green)](https://cran.r-project.org/package=pcev)


R package which implements Principal components of explained variance (PCEV).

PCEV is a statistical tool for the analysis of a mutivariate response vector. It is a dimension-reduction technique, similar to Principal Components Analysis (PCA), that seeks to maximize the proportion of variance (in the response vector) being explained by a set of covariates. It implements three versions:

 - the classic version, when p < n;
 - the singular version, when p > n;
 - the block version, our extension of the algorithm for the case of a high number of data points (p>>n).
 
For the first two versions, we provide hypothesis testing based on Roy's largest root.

For more information you can look at the [vignette](https://cran.r-project.org/package=pcev/vignettes/pcev.pdf). Alternatively, if you have already installed the package along with the vignette, you can access the vignette from within ```R``` by using the following command:

``` r
vignette("pcev")
```

## Installation

This package is available on [CRAN](https://cran.r-project.org/package=pcev). Alternatively, you can install from GitHub using the [devtools](https://cran.r-project.org/package=devtools) package:

``` r
library(devtools)
devtools::install_github('GreenwoodLab/pcev', build_vignettes = TRUE)
```

The main function is ```computePCEV```, and indeed most users will only need this one function. See the documentation for more information about its parameters and for some examples.

## References

  - Turgeon, M., Oualkacha, K., Ciampi, A., Miftah, H., Dehghan, G., Zanke, B.W., Benedet, A.L., Rosa-Neto, P., Greenwood, C.M.T., Labbe, A., for the Alzheimer’s Disease Neuroimaging Initiative. “Principal component of explained variance: an efficient and optimal data dimension reduction framework for association studies”. To appear in [Statistical Methods in Medical Research](http://dx.doi.org/10.1177/0962280216660128).