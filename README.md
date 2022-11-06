# <img src="https://github.com/stamats/MKinfer/raw/master/hex-MKinfer.png" alt="MKinfer" width="120"/> &emsp; MKinfer
The repository includes the development version of R package MKinfer

[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/MKinfer)](http://cran.r-project.org/package=MKinfer)
[![cran checks](https://cranchecks.info/badges/summary/MKinfer)](https://cran.r-project.org/web/checks/check_results_MKinfer.html)

## Description
Computation of various confidence intervals (Altman et al. (2000), ISBN:978-0-727-91375-3; 
Hedderich and Sachs (2018), ISBN:978-3-662-56657-2) including bootstrapped versions 
(Davison and Hinkley (1997), ISBN:978-0-511-80284-3) as well as 
Hsu (Hedderich and Sachs (2018), ISBN:978-3-662-56657-2), 
permutation (Janssen (1997), <doi:10.1016/S0167-7152(97)00043-6>), 
bootstrap (Davison and Hinkley (1997), ISBN:978-0-511-80284-3) and 
multiple imputation (Barnard and Rubin (1999), <doi:10.1093/biomet/86.4.948>) t-test.
Graphical visualization by volcano plots.

## Installation

```{r, eval = FALSE}
## Installation of CRAN version
install.packages("MKinfer")

## Development version from GitHub
# install.packages("remotes")
remotes::install_github("stamats/MKinfer")
```

## Getting started

```{r}
library(MKinfer)
```
