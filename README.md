
<!-- README.md is generated from README.Rmd. Please edit that file -->

# baydem

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Travis build
status](https://travis-ci.org/eehh-stanford/baydem.svg?branch=master)](https://travis-ci.org/eehh-stanford/baydem)
[![Codecov test
coverage](https://codecov.io/gh/eehh-stanford/baydem/branch/master/graph/badge.svg)](https://codecov.io/gh/eehh-stanford/baydem?branch=master)
<!-- badges: end -->

Bayesian tools for reconstructing past and present demography. This is
the package used in Price et al. (2020):

> Price, M.H., J.M. Capriles, J. Hoggarth, R.K. Bocinsky, C.E. Ebert,
> and J.H. Jones (2020). *End-to-end Bayesian analysis of radiocarbon
> dates reveals new insights into lowland Maya demography*. **In
> review.**

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("eehh-stanford/baydem")
```

## Analysis

The [:file\_folder: analysis](/analysis) directory contains the code,
data, and output for Price et al. 2020.

To re-create the analysis, please run the following in *R*:

``` r
# install.packages("devtools")
devtools::install_github("eehh-stanford/baydem")

library("baydem")

list.files("analysis", 
           full.names = TRUE, 
           pattern = "FINAL") %>%
           purrr::walk(source)
```

## Contributor Code of Conduct

Please note that the ‘baydem’ project is released with a [Contributor
Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project,
you agree to abide by its terms.
