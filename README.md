
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SANvi v0.0.1 <img src="man/figures/sanvi_draft.png" align="right" width="120" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/Fradenti/SANvi/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Fradenti/SANvi/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of SANvi is to …

## Installation

You can install the development version of SANvi from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Fradenti/SANvi")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(SANvi)
#> Loading required package: scales
#> Loading required package: RColorBrewer
# Generate example data
set.seed(1232)
y <- c(rnorm(100),rnorm(100,5))
g <- rep(1:2,rep(100,2))

# Fitting fiSAN via variational inference
est <- SANvi:::variational_fiSAN(y,g,verbose = FALSE)

# Estimate posterior atoms and weights
cl <- estimate_atoms_weights_vi(est)
plot(cl)
```

<img src="man/figures/README-example-1.png" width="100%" />
