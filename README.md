
<!-- README.md is generated from README.Rmd. Please edit that file -->


# jcglasso

<!-- badges: start -->
[![R-CMD-check](https://github.com/gianluca-sottile/jcglasso/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gianluca-sottile/jcglasso/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

'jcglasso' implements the joint conditional graphical lasso of Sottile et al.
(2025): for `K` sub-populations with censored and/or missing values, it
estimates sparse networks among the predictors (Omega_k), sparse effects
of the predictors on the responses (B_k) and sparse networks among the
responses given the predictors (Theta_k), encouraging similar structures
across sub-populations.

## Installation

```r
# development version
remotes::install_github("gianluca-sottile/jcglasso")
```

## Example

``` r
library(jcglasso)
set.seed(1234)
sim <- rjcggm(n = 100, p = 10, q = 4, K = 2, ncompPrecision = 2, ncompB = 2)
fitx <- jcglasso(data = sim$ZX, nrho = 10, rho.min.ratio = 0.1)
nu <- select_jcglasso(fitx, GoF = BIC)$rho
fit <- jcglasso(data = sim$Z, nrho = 5, nlambda = 5, rho.min.ratio = 0.1,
                lambda.min.ratio = 0.1, nu = nu)
best <- select_jcglasso(fit, GoF = BIC)
mle <- jcggm(best)
round(coef(mle, type = "B", class.id = 1)[, 1:5], 2)
#>                Y1    Y2    Y3    Y4    Y5
#> (Intercept) 33.90 34.15 34.17 39.19 33.63
#> X1           0.00  0.41  0.00  0.63  0.00
#> X2           0.54  0.55  0.00 -0.07  0.65
#> X3           0.00  0.00  0.90 -0.31 -0.15
#> X4           0.56  0.20  0.89  0.46  0.66
```

## Reference

Sottile, G., Augugliaro, L., Vinciotti, V., Arancio, W. and Coronnello, C. (2025).
Sparse inference of the human haematopoietic system from heterogeneous and
partially observed genomic data. *Journal of the Royal Statistical Society Series
C: Applied Statistics*, 74(1), 204-228. <https://doi.org/10.1093/jrsssc/qlae052>
