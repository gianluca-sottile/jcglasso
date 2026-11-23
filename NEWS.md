# jcglasso 1.0.0

* First CRAN release. The package implements the joint conditional graphical
  lasso estimator of Sottile et al. (2025, JRSS-C) and is a stand-alone
  evolution of the development scripts built on top of 'cglasso':
  - numerical kernels rewritten in C (`.Call`, BLAS/LAPACK of R): E-step,
    accelerated proximal gradient for the regression coefficients, joint
    graphical lasso (ADMM on the connected components of the screening rule);
  - the M-step uses the expected second moments of the E-step also for the
    regression coefficients (the scripts used the products of the imputed
    predictors), and the joint mean of responses and predictors in the E-step;
  - the tolerance of the M-step is tied to the EM tolerance
    (`min(thr.bcd, thr.em / 10)`), which removes the stalls of the EM algorithm
    observed when both tolerances were equal;
  - log-scale Mills ratios in the E-step (no 0/0 in the tails), missing values
    coded as `NA` (no more `-9999`);
  - compact output (estimates, Q-function values and degrees of freedom along the
    path); imputed data, fitted values and residuals are computed on request
    (`impute()`, `fitted()`, `residuals()`);
  - Q-function with conditional and joint components, AIC/BIC/eBIC, model
    selection, maximum likelihood refit (`jcggm()`), KKT check (`kkt()`), graphs
    (`to_graph()`, with 'igraph' in Suggests), simulation (`rjcggm()`, no 'MASS');
  - accessors are S3 generics whose default methods call 'cglasso' (if
    installed), so that both packages can be used together.
