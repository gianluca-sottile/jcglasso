#' jcglasso: Joint Conditional Graphical Lasso with Censored and Missing Values
#'
#' The package implements the joint conditional graphical lasso estimator of
#' Sottile et al. (2025) for `K` sub-populations of a conditional Gaussian
#' graphical model with censored and/or missing values: for each sub-population,
#' a sparse precision matrix of the predictors, sparse regression coefficients and
#' a sparse precision matrix of the responses given the predictors are estimated,
#' with group lasso penalties that encourage similar structures across the
#' sub-populations.
#'
#' The main functions are [datajcggm()] (data), [jcglasso()] (estimation along a
#' path of tuning parameters), [BIC.jcglasso()] and [select_jcglasso()] (model
#' selection), [jcggm()] (maximum likelihood refit), [to_graph()] (graphs) and
#' [rjcggm()] (simulation). The model with `K = 1` is the conditional graphical
#' lasso with censored and missing values of package 'cglasso' (Augugliaro et al.,
#' 2023). Functions with the same names in the two packages (e.g., [getMatrix()],
#' [QFun()], [to_graph()]) are S3 generics in 'jcglasso' whose default methods call
#' the functions of 'cglasso': to use both packages, attach 'jcglasso' after
#' 'cglasso'.
#'
#' @references Sottile, G., Augugliaro, L., Vinciotti, V., Arancio, W. and
#'   Coronnello, C. (2025). Sparse inference of the human haematopoietic system
#'   from heterogeneous and partially observed genomic data. *Journal of the Royal
#'   Statistical Society Series C: Applied Statistics*, 74(1), 204-228.
#'   \doi{10.1093/jrsssc/qlae052}
#'
#'   Augugliaro, L., Sottile, G., Wit, E. C. and Vinciotti, V. (2023). cglasso:
#'   An R package for conditional graphical lasso inference with censored and
#'   missing values. *Journal of Statistical Software*, 105(1), 1-58.
#'   \doi{10.18637/jss.v105.i01}
#' @useDynLib jcglasso, .registration = TRUE, .fixes = ""
#' @importFrom stats AIC BIC coef fitted nobs predict residuals
#' @keywords internal
"_PACKAGE"
