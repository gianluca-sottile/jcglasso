#' Joint Conditional Graphical Lasso with Censored and Missing Values
#'
#' Fits the joint conditional graphical lasso estimator of Sottile et al. (2025)
#' to a dataset with `K` classes and censored and/or missing values, over a path
#' of tuning parameters. For each class `k`, the predictors follow a Gaussian
#' graphical model with precision matrix \eqn{\Omega_k} and the responses follow a
#' conditional Gaussian graphical model with regression coefficients \eqn{B_k} and
#' precision matrix \eqn{\Theta_k}; similarity across classes is encouraged by
#' sparse group lasso penalties.
#'
#' @param formula an optional formula `cbind(Y1, Y2, ...) ~ X1 + X2 + ...` to
#'   select responses and predictors of `data`; `.` means all of them (default
#'   `. ~ .`). The right-hand side can only contain names of predictors (no
#'   transformations, since predictors can be censored or missing); use `~ 1` to
#'   fit the model without predictors.
#' @param data an object of class `"datajcggm"`.
#' @param diagonal logical: should the diagonal entries of the precision matrices
#'   be penalized? Default `FALSE`.
#' @param weights.B,weights.Tht,weights.Omega optional arrays (or lists of `K`
#'   matrices) of non-negative weights of the penalties of \eqn{B} (`q x p x K`),
#'   \eqn{\Theta} (`p x p x K`) and \eqn{\Omega} (`q x q x K`). A weight equal to
#'   `Inf` defines a structural zero, a weight equal to `0` an unpenalized
#'   parameter. By default, all weights are equal to one (and the diagonals of the
#'   precision matrices are unpenalized).
#' @param penalty the penalty used across classes for the precision matrices:
#'   currently only `"group"` (the group lasso penalty of Danaher et al., 2014).
#' @param nrho,rho.min.ratio,rho the number of values, the ratio between the
#'   smallest and the largest value, or a user supplied decreasing sequence of
#'   values of the tuning parameter \eqn{\rho} of the precision matrices
#'   \eqn{\Theta_k}.
#' @param nlambda,lambda.min.ratio,lambda the same for the tuning parameter
#'   \eqn{\lambda} of the regression coefficients \eqn{B_k} (only with predictors).
#' @param nu the tuning parameter of the precision matrices \eqn{\Omega_k} of the
#'   predictors: a single value, or `NULL` (default) to use \eqn{\nu = \rho}
#'   along the path. Sottile et al. (2025) suggest to select \eqn{\nu} by fitting
#'   the model of the predictors alone (see the examples).
#' @param alpha1,alpha2,alpha3 the mixing parameters in \eqn{[0, 1]} between the
#'   lasso (`alpha = 1`) and the group lasso (`alpha = 0`) penalties of
#'   \eqn{\Theta}, \eqn{B} and \eqn{\Omega}, respectively.
#' @param maxit.em,thr.em maximum number of iterations and tolerance of the EM
#'   algorithm.
#' @param maxit.bcd,thr.bcd maximum number of iterations and tolerance of the
#'   algorithms of the M-step.
#' @param trace integer: `0` no output, `1` one line for each fitted model, `2`
#'   also the EM iterations.
#'
#' @details
#' Let \eqn{f_k = n_k / n}. For given tuning parameters, the estimator maximizes
#' the observed log-likelihood penalized by
#' \eqn{\lambda P_{\alpha_2}(\{B\}) + \rho \tilde P_{\alpha_1}(\{\Theta\}) + \nu \tilde P_{\alpha_3}(\{\Omega\})},
#' where, e.g.,
#' \deqn{P_\alpha(\{B\}) = \alpha \sum_k \sum_{h,j} w_{hjk}|\beta_{hjk}| + (1 - \alpha)\sum_{h,j} w_{hj1}\Big(\sum_k \beta_{hjk}^2\Big)^{1/2},}
#' and the penalties of the precision matrices exclude the diagonal entries
#' (unless `diagonal = TRUE`). The EM algorithm of Sottile et al. (2025, Section
#' 5) is used: the E-step computes the conditional expected values of the
#' sufficient statistics (with the approximation of Guo et al., 2015, for the
#' second moments of censored values); the M-step updates the means, then
#' \eqn{\Omega} by the joint graphical lasso (ADMM of Danaher et al., 2014, applied
#' to the connected components of the screening rule), \eqn{B} given \eqn{\Theta}
#' by an accelerated proximal gradient algorithm, and \eqn{\Theta} given \eqn{B} by
#' the joint graphical lasso. The regression part of the objective function is
#' scaled by one half, so that the smallest value of \eqn{\lambda} giving
#' \eqn{B = 0} is \eqn{\max_{h,j}\{\sum_k (f_k |S_{xy,k}\Theta_k|_{hj})^2\}^{1/2}}.
#'
#' The path is computed for decreasing values of \eqn{\lambda} (outer loop) and
#' \eqn{\rho} (inner loop), with warm starts. Only the estimated parameters are
#' stored; the imputed data and the fitted values can be computed with
#' [impute()], [fitted.jcglasso()] and [residuals.jcglasso()].
#'
#' @return An object of class `"jcglasso"`, a list with components
#' \item{B}{array `(q + 1) x p x K x nlambda x nrho` of regression coefficients
#'   (first row: intercepts).}
#' \item{Tht}{array `p x p x K x nlambda x nrho` of precision matrices of the responses.}
#' \item{Omega}{array `q x q x K x nlambda x nrho` of precision matrices of the
#'   predictors (`NULL` without predictors).}
#' \item{psi}{array `(p + q) x K x nlambda x nrho` of expected values.}
#' \item{Q}{list of arrays `K x nlambda x nrho` with the values of the
#'   Q-function of the conditional model (`yx`) and of the predictors (`x`).}
#' \item{dfB, dfTht, dfOmg}{numbers of non-zero (off-diagonal) estimates.}
#' \item{nit}{array `2 x nlambda x nrho` with the numbers of EM steps and of
#'   iterations of the M-step algorithms.}
#' \item{conv}{matrix `nlambda x nrho` of convergence codes (`0` converged).}
#' \item{lambda, rho, nu}{the tuning parameters.}
#' and `call`, `Z` (the data used), the weights and the other arguments.
#' @references Sottile, G., Augugliaro, L., Vinciotti, V., Arancio, W. and
#'   Coronnello, C. (2025). Sparse inference of the human haematopoietic system
#'   from heterogeneous and partially observed genomic data. *Journal of the Royal
#'   Statistical Society Series C: Applied Statistics*, 74(1), 204-228.
#'   \doi{10.1093/jrsssc/qlae052}
#'
#'   Danaher, P., Wang, P. and Witten, D. M. (2014). The joint graphical lasso
#'   for inverse covariance estimation across multiple classes. *Journal of the
#'   Royal Statistical Society Series B*, 76(2), 373-397. \doi{10.1111/rssb.12033}
#'
#'   Augugliaro, L., Sottile, G., Wit, E. C. and Vinciotti, V. (2023). cglasso:
#'   An R package for conditional graphical lasso inference with censored and
#'   missing values. *Journal of Statistical Software*, 105(1), 1-58.
#'   \doi{10.18637/jss.v105.i01}
#' @seealso [datajcggm()], [rjcggm()], [jcggm()], [BIC.jcglasso()],
#'   [select_jcglasso()], [coef.jcglasso()], [to_graph()]
#' @examples
#' set.seed(1234)
#' sim <- rjcggm(n = 100, p = 10, q = 4, K = 2, ncompPrecision = 2, ncompB = 2)
#'
#' # 1. model of the predictors (choice of nu)
#' fitx <- jcglasso(data = sim$ZX, nrho = 10, rho.min.ratio = 0.1)
#' bestx <- select_jcglasso(fitx, GoF = BIC)
#'
#' # 2. conditional model, with nu fixed
#' fit <- jcglasso(data = sim$Z, nrho = 5, rho.min.ratio = 0.1, nlambda = 5,
#'                 lambda.min.ratio = 0.1, nu = bestx$rho)
#' fit
#' best <- select_jcglasso(fit, GoF = BIC)
#' coef(best, type = "B")
#' @export
jcglasso <- function(formula, data, diagonal = FALSE, weights.B = NULL, weights.Tht = NULL, weights.Omega = NULL,
                     penalty = "group", nrho = 10L, rho.min.ratio, rho, nlambda = 10L, lambda.min.ratio, lambda,
                     nu = NULL, alpha1 = 0.5, alpha2 = 0.5, alpha3 = 0.5,
                     maxit.em = 1e3, thr.em = 1e-4, maxit.bcd = 1e4, thr.bcd = 1e-4, trace = 0L) {
  this.call <- match.call()
  if (missing(data) || !is.datajcggm(data)) stop(sQuote("data"), " must be an object of class ",
                                                 sQuote("datajcggm"), call. = FALSE)
  penalty <- match.arg(penalty, "group")
  if (missing(formula)) formula <- . ~ .
  Z <- subset_by_formula(data, formula)
  K <- length(Z$classes)
  n <- nobs(Z)
  p <- ncol(Z$Y[[1L]])
  q <- if (is.null(Z$X)) 0L else ncol(Z$X[[1L]])
  if (p < 2L && q == 0L) stop("at least two variables are needed", call. = FALSE)
  check_scalar(alpha1, "alpha1", 0, 1); check_scalar(alpha2, "alpha2", 0, 1); check_scalar(alpha3, "alpha3", 0, 1)
  check_logical(diagonal, "diagonal")
  check_int(maxit.em, "maxit.em"); check_int(maxit.bcd, "maxit.bcd")
  check_scalar(thr.em, "thr.em", 0, Inf, open = TRUE); check_scalar(thr.bcd, "thr.bcd", 0, Inf, open = TRUE)
  if (!trace %in% 0:2) stop(sQuote("trace"), " must be 0, 1 or 2", call. = FALSE)
  ynames <- colnames(Z$Y[[1L]]); xnames <- if (q) colnames(Z$X[[1L]])
  wTht <- make_weights(weights.Tht, c(p, p, K), list(ynames, ynames, Z$classes), "weights.Tht", diagonal, square = TRUE)
  wB <- if (q) make_weights(weights.B, c(q, p, K), list(xnames, ynames, Z$classes), "weights.B")
  wOmg <- if (q) make_weights(weights.Omega, c(q, q, K), list(xnames, xnames, Z$classes), "weights.Omega", diagonal,
                              square = TRUE)
  if (!is.null(nu)) check_scalar(nu, "nu", 0, Inf)
  grid <- list(nrho = nrho, rho = if (!missing(rho)) rho, rho.min.ratio = if (!missing(rho.min.ratio)) rho.min.ratio,
               nlambda = if (q) nlambda else 1L, lambda = if (!missing(lambda)) lambda,
               lambda.min.ratio = if (!missing(lambda.min.ratio)) lambda.min.ratio)
  if (q == 0L && (!is.null(grid$lambda) || !is.null(grid$lambda.min.ratio)))
    stop(sQuote("lambda"), " cannot be used without predictors", call. = FALSE)
  out <- jcglasso.fit(Z, diagonal = diagonal, wB = wB, wTht = wTht, wOmg = wOmg, grid = grid, nu = nu,
                      alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, maxit.em = maxit.em, thr.em = thr.em,
                      maxit.bcd = maxit.bcd, thr.bcd = thr.bcd, trace = trace)
  if (any(out$conv != 0L))
    warning("the algorithm did not converge for ", sum(out$conv != 0L), " of the fitted models", call. = FALSE)
  out$call <- this.call
  out$Z <- Z
  out$weights.B <- wB; out$weights.Tht <- wTht; out$weights.Omega <- wOmg
  out$diagonal <- diagonal; out$penalty <- penalty
  out$alpha1 <- alpha1; out$alpha2 <- alpha2; out$alpha3 <- alpha3
  out$control <- list(maxit.em = maxit.em, thr.em = thr.em, maxit.bcd = maxit.bcd, thr.bcd = thr.bcd)
  out$nobs <- n; out$nresp <- p; out$npred <- q
  out$model <- model_name(Z, q)
  class(out) <- "jcglasso"
  out
}

model_name <- function(Z, q) {
  mar <- any(vapply(Z$R, function(r) any(r == 9L), TRUE))
  cens <- any(vapply(Z$R, function(r) any(r == 1L | r == -1L), TRUE))
  nm <- if (mar && cens) "joint hybrid glasso" else if (mar) "joint missglasso" else if (cens) "joint censored glasso"
        else "joint glasso"
  if (q > 0L) paste("conditional", nm) else nm
}

# Restrict a datajcggm object to the variables in the formula.
subset_by_formula <- function(data, formula) {
  if (!inherits(formula, "formula") || length(formula) != 3L) stop(sQuote("formula"), " must be a two-sided formula",
                                                                   call. = FALSE)
  ynames <- colnames(data$Y[[1L]]); xnames <- if (!is.null(data$X)) colnames(data$X[[1L]])
  lhs <- formula[[2L]]
  ysel <- if (identical(lhs, quote(.))) ynames else {
    if (is.call(lhs) && identical(lhs[[1L]], quote(cbind))) vapply(as.list(lhs)[-1L], deparse, "")
    else deparse(lhs)
  }
  if (anyDuplicated(ysel)) stop("repeated response variables in ", sQuote("formula"), call. = FALSE)
  if (any(!ysel %in% ynames)) stop("not a response variable of ", sQuote("data"), ": ",
                                   paste(setdiff(ysel, ynames), collapse = ", "), call. = FALSE)
  rhs <- formula[[3L]]
  if (identical(rhs, quote(.))) xsel <- xnames
  else {
    tt <- stats::terms(stats::as.formula(call("~", rhs)))
    xsel <- attr(tt, "term.labels")
    if (!attr(tt, "intercept")) stop("models without intercepts are not available", call. = FALSE)
    bad <- !vapply(xsel, function(v) make.names(v) == v || v %in% xnames, TRUE) | !xsel %in% xnames
    if (any(bad)) stop("the right-hand side of ", sQuote("formula"), " can only contain predictors of ",
                       sQuote("data"), ": ", paste(xsel[bad], collapse = ", "), call. = FALSE)
  }
  Z <- data
  ix <- match(c(ysel, xsel), c(ynames, xnames))
  Z$Y <- lapply(data$Y, function(y) y[, ysel, drop = FALSE])
  Z$X <- if (length(xsel)) lapply(data$X, function(x) x[, xsel, drop = FALSE]) else NULL
  Z$R <- lapply(data$R, function(r) r[, ix, drop = FALSE])
  Z$lo <- data$lo[ix]; Z$up <- data$up[ix]
  Z$zm <- data$zm[ix, , drop = FALSE]; Z$zv <- data$zv[ix, , drop = FALSE]
  Z
}

make_weights <- function(w, dims, dn, what, diagonal = FALSE, square = FALSE) {
  K <- dims[3L]
  if (is.null(w)) {
    w <- array(1, dims, dn)
    if (square && !diagonal) for (k in seq_len(K)) diag(w[, , k]) <- 0
    return(w)
  }
  if (is.list(w)) {
    if (length(w) != K) stop(sQuote(what), " must be a list of ", K, " matrices", call. = FALSE)
    w <- array(unlist(lapply(w, as.matrix)), dims)
  }
  if (is.matrix(w) && K == 1L) w <- array(w, dims)
  if (!is.array(w) || length(dim(w)) != 3L || any(dim(w) != dims)) stop(sQuote(what), " must have dimensions ", paste(dims, collapse = " x "),
                                                call. = FALSE)
  if (anyNA(w) || any(w < 0)) stop(sQuote(what), " must contain non-negative weights", call. = FALSE)
  if (square) for (k in seq_len(K)) {
    if (!isSymmetric(unname(w[, , k]))) stop(sQuote(what), " must contain symmetric matrices", call. = FALSE)
    if (!diagonal && any(diag(w[, , k]) != 0)) stop("diagonal = FALSE, but some diagonal entries of ", sQuote(what),
                                                    " are not zero", call. = FALSE)
  }
  dimnames(w) <- dn
  w
}

check_scalar <- function(x, what, lo, hi, open = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && !is.na(x) && (if (open) x > lo else x >= lo) && x <= hi
  if (!ok) stop(sQuote(what), " must be a number in ", if (open) "(" else "[", lo, ", ", hi, "]", call. = FALSE)
}
check_int <- function(x, what) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 1 || x != round(x))
    stop(sQuote(what), " must be a positive integer", call. = FALSE)
}
check_logical <- function(x, what) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) stop(sQuote(what), " must be TRUE or FALSE", call. = FALSE)
}
