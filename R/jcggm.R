#' Maximum Likelihood Refit of a jcglasso Model
#'
#' Computes the maximum likelihood estimates of the parameters of a joint
#' conditional Gaussian graphical model with the sparsity pattern of a model
#' fitted by [jcglasso()] (post-hoc refitting): the zero estimates are treated as
#' structural zeros and the non-zero ones are not penalized.
#'
#' @param object an object of class `"jcglasso"`.
#' @param GoF,lambda.id,rho.id used to select the model, see [select_jcglasso()].
#' @param tp.min the tiny value of the tuning parameters used for the refit
#'   (default `1e-6`).
#' @param trace integer: see [jcglasso()].
#' @param ... further arguments passed to [select_jcglasso()].
#' @details The refit uses the EM algorithm of [jcglasso()] with the support
#'   fixed and penalties equal to `tp.min` (needed for numerical stability when
#'   the number of non-zero parameters is close to the sample size).
#' @return An object of class `c("jcggm", "jcglasso")` with a single model.
#' @seealso [jcglasso()], [QFun()]
#' @examples
#' set.seed(1)
#' sim <- rjcggm(n = 50, p = 6, q = 2, K = 2, ncompPrecision = 2)
#' fit <- jcglasso(data = sim$Z, nrho = 3, nlambda = 3, rho.min.ratio = 0.3,
#'                 lambda.min.ratio = 0.3)
#' mle <- jcggm(fit, GoF = BIC)
#' coef(mle, type = "Theta")
#' @export
jcggm <- function(object, GoF = AIC, lambda.id, rho.id, tp.min = 1e-6, trace = 0L, ...) {
  if (!inherits(object, "jcglasso")) stop(sQuote("object"), " must be an object of class ", sQuote("jcglasso"),
                                          call. = FALSE)
  check_scalar(tp.min, "tp.min", 0, Inf, open = TRUE)
  this.call <- match.call()
  if (object$nlambda > 1L || object$nrho > 1L || !missing(lambda.id) || !missing(rho.id)) {
    args <- list(object = object, GoF = GoF, ...)
    if (!missing(lambda.id)) args$lambda.id <- lambda.id
    if (!missing(rho.id)) args$rho.id <- rho.id
    object <- do.call(select_jcglasso, args)
  }
  q <- object$npred; K <- length(object$nobs)
  wTht <- object$weights.Tht
  Tht <- object$Tht[, , , 1L, 1L, drop = FALSE]; dim(Tht) <- dim(Tht)[1:3]
  for (k in seq_len(K)) {
    z <- Tht[, , k] == 0
    wTht[, , k][z] <- Inf
    wTht[, , k][!z] <- 0
  }
  wB <- wOmg <- NULL; B <- Omg <- NULL
  if (q) {
    B <- object$B[-1L, , , 1L, 1L, drop = FALSE]; dim(B) <- dim(B)[1:3]
    Omg <- object$Omega[, , , 1L, 1L, drop = FALSE]; dim(Omg) <- dim(Omg)[1:3]
    wB <- object$weights.B; wOmg <- object$weights.Omega
    for (k in seq_len(K)) {
      wB[, , k][B[, , k] == 0] <- Inf; wB[, , k][B[, , k] != 0] <- 0
      z <- Omg[, , k] == 0
      wOmg[, , k][z] <- Inf; wOmg[, , k][!z] <- 0
    }
    dimnames(B) <- dimnames(object$weights.B); dimnames(Omg) <- dimnames(object$weights.Omega)
  }
  dimnames(Tht) <- dimnames(object$weights.Tht)
  psi <- object$psi[, , 1L, 1L]; dim(psi) <- dim(object$psi)[1:2]; dimnames(psi) <- dimnames(object$psi)[1:2]
  grid <- list(rho = tp.min, nrho = 1L, lambda = if (q) tp.min, nlambda = 1L)
  out <- jcglasso.fit(object$Z, diagonal = FALSE, wB = wB, wTht = wTht, wOmg = wOmg, grid = grid,
                      nu = if (q) tp.min, alpha1 = 1, alpha2 = 1, alpha3 = 1,
                      maxit.em = object$control$maxit.em, thr.em = object$control$thr.em,
                      maxit.bcd = object$control$maxit.bcd, thr.bcd = object$control$thr.bcd, trace = trace,
                      start = list(psi = psi, Tht = Tht, B = B, Omega = Omg))
  res <- object
  for (nm in names(out)) res[nm] <- list(out[[nm]])
  res$lambda <- object$lambda; res$rho <- object$rho; res$nu <- object$nu
  res$refit <- list(tp.min = tp.min, weights.B = wB, weights.Tht = wTht, weights.Omega = wOmg)
  res$call <- this.call
  res$model <- paste("MLE of the", object$model)
  if (any(out$conv != 0L)) warning("the algorithm did not converge", call. = FALSE)
  class(res) <- c("jcggm", "jcglasso")
  res
}
