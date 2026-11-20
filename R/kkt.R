#' Optimality Conditions of a Fitted jcglasso Model
#'
#' Checks the Karush-Kuhn-Tucker conditions of the M-step at a fitted model:
#' given the expected sufficient statistics computed by an E-step at the
#' solution, returns the largest violation of the conditions for \eqn{B},
#' \eqn{\Theta} and \eqn{\Omega}. At convergence of the EM algorithm the
#' violations are small (of the order of the tolerances `thr.em` and `thr.bcd`).
#'
#' @param object an object of class `"jcglasso"`.
#' @param lambda.id,rho.id indices of the model.
#' @return A named vector with the largest violations for `B`, `Theta` and `Omega`.
#' @seealso [jcglasso()]
#' @examples
#' set.seed(1)
#' sim <- rjcggm(n = 50, p = 6, q = 2, K = 2, ncompPrecision = 2)
#' fit <- jcglasso(data = sim$Z, nrho = 3, nlambda = 3, rho.min.ratio = 0.3,
#'                 lambda.min.ratio = 0.3)
#' kkt(fit, lambda.id = 3, rho.id = 3)
#' @export
kkt <- function(object, lambda.id, rho.id) {
  if (!inherits(object, "jcglasso")) stop(sQuote("object"), " must be an object of class ", sQuote("jcglasso"),
                                          call. = FALSE)
  ids <- path_ids(object, lambda.id, rho.id)
  h <- ids[1L]; o <- ids[2L]
  Z <- object$Z; K <- length(Z$classes); n <- object$nobs; f <- n / sum(n)
  p <- object$nresp; q <- object$npred; idY <- seq_len(p); idX <- if (q) p + seq_len(q)
  a1 <- if (K == 1L) 1 else object$alpha1; a2 <- if (K == 1L) 1 else object$alpha2
  a3 <- if (K == 1L) 1 else object$alpha3
  Tht <- object$Tht[, , , h, o, drop = FALSE]; dim(Tht) <- dim(Tht)[1:3]
  B <- if (q) { b <- object$B[-1L, , , h, o, drop = FALSE]; dim(b) <- dim(b)[1:3]; b }
  Omg <- if (q) { w <- object$Omega[, , , h, o, drop = FALSE]; dim(w) <- dim(w)[1:3]; w }
  S <- array(0, c(p + q, p + q, K))
  for (k in seq_len(K)) {
    Zk <- if (q) cbind(Z$Y[[k]], Z$X[[k]]) else Z$Y[[k]]
    Xi <- build_Xi(Tht[, , k], if (q) B[, , k], if (q) Omg[, , k])
    es <- k_estep(Zk, Z$R[[k]], object$psi[, k, h, o], Xi, Z$lo, Z$up, Z$zm[, k] - 3 * sqrt(Z$zv[, k]),
                  Z$zm[, k] + 3 * sqrt(Z$zv[, k]))
    m <- es$T1 / n[k]
    S[, , k] <- es$T2 / n[k] - tcrossprod(m)
  }
  wT <- object$weights.Tht; wB <- object$weights.B; wO <- object$weights.Omega
  lam <- object$lambda[h]; rho <- object$rho[o]; nu <- object$nu[o]
  if (inherits(object, "jcggm")) {
    wT <- object$refit$weights.Tht; wB <- object$refit$weights.B; wO <- object$refit$weights.Omega
    lam <- rho <- nu <- object$refit$tp.min; a1 <- a2 <- a3 <- 1
  }
  Sc <- array(0, c(p, p, K))
  for (k in seq_len(K)) Sc[, , k] <- S_cond(S[, , k], idY, idX, if (q) B[, , k])
  vT <- kkt_prec(Tht, Sc, f, rho, a1, wT, object$diagonal)
  vB <- vO <- NA_real_
  if (q) {
    G <- grad_B_R(B, S[idX, idX, , drop = FALSE], S[idX, idY, , drop = FALSE], Tht, f)
    vB <- kkt_sgl(B, G, lam, a2, wB)
    vO <- kkt_prec(Omg, S[idX, idX, , drop = FALSE], f, nu, a3, wO, object$diagonal)
  }
  c(B = vB, Theta = vT, Omega = vO)
}

# Largest violation of the optimality conditions of sum_k g_k(x_k) + lam * SGL(x),
# where G is the gradient of the smooth part.
kkt_sgl <- function(x, G, lam, alpha, w, offdiag = FALSE) {
  K <- dim(x)[3L]
  nrm <- sqrt(apply(x^2, c(1L, 2L), sum))
  viol <- 0
  for (i in seq_len(dim(x)[1L])) for (j in seq_len(dim(x)[2L])) {
    if (offdiag && i == j) next
    wij <- w[i, j, ]
    if (any(is.infinite(wij))) next
    if (nrm[i, j] > 0) {
      for (k in seq_len(K)) if (x[i, j, k] != 0)
        viol <- max(viol, abs(G[i, j, k] + lam * alpha * wij[k] * sign(x[i, j, k]) +
                                lam * (1 - alpha) * wij[k] * x[i, j, k] / nrm[i, j]))
    } else {
      ex <- sqrt(sum(pmax(abs(G[i, j, ]) - lam * alpha * wij, 0)^2)) - lam * (1 - alpha) * max(wij)
      viol <- max(viol, ex)
    }
  }
  viol
}

kkt_prec <- function(Tht, S, f, rho, alpha, w, diagonal) {
  K <- dim(Tht)[3L]
  G <- Tht
  for (k in seq_len(K)) G[, , k] <- f[k] * (S[, , k] - solve(Tht[, , k]))
  v <- kkt_sgl(Tht, G, rho, alpha, w, offdiag = TRUE)
  if (!diagonal) v <- max(v, max(abs(vapply(seq_len(K), function(k) diag(G[, , k]), numeric(dim(Tht)[1L])))))
  v
}
