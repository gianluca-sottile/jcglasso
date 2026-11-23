# Reference implementations in R of the numerical kernels (used for testing the C
# code and as a fallback). Notation: K classes, q predictors, p responses,
# f = class weights n_k / n.

# Group weights: for each entry, the mean of the finite weights across the classes
# (structural zeros, i.e. weights >= BIG, are excluded; they are enforced by the lasso part).
group_weights <- function(w, big = 1e300) {
  W <- w; W[W >= big] <- NA
  gw <- apply(W, c(1L, 2L), mean, na.rm = TRUE)
  gw[is.na(gw)] <- big
  gw
}

# Sparse group lasso proximal operator on an array x (q x p x K): element-wise soft
# thresholding with thresholds a * w, followed by group shrinkage of each (h, j) group
# across the K classes with thresholds g * gw (gw: matrix of group weights).
prox_sgl_R <- function(x, a, g, w, gw) {
  x <- sign(x) * pmax(abs(x) - a * w, 0)
  if (dim(x)[3L] > 1L && any(g > 0)) {
    nrm <- sqrt(apply(x^2, c(1L, 2L), sum))
    fac <- pmax(1 - g * gw / pmax(nrm, .Machine$double.xmin), 0)
    fac[nrm == 0] <- 0
    x <- x * array(fac, dim(x))
  }
  x
}

# Gradient of (1/2) sum_k f_k tr(Tht_k (B_k' Sxx_k B_k - 2 B_k' Sxy_k)) w.r.t. B.
grad_B_R <- function(B, Sxx, Sxy, Tht, f) {
  G <- B
  for (k in seq_along(f)) G[, , k] <- f[k] * (Sxx[, , k] %*% B[, , k] - Sxy[, , k]) %*% Tht[, , k]
  G
}

# Objective of the B sub-problem.
obj_B_R <- function(B, Sxx, Sxy, Tht, f, lambda, alpha, w) {
  val <- 0
  for (k in seq_along(f)) {
    Bk <- B[, , k]
    val <- val + 0.5 * f[k] * sum(Tht[, , k] * (crossprod(Bk, Sxx[, , k] %*% Bk) - 2 * crossprod(Bk, Sxy[, , k])))
  }
  pen <- alpha * sum(w * abs(B))
  if (dim(B)[3L] > 1L) pen <- pen + (1 - alpha) * sum(group_weights(w) * sqrt(apply(B^2, c(1L, 2L), sum)))
  val + lambda * pen
}

# Accelerated proximal gradient (Barzilai-Borwein initial step, adaptive restart,
# TFOCS-style step update), as in the Fortran subroutine 'apg'.
apg_B_R <- function(B, Sxx, Sxy, Tht, f, lambda, alpha, w, maxit = 1e4, thr = 1e-6, gw = group_weights(w)) {
  x <- B; y <- x
  g <- grad_B_R(y, Sxx, Sxy, Tht, f)
  nrmg <- sqrt(sum(g^2))
  if (nrmg == 0) return(list(B = x, nit = 0L, conv = 0L))
  t <- 1 / nrmg
  xh <- x - t * g
  gh <- grad_B_R(xh, Sxx, Sxy, Tht, f)
  den <- sum((g - gh)^2)
  t <- if (den > 0) abs(sum((x - xh) * (g - gh)) / den) else 1 / nrmg
  scale <- max(nrmg, lambda)
  theta <- 1
  conv <- 1L
  for (it in seq_len(maxit)) {
    xold <- x; yold <- y
    x <- prox_sgl_R(y - t * g, t * lambda * alpha, t * lambda * (1 - alpha), w, gw)
    # stopping rule on the gradient mapping, relative to the scale of the problem
    dn <- sqrt(sum((y - x)^2))
    if (dn / t <= thr * scale || dn <= 1e-14 * max(1, sqrt(sum(x^2)))) { conv <- 0L; break }
    theta <- 2 / (1 + sqrt(1 + 4 / theta^2))
    if (sum((y - x) * (x - xold)) > 0) { x <- xold; y <- x; theta <- 1 }
    else y <- x + (1 - theta) * (x - xold)
    gold <- g
    g <- grad_B_R(y, Sxx, Sxy, Tht, f)
    dy <- y - yold
    den <- abs(sum(dy * (gold - g)))
    that <- if (den > 0) 0.5 * sum(dy^2) / den else t
    t <- min(1.01 * t, max(0.5 * t, that))
  }
  list(B = x, nit = it, conv = conv)
}

# ADMM for the B sub-problem (Sottile et al., 2025, Section 5.2): the B-update solves
# f_k Sxx_k B Tht_k + tau B = f_k Sxy_k Tht_k + tau (Gamma - U) with the eigen-decompositions
# of Sxx_k and Tht_k.
admm_B_R <- function(B, Sxx, Sxy, Tht, f, lambda, alpha, w, tau = 2, maxit = 1e4, thr = 1e-6, gw = group_weights(w)) {
  K <- length(f)
  ex <- lapply(seq_len(K), function(k) eigen(Sxx[, , k], symmetric = TRUE))
  et <- lapply(seq_len(K), function(k) eigen(Tht[, , k], symmetric = TRUE))
  C <- B
  for (k in seq_len(K)) C[, , k] <- f[k] * Sxy[, , k] %*% Tht[, , k]
  G <- B; U <- 0 * B
  conv <- 1L
  for (it in seq_len(maxit)) {
    for (k in seq_len(K)) {
      Vx <- ex[[k]]$vectors; Vt <- et[[k]]$vectors
      rhs <- crossprod(Vx, C[, , k] + tau * (G[, , k] - U[, , k])) %*% Vt
      B[, , k] <- Vx %*% (rhs / (f[k] * outer(ex[[k]]$values, et[[k]]$values) + tau)) %*% t(Vt)
    }
    Gold <- G
    G <- prox_sgl_R(B + U, lambda * alpha / tau, lambda * (1 - alpha) / tau, w, gw)
    U <- U + B - G
    r <- sqrt(sum((B - G)^2)); s <- tau * sqrt(sum((G - Gold)^2))
    if (r <= thr * max(1, sqrt(sum(G^2))) && s <= thr * max(1, tau * sqrt(sum(U^2)))) { conv <- 0L; break }
  }
  list(B = G, nit = it, conv = conv)
}

# Joint graphical lasso with group penalty (Danaher et al., 2014):
#   min sum_k f_k (tr(S_k Tht_k) - log det Tht_k)
#       + rho * [alpha sum_k sum_{i != j} w_ijk |tht_ijk| + (1 - alpha) sum_{i != j} gw_ij (sum_k tht_ijk^2)^(1/2)]
# solved by ADMM with adaptive penalty parameter; diagonal entries unpenalized unless
# 'pendiag = TRUE' (then they are penalized with the lasso part only).
# Before the ADMM, the variables are split into connected components (Danaher et al.,
# 2014, Theorem 2): i and j are disconnected if || (f_k |S_ijk| - alpha rho w_ijk)_+ ||_2 <=
# (1 - alpha) rho gw_ij.
jgl_components_R <- function(S, f, rho, alpha, w, gw) {
  p <- dim(S)[1L]; K <- dim(S)[3L]
  A <- matrix(0, p, p)
  for (k in seq_len(K)) A <- A + pmax(f[k] * abs(S[, , k]) - alpha * rho * w[, , k], 0)^2
  adj <- sqrt(A) > (1 - alpha) * rho * gw
  diag(adj) <- FALSE
  comp <- integer(p); nc <- 0L
  for (v in seq_len(p)) if (comp[v] == 0L) {
    nc <- nc + 1L; stack <- v; comp[v] <- nc
    while (length(stack)) {
      u <- stack[1L]; stack <- stack[-1L]
      nb <- which(adj[u, ] & comp == 0L)
      comp[nb] <- nc; stack <- c(stack, nb)
    }
  }
  comp
}

jgl_admm_block_R <- function(S, f, rho, alpha, w, gw, Tht, pendiag = FALSE, maxit = 1e4, thr = 1e-5, tau = 2) {
  p <- dim(S)[1L]; K <- dim(S)[3L]
  Z <- Tht; U <- 0 * Tht
  offd <- !diag(p)
  conv <- 1L
  for (it in seq_len(maxit)) {
    Zold <- Z
    for (k in seq_len(K)) {
      e <- eigen(S[, , k] - (tau * Z[, , k] - U[, , k]) / f[k], symmetric = TRUE)
      d <- f[k] * (-e$values + sqrt(e$values^2 + 4 * tau / f[k])) / (2 * tau)
      Tht[, , k] <- e$vectors %*% (d * t(e$vectors))
    }
    A <- Tht + U / tau
    lam1 <- alpha * rho * w / tau
    if (!pendiag) for (k in seq_len(K)) diag(lam1[, , k]) <- 0
    gw0 <- gw; diag(gw0) <- 0           # the group penalty does not involve the diagonal
    Z <- prox_sgl_R(A, lam1, (1 - alpha) * rho / tau, 1, gw0)
    U <- U + tau * (Tht - Z)
    r <- sum(f * sapply(seq_len(K), function(k) sqrt(sum((Tht[, , k] - Z[, , k])^2))))
    s <- sum(f * sapply(seq_len(K), function(k) sqrt(sum((tau * (Z[, , k] - Zold[, , k]))^2))))
    e1 <- sum(f * sapply(seq_len(K), function(k) p * thr + thr * max(sqrt(sum(Tht[, , k]^2)), sqrt(sum(Z[, , k]^2)))))
    e2 <- sum(f * sapply(seq_len(K), function(k) p * thr + thr * sqrt(sum(U[, , k]^2))))
    if (r < e1 && s < e2) { conv <- 0L; break }
    if (r > 10 * s) tau <- 2 * tau else if (s > 10 * r) tau <- tau / 2
  }
  list(Tht = Z, nit = it, conv = conv)
}

jgl_R <- function(S, f, rho, alpha, w, Tht, pendiag = FALSE, maxit = 1e4, thr = 1e-5, gw = group_weights(w)) {
  p <- dim(S)[1L]; K <- dim(S)[3L]
  comp <- jgl_components_R(S, f, rho, alpha, w, gw)
  out <- array(0, dim(S), dimnames(S))
  nit <- 0L; conv <- 0L
  for (cc in unique(comp)) {
    id <- which(comp == cc)
    if (length(id) == 1L) {
      for (k in seq_len(K)) {
        s <- S[id, id, k] + if (pendiag) alpha * rho * w[id, id, k] / f[k] else 0
        out[id, id, k] <- 1 / s
      }
    } else {
      r <- jgl_admm_block_R(S[id, id, , drop = FALSE], f, rho, alpha, w[id, id, , drop = FALSE],
                            gw[id, id, drop = FALSE], Tht[id, id, , drop = FALSE], pendiag, maxit, thr)
      out[id, id, ] <- r$Tht
      nit <- nit + r$nit; conv <- max(conv, r$conv)
    }
  }
  list(Tht = out, comp = comp, nit = nit, conv = conv)
}

# E-step for one class: conditional means of the unobserved entries given the observed
# ones (truncated normal for censored entries, normal for missing-at-random entries) and
# the expected second moments, with the approximation of Guo et al. (2015).
#   Z: n x d data (values of censored entries are irrelevant; NA for MAR)
#   R: n x d codes (0 observed, -1 left-censored, 1 right-censored, 9 missing at random)
#   psi: d-vector of means; Xi: d x d precision matrix; lo, up: censoring limits
#   zlo, zup: bounds used to clip the imputed MAR values
estep_R <- function(Z, R, psi, Xi, lo, up, zlo, zup, pat = make_patterns(R)) {
  n <- nrow(Z); d <- ncol(Z)
  Zi <- Z; Zi[R != 0] <- 0
  addvar <- matrix(0, n, d)
  for (h in seq_len(length(pat$start) - 1L)) {
    rows <- pat$row[(pat$start[h] + 1L):pat$start[h + 1L]] + 1L
    r <- R[rows[1L], ]
    m <- which(r != 0); o <- which(r == 0)
    Sm <- solve(Xi[m, m, drop = FALSE])
    sd <- sqrt(diag(Sm))
    M <- matrix(psi[m], length(rows), length(m), byrow = TRUE)
    if (length(o)) {
      D <- sweep(Z[rows, o, drop = FALSE], 2L, psi[o])
      M <- M - D %*% t(Sm %*% Xi[m, o, drop = FALSE])
    }
    for (jj in seq_along(m)) {
      j <- m[jj]; mu <- M[, jj]; s <- sd[jj]
      if (r[j] == 9L) {
        Zi[rows, j] <- pmin(pmax(mu, zlo[j]), zup[j]); addvar[rows, j] <- s^2
      } else {
        if (r[j] == -1L) { z <- (lo[j] - mu) / s; lr <- stats::dnorm(z, log = TRUE) - stats::pnorm(z, log.p = TRUE); sg <- -1 }
        else { z <- (up[j] - mu) / s; lr <- stats::dnorm(z, log = TRUE) - stats::pnorm(z, lower.tail = FALSE, log.p = TRUE); sg <- 1 }
        ratio <- exp(lr)
        tm <- mu + sg * s * ratio
        tv <- s^2 * (1 + sg * z * ratio - ratio^2)
        # as in cglasso: extreme imputations are clipped at the bounds zlo / zup (by
        # default, marginal mean -/+ 3 standard deviations) with the conditional variance
        bnd <- if (r[j] == -1L) min(zlo[j], lo[j]) else max(zup[j], up[j])
        bad <- !is.finite(tm) | !is.finite(tv) | tv <= 0 | (if (r[j] == -1L) tm < bnd else tm > bnd)
        tm[bad] <- bnd
        tv[bad] <- s^2
        Zi[rows, j] <- tm; addvar[rows, j] <- tv
      }
    }
  }
  list(Z = Zi, T1 = colSums(Zi), T2 = crossprod(Zi) + diag(colSums(addvar), d), addvar = addvar)
}
