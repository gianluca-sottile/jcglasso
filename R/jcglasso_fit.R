# EM algorithm for the joint conditional graphical lasso (internal).

BIG <- 1e300   # weights equal to Inf (structural zeros)

# Joint precision matrix of Z = (Y, X): [[Tht, -Tht B'], [-B Tht, Omega + B Tht B']].
build_Xi <- function(Tht, B, Omg) {
  if (is.null(B)) return(Tht)
  TB <- -Tht %*% t(B)
  rbind(cbind(Tht, TB), cbind(t(TB), Omg + B %*% Tht %*% t(B)))
}

# Conditional covariance of Y given X: Syy - Syx B - B' Sxy + B' Sxx B.
S_cond <- function(S, idY, idX, B) {
  if (is.null(B)) return(S[idY, idY, drop = FALSE])
  SxyB <- crossprod(B, S[idX, idY, drop = FALSE])
  S[idY, idY, drop = FALSE] - SxyB - t(SxyB) + crossprod(B, S[idX, idX, drop = FALSE] %*% B)
}

Qvalue <- function(Tht, S, n) {
  ld <- determinant(Tht, logarithm = TRUE)
  ld <- if (ld$sign > 0) as.numeric(ld$modulus) else -Inf
  n / 2 * (ld - sum(Tht * S)) - n * nrow(Tht) / 2 * log(2 * pi)
}

jcglasso.fit <- function(Z, diagonal, wB, wTht, wOmg, grid, nu, alpha1, alpha2, alpha3,
                         maxit.em, thr.em, maxit.bcd, thr.bcd, trace,
                         start = NULL, fixed.grid = FALSE) {
  K <- length(Z$classes); n <- nobs(Z); f <- n / sum(n)
  p <- ncol(Z$Y[[1L]]); q <- if (is.null(Z$X)) 0L else ncol(Z$X[[1L]]); d <- p + q
  idY <- seq_len(p); idX <- if (q) p + seq_len(q)
  ynames <- colnames(Z$Y[[1L]]); xnames <- if (q) colnames(Z$X[[1L]]); znames <- c(ynames, xnames)
  if (K == 1L) alpha1 <- alpha2 <- alpha3 <- 1
  wTht[is.infinite(wTht)] <- BIG
  if (q) { wB[is.infinite(wB)] <- BIG; wOmg[is.infinite(wOmg)] <- BIG }

  Zm <- lapply(seq_len(K), function(k) if (q) cbind(Z$Y[[k]], Z$X[[k]]) else Z$Y[[k]])
  zlo <- Z$zm - 3 * sqrt(Z$zv); zup <- Z$zm + 3 * sqrt(Z$zv)
  pats <- lapply(Z$R, make_patterns)

  # starting values
  psi <- Z$zm
  Tht <- array(0, c(p, p, K), list(ynames, ynames, Z$classes))
  for (k in seq_len(K)) Tht[, , k] <- diag(1 / Z$zv[idY, k], p)
  B <- if (q) array(0, c(q, p, K), list(xnames, ynames, Z$classes))
  Omg <- if (q) { o <- array(0, c(q, q, K), list(xnames, xnames, Z$classes))
                  for (k in seq_len(K)) o[, , k] <- diag(1 / Z$zv[idX, k], q); o }
  if (!is.null(start)) { Tht <- start$Tht; if (q) { B <- start$B; Omg <- start$Omega }; psi <- start$psi }

  estep_all <- function(psi, Tht, B, Omg) {
    S <- array(0, c(d, d, K)); psi_new <- psi; Zi <- vector("list", K)
    for (k in seq_len(K)) {
      Xi <- build_Xi(Tht[, , k], if (q) B[, , k], if (q) Omg[, , k])
      es <- k_estep(Zm[[k]], Z$R[[k]], psi[, k], Xi, Z$lo, Z$up, zlo[, k], zup[, k], pats[[k]])
      psi_new[, k] <- es$T1 / n[k]
      S[, , k] <- es$T2 / n[k] - tcrossprod(psi_new[, k])
      Zi[[k]] <- es$Z
    }
    list(psi = psi_new, S = S, Z = Zi)
  }

  # grids of tuning parameters (from the starting E-step)
  E0 <- estep_all(psi, Tht, B, Omg)
  S0 <- E0$S
  if (!is.null(grid$rho)) {
    rho <- sort(grid$rho, decreasing = TRUE)
    if (any(rho < 0)) stop("negative values in ", sQuote("rho"), call. = FALSE)
  } else {
    U <- upper.tri(diag(p))
    g <- vapply(seq_len(K), function(k) f[k] * abs(S0[idY, idY, k][U]), numeric(sum(U)))
    g <- matrix(g, ncol = K)
    rho_max <- if (alpha1 < 1) max(sqrt(rowSums(g^2))) else max(g)
    rmr <- if (!is.null(grid$rho.min.ratio)) grid$rho.min.ratio else if (all(p < n)) 1e-2 else 5e-2
    check_scalar(rmr, "rho.min.ratio", 0, 1)
    check_int(grid$nrho, "nrho")
    rho <- if (grid$nrho == 1L) rho_max else exp(seq(log(rho_max), log(max(rho_max * rmr, 1e-10)),
                                                     length.out = grid$nrho))
  }
  nrho <- length(rho)
  if (q) {
    if (!is.null(grid$lambda)) {
      lambda <- sort(grid$lambda, decreasing = TRUE)
      if (any(lambda < 0)) stop("negative values in ", sQuote("lambda"), call. = FALSE)
    } else {
      G <- vapply(seq_len(K), function(k) f[k] * abs(S0[idX, idY, k] %*% Tht[, , k]), matrix(0, q, p))
      G <- array(G, c(q, p, K))
      lambda_max <- if (alpha2 < 1) max(sqrt(apply(G^2, c(1L, 2L), sum))) else max(G)
      lmr <- if (!is.null(grid$lambda.min.ratio)) grid$lambda.min.ratio else if (all(q < n)) 1e-2 else 5e-2
      check_scalar(lmr, "lambda.min.ratio", 0, 1)
      check_int(grid$nlambda, "nlambda")
      lambda <- if (grid$nlambda == 1L) lambda_max else exp(seq(log(lambda_max), log(max(lambda_max * lmr, 1e-10)),
                                                                length.out = grid$nlambda))
    }
  } else lambda <- 0
  nlambda <- length(lambda)
  nuv <- if (q) { if (is.null(nu)) rho else rep(nu, nrho) } else rep(0, nrho)
  # structural zeros need a positive penalty
  eps_pen <- 1e-12
  # the M-step must be solved more accurately than the EM tolerance, otherwise the EM
  # algorithm stalls at the accuracy of the inner solvers
  thr.in <- min(thr.bcd, thr.em / 10)
  lab_l <- paste0("lambda_", seq_len(nlambda)); lab_r <- paste0("rho_", seq_len(nrho))

  # output
  Bout <- if (q) array(0, c(q + 1L, p, K, nlambda, nrho), list(c("(Intercept)", xnames), ynames, Z$classes, lab_l, lab_r))
          else array(0, c(1L, p, K, nlambda, nrho), list("(Intercept)", ynames, Z$classes, lab_l, lab_r))
  ThtOut <- array(0, c(p, p, K, nlambda, nrho), list(ynames, ynames, Z$classes, lab_l, lab_r))
  OmgOut <- if (q) array(0, c(q, q, K, nlambda, nrho), list(xnames, xnames, Z$classes, lab_l, lab_r))
  psiOut <- array(0, c(d, K, nlambda, nrho), list(znames, Z$classes, lab_l, lab_r))
  Qyx <- Qx <- dfB <- dfTht <- dfOmg <- array(0, c(K, nlambda, nrho), list(Z$classes, lab_l, lab_r))
  nit <- array(0L, c(2L, nlambda, nrho), list(c("EM", "M-step"), lab_l, lab_r))
  conv <- matrix(0L, nlambda, nrho, dimnames = list(lab_l, lab_r))

  cur <- list(psi = psi, Tht = Tht, B = B, Omg = Omg)
  first_of_row <- cur
  for (h in seq_len(nlambda)) {
    if (h > 1L) cur <- first_of_row
    for (o in seq_len(nrho)) {
      lam <- max(lambda[h], if (q && any(wB >= BIG)) eps_pen else 0)
      rh <- max(rho[o], if (any(wTht >= BIG)) eps_pen else 0)
      nuo <- max(nuv[o], if (q && any(wOmg >= BIG)) eps_pen else 0)
      psi <- cur$psi; Tht <- cur$Tht; B <- cur$B; Omg <- cur$Omg
      ninner <- 0L; code <- 1L
      for (it in seq_len(maxit.em)) {
        Tht_o <- Tht; B_o <- B; Omg_o <- Omg
        E <- estep_all(psi, Tht, B, Omg)
        psi <- E$psi; S <- E$S
        if (q) {
          r <- k_jgl(S[idX, idX, , drop = FALSE], f, nuo, alpha3, wOmg, Omg, diagonal, maxit.bcd, thr.in)
          Omg[] <- r$Tht; ninner <- ninner + r$nit
          r <- k_apg(B, S[idX, idX, , drop = FALSE], S[idX, idY, , drop = FALSE], Tht, f, lam, alpha2, wB,
                     maxit.bcd, thr.in)
          B[] <- r$B; ninner <- ninner + r$nit
        }
        Sc <- array(0, c(p, p, K))
        for (k in seq_len(K)) Sc[, , k] <- S_cond(S[, , k], idY, idX, if (q) B[, , k])
        r <- k_jgl(Sc, f, rh, alpha1, wTht, Tht, diagonal, maxit.bcd, thr.in)
        Tht[] <- r$Tht; ninner <- ninner + r$nit
        dlt <- max(rel_change(Tht, Tht_o, f), if (q) rel_change(B, B_o, f) else 0,
                   if (q) rel_change(Omg, Omg_o, f) else 0)
        if (trace == 2L) cat(sprintf("    EM step %4d: max relative change = %.3e\n", it, dlt))
        if (dlt <= thr.em) { code <- 0L; break }
      }
      # store the solution
      for (k in seq_len(K)) {
        if (q) {
          Bout[-1L, , k, h, o] <- B[, , k]
          Bout[1L, , k, h, o] <- psi[idY, k] - drop(crossprod(B[, , k], psi[idX, k]))
          OmgOut[, , k, h, o] <- Omg[, , k]
          dfB[k, h, o] <- sum(B[, , k] != 0)
          dfOmg[k, h, o] <- sum(Omg[, , k][upper.tri(Omg[, , k])] != 0)
          Qx[k, h, o] <- Qvalue(Omg[, , k], S[idX, idX, k], n[k])
        } else Bout[1L, , k, h, o] <- psi[idY, k]
        ThtOut[, , k, h, o] <- Tht[, , k]
        dfTht[k, h, o] <- sum(Tht[, , k][upper.tri(Tht[, , k])] != 0)
        Qyx[k, h, o] <- Qvalue(Tht[, , k], Sc[, , k], n[k])
      }
      psiOut[, , h, o] <- psi
      nit[, h, o] <- c(it, ninner)
      conv[h, o] <- code
      if (trace >= 1L) cat(sprintf("model (lambda %d, rho %d): lambda = %.4g, rho = %.4g, nu = %.4g, EM steps = %d%s\n",
                                   h, o, lambda[h], rho[o], nuv[o], it, if (code) " (not converged)" else ""))
      cur <- list(psi = psi, Tht = Tht, B = B, Omg = Omg)
      if (o == 1L) first_of_row <- cur
    }
  }
  list(B = Bout, Tht = ThtOut, Omega = OmgOut, psi = psiOut, Q = list(yx = Qyx, x = if (q) Qx),
       dfB = if (q) dfB, dfTht = dfTht, dfOmg = if (q) dfOmg, nit = nit, conv = conv,
       lambda = lambda, rho = rho, nu = nuv, nlambda = nlambda, nrho = nrho)
}

rel_change <- function(new, old, f) {
  K <- length(f)
  sum(f * vapply(seq_len(K), function(k) {
    a <- new[, , k]; b <- old[, , k]
    sqrt(sum((a - b)^2)) / max(sqrt(sum(a^2)), 1e-10)
  }, 0))
}
