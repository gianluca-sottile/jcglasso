test_that("E-step moments of a censored bivariate normal", {
  Sg <- matrix(c(1, 0.6, 0.6, 1), 2); Xi <- solve(Sg); psi <- c(0, 0)
  Z <- matrix(c(1.5, 0.3), 1)        # first variable right censored at 1.5
  R <- matrix(c(1L, 0L), 1)
  es <- ns$estep_R(Z, R, psi, Xi, c(-1e300, -1e300), c(1.5, 1e300), c(-10, -10), c(10, 10))
  m <- 0.6 * 0.3; s <- sqrt(1 - 0.36); a <- (1.5 - m) / s
  expect_equal(es$Z[1, 1], m + s * dnorm(a) / pnorm(a, lower.tail = FALSE), tolerance = 1e-12)
  r <- dnorm(a) / pnorm(a, lower.tail = FALSE)
  expect_equal(es$T2[1, 1] - es$Z[1, 1]^2, s^2 * (1 + a * r - r^2), tolerance = 1e-10)
})

test_that("APG and ADMM solve the same B sub-problem", {
  set.seed(3)
  q <- 4; p <- 3; K <- 2; f <- c(0.4, 0.6)
  Sxx <- Sxy <- Tht <- NULL
  Sxx <- array(0, c(q, q, K)); Sxy <- array(0, c(q, p, K)); Tht <- array(0, c(p, p, K))
  for (k in 1:K) {
    X <- matrix(rnorm(50 * q), 50); Y <- X %*% matrix(rnorm(q * p), q) + matrix(rnorm(50 * p), 50)
    Sxx[, , k] <- crossprod(X) / 50; Sxy[, , k] <- crossprod(X, Y) / 50; Tht[, , k] <- solve(cov(Y))
  }
  w <- array(1, c(q, p, K))
  a <- ns$apg_B_R(array(0, c(q, p, K)), Sxx, Sxy, Tht, f, 0.05, 0.5, w, thr = 1e-10)
  b <- ns$admm_B_R(array(0, c(q, p, K)), Sxx, Sxy, Tht, f, 0.05, 0.5, w, thr = 1e-10)
  expect_equal(a$B, b$B, tolerance = 1e-6)
  G <- ns$grad_B_R(a$B, Sxx, Sxy, Tht, f)
  expect_lt(ns$kkt_sgl(a$B, G, 0.05, 0.5, w), 1e-6)
})

test_that("joint graphical lasso satisfies the optimality conditions", {
  set.seed(4)
  p <- 6; K <- 2; f <- c(0.5, 0.5); S <- array(0, c(p, p, K))
  for (k in 1:K) S[, , k] <- cov(matrix(rnorm(40 * p), 40) %*% chol(0.5^abs(outer(1:p, 1:p, "-"))))
  w <- array(1, c(p, p, K)); for (k in 1:K) diag(w[, , k]) <- 0
  T0 <- array(0, c(p, p, K)); for (k in 1:K) T0[, , k] <- diag(1 / diag(S[, , k]))
  for (alpha in c(0, 0.5, 1)) {
    r <- ns$jgl_R(S, f, 0.05, alpha, w, T0, FALSE, 1e4, 1e-9)
    expect_lt(ns$kkt_prec(r$Tht, S, f, 0.05, alpha, w, FALSE), 1e-5)
  }
})

test_that("compiled kernels agree with the R implementations", {
  skip_if_not(has_C)
  sim <- small_sim()
  old <- options(jcglasso.engine = "R"); on.exit(options(old))
  f1 <- jcglasso(data = sim$Z, nrho = 2, nlambda = 2, rho.min.ratio = 0.3, lambda.min.ratio = 0.3)
  options(jcglasso.engine = "C")
  f2 <- jcglasso(data = sim$Z, nrho = 2, nlambda = 2, rho.min.ratio = 0.3, lambda.min.ratio = 0.3)
  expect_equal(f1$B, f2$B, tolerance = 1e-8)
  expect_equal(f1$Tht, f2$Tht, tolerance = 1e-8)
  expect_equal(f1$Omega, f2$Omega, tolerance = 1e-8)
})
