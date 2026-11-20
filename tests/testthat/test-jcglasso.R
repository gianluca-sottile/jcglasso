test_that("path dimensions, maximal tuning parameters and optimality", {
  sim <- small_sim()
  fit <- jcglasso(data = sim$Z, nrho = 3, nlambda = 3, rho.min.ratio = 0.3, lambda.min.ratio = 0.3, nu = 0.1)
  expect_s3_class(fit, "jcglasso")
  expect_equal(dim(fit$B), c(4L, 6L, 2L, 3L, 3L))
  expect_equal(dim(fit$Tht), c(6L, 6L, 2L, 3L, 3L))
  expect_equal(dim(fit$Omega), c(3L, 3L, 2L, 3L, 3L))
  expect_true(all(fit$conv == 0L))
  expect_true(all(fit$dfB[, 1L, ] == 0))
  expect_true(all(fit$dfTht[, , 1L][, 1L] == 0))
  v <- kkt(fit, lambda.id = 3, rho.id = 3)
  expect_lt(max(v), 1e-3)
  expect_output(print(fit), "conditional")
})

test_that("model without predictors and K = 1", {
  sim <- small_sim(K = 1)
  fit <- jcglasso(data = sim$ZX, nrho = 4, rho.min.ratio = 0.2)
  expect_equal(fit$npred, 0L)
  expect_null(fit$Omega)
  expect_lt(max(kkt(fit, rho.id = 4), na.rm = TRUE), 1e-3)
})

test_that("formula selects responses and predictors", {
  sim <- small_sim()
  fit <- jcglasso(cbind(Y1, Y2, Y3) ~ X1 + X2, data = sim$Z, nrho = 2, nlambda = 2)
  expect_equal(fit$nresp, 3L); expect_equal(fit$npred, 2L)
  fit0 <- jcglasso(. ~ 1, data = sim$Z, nrho = 2)
  expect_equal(fit0$npred, 0L)
  expect_error(jcglasso(. ~ log(X1), data = sim$Z), "predictors")
  expect_error(jcglasso(cbind(Y1, Z9) ~ ., data = sim$Z), "not a response")
})

test_that("weights define structural zeros and unpenalized parameters", {
  sim <- small_sim()
  wB <- array(1, c(3, 6, 2)); wB[1, 1, ] <- Inf; wB[2, 2, ] <- 0
  fit <- jcglasso(data = sim$Z, weights.B = wB, nrho = 1, nlambda = 2, lambda.min.ratio = 0.5, nu = 0.1)
  expect_true(all(fit$B[2, 1, , , ] == 0))
  expect_true(all(fit$B[3, 2, , 2, ] != 0))
  expect_error(jcglasso(data = sim$Z, weights.B = wB[, , 1]), "dimensions")
})

test_that("argument checks", {
  sim <- small_sim()
  expect_error(jcglasso(data = sim$Z, alpha1 = 2), "alpha1")
  expect_error(jcglasso(data = sim$Z$Y), "datajcggm")
  expect_error(jcglasso(data = sim$ZX, lambda = 1), "predictors")
})
