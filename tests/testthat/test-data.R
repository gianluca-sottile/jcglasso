test_that("datajcggm codes censored and missing values", {
  set.seed(1)
  Y <- matrix(rnorm(40 * 3), 40, 3)
  Y[Y > 1] <- 1; Y[Y < -1.5] <- -1.5
  Y[1:3, 2] <- NA
  X <- matrix(rnorm(40 * 2), 40, 2)
  Z <- datajcggm(Y = list(a = Y, b = Y), lo = c(rep(-1.5, 3), -Inf, -Inf), up = c(rep(1, 3), Inf, Inf), X = list(X, X))
  expect_true(is.datajcggm(Z))
  expect_equal(nobs(Z), c(a = 40L, b = 40L))
  expect_equal(nresp(Z), 3L); expect_equal(npred(Z), 2L)
  R <- event(Z)$a
  expect_equal(sum(R == 9L), 3L)
  expect_equal(sum(R[, 1:3] == 1L), sum(Y >= 1, na.rm = TRUE))
  expect_equal(sum(R[, 1:3] == -1L), sum(Y <= -1.5, na.rm = TRUE))
  expect_output(print(Z), "datajcggm")
  expect_output(print(summary(Z)), "Responses")
  expect_equal(dim(ColMeans(Z)), c(5L, 2L))
  expect_error(datajcggm(Y = list(Y, Y[, 1:2])), "same columns")
  expect_error(datajcggm(Y = Y, lo = 2, up = 1), "smaller")
})

test_that("marginal estimates are the censored-normal MLE", {
  set.seed(2)
  y <- rnorm(300, 1, 2); up <- 2; y[y > up] <- up
  Z <- datajcggm(Y = cbind(y, rnorm(300)), up = c(up, Inf))
  nll <- function(t) { m <- t[1]; s <- exp(t[2])
    -sum(dnorm(y[y < up], m, s, log = TRUE)) - sum(y >= up) * pnorm(up, m, s, lower.tail = FALSE, log.p = TRUE) }
  o <- optim(c(0, 0), nll, control = list(reltol = 1e-12))
  expect_equal(unname(ColMeans(Z)[1, 1]), o$par[1], tolerance = 1e-3)
  expect_equal(unname(ColVars(Z)[1, 1]), exp(2 * o$par[2]), tolerance = 1e-3)
})

test_that("rjcggm returns datasets and true parameters", {
  sim <- small_sim()
  expect_true(is.datajcggm(sim$Z)); expect_true(is.datajcggm(sim$ZX))
  expect_equal(npred(sim$Z), 3L); expect_equal(npred(sim$ZX), 0L)
  expect_length(sim$true$Theta, 2L)
  expect_true(all(vapply(sim$true$Theta, function(M) min(eigen(M)$values) > 0, TRUE)))
})
