#' Simulate Data from Joint Conditional Gaussian Graphical Models
#'
#' Generates a dataset with `K` classes from the model of Sottile et al. (2025),
#' with sparse precision matrices sharing part of their structure across
#' classes, right-censored responses and missing-at-random values, as in the
#' simulation study of the paper.
#'
#' @param n sample size of each class.
#' @param p,q numbers of responses and predictors.
#' @param K number of classes.
#' @param ncompPrecision size of the star-shaped blocks of the precision matrices
#'   \eqn{\Theta_k} and \eqn{\Omega_k} (it should divide `p` and `q`).
#' @param ncompB number of non-zero regression coefficients for each response.
#' @param percConnectedCompK proportion of the edges of the first class that are
#'   removed in the other classes (a number, or a vector of length `K - 1`).
#' @param perc.NA proportion of missing values in the columns with missing values.
#' @param perc.X.na,perc.Y.na proportions of predictors and responses with missing
#'   values.
#' @param up right censoring value of the responses.
#' @param perc.cens expected proportion of censored values in the censored
#'   responses.
#' @param perc.Y.cens proportion of censored responses.
#' @param Thtmin,Thtmax,Omgmin,Omgmax,Bmin,Bmax ranges of the non-zero
#'   off-diagonal entries of \eqn{\Theta_k}, \eqn{\Omega_k} and of the non-zero
#'   regression coefficients.
#'
#' @return A list with components `ZX` (a `"datajcggm"` object with the predictors
#'   as responses, to be used to select \eqn{\nu}), `Z` (the `"datajcggm"` object
#'   with responses and predictors), and `true`, a list with the true parameters
#'   `B`, `Theta`, `Omega` and the complete data `X` and `Y`.
#' @seealso [datajcggm()], [jcglasso()]
#' @examples
#' set.seed(1)
#' sim <- rjcggm(n = 50, p = 10, q = 4, K = 2, ncompPrecision = 2)
#' summary(sim$Z)
#' @export
rjcggm <- function(n = 300L, p = 25L, q = 10L, K = 2L, ncompPrecision = 5L, ncompB = 2L, percConnectedCompK = 0.5,
                   perc.NA = 0.25, perc.X.na = 0.25, perc.Y.na = 0.25, up = 40, perc.cens = 0.25, perc.Y.cens = 0.25,
                   Thtmin = 0.3, Thtmax = 0.5, Omgmin = 0.3, Omgmax = 0.5, Bmin = 0.5, Bmax = 0.7) {
  for (a in c("n", "p", "q", "K", "ncompPrecision", "ncompB")) check_int(get(a), a)
  if (ncompB > q) stop(sQuote("ncompB"), " cannot be larger than ", sQuote("q"), call. = FALSE)
  if (length(percConnectedCompK) == 1L) percConnectedCompK <- rep(percConnectedCompK, K - 1L)
  if (length(percConnectedCompK) != K - 1L) stop(sQuote("percConnectedCompK"), " must have length 1 or K - 1",
                                                 call. = FALSE)
  unconnected <- c(NA, percConnectedCompK)
  perc.cens <- max(perc.cens, 1e-6)
  indexB <- vapply(seq_len(p), function(j) sort(sample(q, ncompB)), integer(ncompB))
  indexB <- matrix(indexB, nrow = ncompB)
  Om <- Th <- Bl <- X <- Y <- Xc <- Yc <- vector("list", K)
  cens <- rep(1e-6, p); cens[sample(p, floor(p * perc.Y.cens))] <- perc.cens
  for (k in seq_len(K)) {
    Om[[k]] <- sim_precision(q, ncompPrecision, unconnected[k], Omgmin, Omgmax)
    X[[k]] <- rmvn(n, solve(Om[[k]]))
    Bl[[k]] <- sim_B(X[[k]], cens, up, q, p, Bmin, Bmax, indexB, unconnected[k])
    Th[[k]] <- sim_precision(p, ncompPrecision, unconnected[k], Thtmin, Thtmax)
    Y[[k]] <- cbind(1, X[[k]]) %*% Bl[[k]] + rmvn(n, solve(Th[[k]]))
    colnames(Y[[k]]) <- paste0("Y", seq_len(p)); colnames(X[[k]]) <- paste0("X", seq_len(q))
    Xc[[k]] <- X[[k]]; Yc[[k]] <- Y[[k]]
    nx <- floor(q * perc.X.na); ny <- floor(p * perc.Y.na); nr <- floor(n * perc.NA)
    if (nx && nr) for (j in sample(q, nx)) X[[k]][sample(n, nr), j] <- NA
    if (ny && nr) for (j in sample(p, ny)) Y[[k]][sample(n, nr), j] <- NA
  }
  names(X) <- names(Y) <- paste0("class", seq_len(K))
  list(ZX = suppressMessages(datajcggm(Y = X, up = up)), Z = suppressMessages(datajcggm(Y = Y, X = X, up = up)),
       true = list(B = Bl, Theta = Th, Omega = Om, X = Xc, Y = Yc))
}

rmvn <- function(n, Sigma) {
  matrix(stats::rnorm(n * ncol(Sigma)), n) %*% chol(Sigma)
}

# Block diagonal precision matrix with star-shaped blocks; a proportion 'unconnected'
# of its edges is removed.
sim_precision <- function(p, size, unconnected = NA, tmin = 0.3, tmax = 0.5) {
  size <- min(size, p)
  nb <- ceiling(p / size)
  Tht <- diag(p)
  for (b in seq_len(nb)) {
    id <- ((b - 1L) * size + 1L):min(b * size, p)
    if (length(id) > 1L) {
      v <- sort(stats::runif(length(id) - 1L, tmin, tmax), decreasing = TRUE)
      Tht[id[1L], id[-1L]] <- Tht[id[-1L], id[1L]] <- v
    }
  }
  if (!is.na(unconnected)) {
    nz <- which(Tht[upper.tri(Tht)] != 0)
    off <- nz[sample.int(length(nz), floor(unconnected * length(nz)))]
    U <- Tht[upper.tri(Tht)]; U[off] <- 0
    Tht[upper.tri(Tht)] <- U
    Tht[lower.tri(Tht)] <- t(Tht)[lower.tri(Tht)]
  }
  # ensure positive definiteness
  ev <- min(eigen(Tht, symmetric = TRUE, only.values = TRUE)$values)
  if (ev < 0.1) diag(Tht) <- diag(Tht) + 0.1 - ev
  Tht
}

# Regression coefficients: the intercepts give the requested probability of censoring.
sim_B <- function(X, cens, up, q, p, Bmin, Bmax, indexB, unconnected = NA) {
  B <- matrix(0, q + 1L, p)
  for (j in seq_len(p)) B[1L + indexB[, j], j] <- stats::runif(nrow(indexB), Bmin, Bmax)
  if (!is.na(unconnected)) {
    nz <- which(B[-1L, ] != 0)
    off <- nz[sample.int(length(nz), floor(unconnected * length(nz)))]
    Bx <- B[-1L, , drop = FALSE]; Bx[off] <- 0; B[-1L, ] <- Bx
  }
  eta <- X %*% B[-1L, , drop = FALSE]
  if (is.finite(up)) {
    for (j in seq_len(p)) {
      fn <- function(b0) mean(stats::pnorm(up - b0 - eta[, j], lower.tail = FALSE)) - cens[j]
      B[1L, j] <- stats::uniroot(fn, c(up - 100, up + 100), tol = 1e-8)$root
    }
  }
  rownames(B) <- c("(Intercept)", paste0("X", seq_len(q)))
  colnames(B) <- paste0("Y", seq_len(p))
  B
}
