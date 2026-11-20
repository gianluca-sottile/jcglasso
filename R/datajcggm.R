#' Datasets for Joint Conditional Gaussian Graphical Models
#'
#' Creates an object of class `"datajcggm"`, that is, a dataset with `K`
#' classes (sub-populations, or experimental conditions) of observations drawn
#' from conditional Gaussian graphical models, with censored and/or missing
#' values. Each class contains a matrix of response variables `Y` and,
#' optionally, a matrix of predictors `X` (with the same variables in all the
#' classes).
#'
#' @param Y a list of `K` numeric matrices with the same columns (one for each
#'   class), or a single matrix (`K = 1`). Missing values are coded as `NA`.
#' @param lo,up the lower and upper censoring values: a number (common to all
#'   variables) or a vector with one value for the response variables followed by
#'   one value for the predictors. Values equal to `lo` (`up`) are treated as
#'   left (right) censored; values below `lo` (above `up`) are replaced by `lo`
#'   (`up`). Use `-Inf`/`Inf` for no censoring.
#' @param X an optional list of `K` numeric matrices of predictors (or a single
#'   matrix), with the same numbers of rows as `Y`. Predictors can also be
#'   censored or missing.
#' @param control a list with the maximum number of iterations (`maxit`) and the
#'   tolerance (`thr`) of the EM algorithm used to estimate the marginal means and
#'   variances of the variables.
#'
#' @details
#' For each class and each variable, the marginal mean and variance are estimated
#' by maximum likelihood under a univariate censored Gaussian model; they are used
#' as starting values by [jcglasso()] and are returned by [ColMeans()] and
#' [ColVars()].
#'
#' @return An object of class `"datajcggm"`, a list with components `Y` and `X`
#'   (lists of `K` matrices; `X` is `NULL` without predictors), `lo` and `up`
#'   (censoring values), `R` (a list of `K` matrices of codes: `0` observed, `-1`
#'   left censored, `1` right censored, `9` missing at random), `zm` and `zv`
#'   (marginal means and variances, matrices with one column for each class), and
#'   `classes` (class labels).
#' @seealso [jcglasso()], [rjcggm()], [getMatrix()], [event()]
#' @references Sottile, G., Augugliaro, L., Vinciotti, V., Arancio, W. and
#'   Coronnello, C. (2025). Sparse inference of the human haematopoietic system
#'   from heterogeneous and partially observed genomic data. *Journal of the Royal
#'   Statistical Society Series C: Applied Statistics*, 74(1), 204-228.
#'   \doi{10.1093/jrsssc/qlae052}
#' @examples
#' set.seed(1)
#' Y1 <- matrix(rnorm(50 * 3), 50, 3); Y2 <- matrix(rnorm(40 * 3), 40, 3)
#' Y1[Y1 > 1.5] <- 1.5; Y2[Y2 > 1.5] <- 1.5   # right censoring at 1.5
#' Y2[1:4, 1] <- NA                           # missing at random
#' Z <- datajcggm(Y = list(A = Y1, B = Y2), up = 1.5)
#' Z
#' summary(Z)
#' @export
datajcggm <- function(Y, lo = -Inf, up = Inf, X = NULL, control = list(maxit = 1e4, thr = 1e-4)) {
  if (missing(Y)) stop(sQuote("Y"), " is missing", call. = FALSE)
  if (is.matrix(Y) || is.data.frame(Y)) Y <- list(Y)
  if (!is.list(Y) || !length(Y)) stop(sQuote("Y"), " must be a matrix or a list of matrices", call. = FALSE)
  K <- length(Y)
  classes <- names(Y)
  if (is.null(classes) || any(classes == "")) classes <- paste0("class", seq_len(K))
  Y <- lapply(Y, as_num_matrix, what = "Y")
  p <- ncol(Y[[1L]])
  ynames <- colnames(Y[[1L]]); if (is.null(ynames)) ynames <- paste0("Y", seq_len(p))
  for (k in seq_len(K)) {
    if (ncol(Y[[k]]) != p) stop("all the matrices in ", sQuote("Y"), " must have the same columns", call. = FALSE)
    if (!is.null(colnames(Y[[k]])) && !identical(colnames(Y[[k]]), ynames))
      stop("all the matrices in ", sQuote("Y"), " must have the same column names", call. = FALSE)
    colnames(Y[[k]]) <- ynames
    if (is.null(rownames(Y[[k]]))) rownames(Y[[k]]) <- seq_len(nrow(Y[[k]]))
  }
  q <- 0L; xnames <- NULL
  if (!is.null(X)) {
    if (is.matrix(X) || is.data.frame(X)) X <- list(X)
    if (length(X) != K) stop(sQuote("X"), " and ", sQuote("Y"), " have a different number of classes", call. = FALSE)
    X <- lapply(X, as_num_matrix, what = "X")
    q <- ncol(X[[1L]])
    xnames <- colnames(X[[1L]]); if (is.null(xnames)) xnames <- paste0("X", seq_len(q))
    for (k in seq_len(K)) {
      if (ncol(X[[k]]) != q) stop("all the matrices in ", sQuote("X"), " must have the same columns", call. = FALSE)
      if (nrow(X[[k]]) != nrow(Y[[k]])) stop(sQuote("X"), " and ", sQuote("Y"), " have a different number of rows",
                                             call. = FALSE)
      colnames(X[[k]]) <- xnames
      rownames(X[[k]]) <- rownames(Y[[k]])
    }
    if (any(xnames %in% ynames)) stop("predictors and responses must have different names", call. = FALSE)
  }
  d <- p + q
  znames <- c(ynames, xnames)
  lo <- expand_limits(lo, d, znames, "lo")
  up <- expand_limits(up, d, znames, "up")
  if (any(lo >= up)) stop(sQuote("lo"), " must be smaller than ", sQuote("up"), call. = FALSE)
  if (!is.list(control)) stop(sQuote("control"), " must be a list", call. = FALSE)
  maxit <- if (is.null(control$maxit)) 1e4 else control$maxit
  thr <- if (is.null(control$thr)) 1e-4 else control$thr

  R <- vector("list", K)
  zm <- zv <- matrix(0, d, K, dimnames = list(znames, classes))
  for (k in seq_len(K)) {
    Zk <- if (q) cbind(Y[[k]], X[[k]]) else Y[[k]]
    for (m in seq_len(d)) {
      below <- which(Zk[, m] < lo[m]); above <- which(Zk[, m] > up[m])
      if (length(below)) {
        message("some entries of ", sQuote(znames[m]), " (class ", classes[k],
                ") are below the lower censoring value and are replaced by it")
        Zk[below, m] <- lo[m]
      }
      if (length(above)) {
        message("some entries of ", sQuote(znames[m]), " (class ", classes[k],
                ") are above the upper censoring value and are replaced by it")
        Zk[above, m] <- up[m]
      }
    }
    Rk <- matrix(0L, nrow(Zk), d, dimnames = dimnames(Zk))
    Rk[!is.na(Zk) & Zk <= rep(lo, each = nrow(Zk))] <- -1L
    Rk[!is.na(Zk) & Zk >= rep(up, each = nrow(Zk))] <- 1L
    Rk[is.na(Zk)] <- 9L
    if (any(colSums(Rk == 0L) < 2L))
      stop("each variable must have at least two observed values in each class", call. = FALSE)
    R[[k]] <- Rk
    Y[[k]] <- Zk[, seq_len(p), drop = FALSE]
    if (q) X[[k]] <- Zk[, p + seq_len(q), drop = FALSE]
    mv <- marginal_censored(Zk, Rk, lo, up, maxit, thr)
    zm[, k] <- mv$m; zv[, k] <- mv$v
  }
  names(Y) <- names(R) <- classes
  if (q) names(X) <- classes
  structure(list(Y = Y, X = if (q) X else NULL, lo = lo, up = up, R = R, zm = zm, zv = zv, classes = classes),
            class = "datajcggm")
}

as_num_matrix <- function(x, what) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.vector(x) && is.numeric(x)) x <- matrix(x, ncol = 1L)
  if (!is.matrix(x) || !is.numeric(x)) stop(sQuote(what), " must contain numeric matrices", call. = FALSE)
  if (any(is.infinite(x))) stop(sQuote(what), " contains infinite values", call. = FALSE)
  storage.mode(x) <- "double"
  x
}

expand_limits <- function(x, d, nms, what) {
  if (!is.numeric(x) || anyNA(x)) stop(sQuote(what), " must be a numeric vector", call. = FALSE)
  if (length(x) == 1L) x <- rep(x, d)
  if (length(x) != d) stop(sQuote(what), " must have length 1 or ", d, call. = FALSE)
  names(x) <- nms
  x
}

# Univariate censored Gaussian MLE (EM) for each column: missing-at-random entries are
# ignored, censored entries contribute through the truncated normal moments.
marginal_censored <- function(Z, R, lo, up, maxit = 1e4, thr = 1e-4) {
  d <- ncol(Z)
  m <- v <- numeric(d)
  for (j in seq_len(d)) {
    r <- R[, j]; z <- Z[r != 9L, j]; r <- r[r != 9L]
    obs <- z[r == 0L]
    mj <- mean(obs); vj <- max(stats::var(obs), 1e-8 * max(1, mj^2))
    if (any(r != 0L)) {
      nl <- sum(r == -1L); nr <- sum(r == 1L); n <- length(r)
      for (it in seq_len(maxit)) {
        s <- sqrt(vj)
        e1 <- e2 <- 0
        if (nl) {
          a <- (lo[j] - mj) / s
          ratio <- exp(stats::dnorm(a, log = TRUE) - stats::pnorm(a, log.p = TRUE))
          tm <- mj - s * ratio; tv <- vj * (1 - a * ratio - ratio^2)
          e1 <- e1 + nl * tm; e2 <- e2 + nl * (max(tv, 0) + tm^2)
        }
        if (nr) {
          b <- (up[j] - mj) / s
          ratio <- exp(stats::dnorm(b, log = TRUE) - stats::pnorm(b, lower.tail = FALSE, log.p = TRUE))
          tm <- mj + s * ratio; tv <- vj * (1 + b * ratio - ratio^2)
          e1 <- e1 + nr * tm; e2 <- e2 + nr * (max(tv, 0) + tm^2)
        }
        mn <- (sum(obs) + e1) / n
        vn <- max((sum(obs^2) + e2) / n - mn^2, 1e-8 * max(1, mn^2))
        dm <- abs(mn - mj) / max(1, abs(mj)); dv <- abs(vn - vj) / max(1, vj)
        mj <- mn; vj <- vn
        if (max(dm, dv) < thr) break
      }
    }
    m[j] <- mj; v[j] <- vj
  }
  list(m = m, v = v)
}

#' @rdname datajcggm
#' @param x an object of class `"datajcggm"`.
#' @export
is.datajcggm <- function(x) inherits(x, "datajcggm")

#' @export
nobs.datajcggm <- function(object, ...) vapply(object$Y, nrow, 1L)

#' @export
print.datajcggm <- function(x, n = 6L, digits = 3L, ...) {
  nk <- nobs(x)
  cat("Dataset of class", sQuote("datajcggm"), "with", length(nk), "classes,",
      nresp(x), "responses and", npred(x), "predictors\n")
  for (k in seq_along(nk)) {
    cat("\nClass ", sQuote(x$classes[k]), " (n = ", nk[k], ")\n", sep = "")
    Zc <- format_censored(x, k, digits)
    print(utils::head(Zc, n), quote = FALSE, right = TRUE)
    if (nk[k] > n) cat("# ...", nk[k] - n, "more rows\n")
  }
  cat("\n'+' right censored, '-' left censored, NA missing\n")
  invisible(x)
}

format_censored <- function(x, k, digits) {
  Z <- cbind(x$Y[[k]], if (!is.null(x$X)) x$X[[k]])
  out <- format(Z, digits = digits)
  out[x$R[[k]] == 1L] <- paste0(out[x$R[[k]] == 1L], "+")
  out[x$R[[k]] == -1L] <- paste0(out[x$R[[k]] == -1L], "-")
  out
}

#' @rdname datajcggm
#' @param object an object of class `"datajcggm"`.
#' @param ... not used.
#' @export
summary.datajcggm <- function(object, ...) {
  out <- lapply(seq_along(object$classes), function(k) {
    Z <- cbind(object$Y[[k]], if (!is.null(object$X)) object$X[[k]])
    R <- object$R[[k]]
    tab <- t(vapply(seq_len(ncol(Z)), function(j) {
      o <- Z[R[, j] == 0L, j]
      c(Lower = unname(object$lo[j]), Min = min(o), Median = stats::median(o), Mean = mean(o), Max = max(o),
        Upper = unname(object$up[j]), `NA%` = 100 * mean(R[, j] == 9L), `LC%` = 100 * mean(R[, j] == -1L),
        `RC%` = 100 * mean(R[, j] == 1L))
    }, numeric(9L)))
    rownames(tab) <- colnames(Z)
    tab
  })
  names(out) <- object$classes
  structure(out, class = "summary.datajcggm", p = ncol(object$Y[[1L]]))
}

#' @export
print.summary.datajcggm <- function(x, digits = 3L, ...) {
  p <- attr(x, "p")
  for (k in names(x)) {
    cat("Class", sQuote(k), "\n")
    tab <- x[[k]]
    cat(" Responses:\n"); print(round(tab[seq_len(p), , drop = FALSE], digits))
    if (nrow(tab) > p) { cat(" Predictors:\n"); print(round(tab[-seq_len(p), , drop = FALSE], digits)) }
    cat("\n")
  }
  invisible(x)
}
