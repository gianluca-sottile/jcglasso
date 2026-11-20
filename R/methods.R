# Methods for jcglasso objects.

#' @export
print.jcglasso <- function(x, digits = 3L, ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Model:", x$model, "\nClasses:", paste(x$Z$classes, collapse = ", "),
      "\nn =", paste(x$nobs, collapse = ", "), " p =", x$nresp, " q =", x$npred, "\n\n")
  tab <- data.frame(lambda = rep(x$lambda, times = x$nrho), rho = rep(x$rho, each = x$nlambda),
                    nu = rep(x$nu, each = x$nlambda))
  if (x$npred == 0L) tab$lambda <- tab$nu <- NULL
  if (x$npred > 0L) tab$dfB <- as.vector(apply(x$dfB, c(2L, 3L), sum))
  tab$dfTht <- as.vector(apply(x$dfTht, c(2L, 3L), sum))
  if (x$npred > 0L) tab$dfOmg <- as.vector(apply(x$dfOmg, c(2L, 3L), sum))
  tab$EM <- as.vector(x$nit[1L, , ])
  print(format(tab, digits = digits), row.names = FALSE)
  cat("\n(alpha1 = ", x$alpha1, if (x$npred > 0L) paste0(", alpha2 = ", x$alpha2, ", alpha3 = ", x$alpha3),
      "; df: numbers of non-zero off-diagonal estimates summed over the classes)\n", sep = "")
  invisible(x)
}

#' Summarizing jcglasso Fits
#'
#' Prints the goodness-of-fit measure of the fitted models and the selected
#' model.
#'
#' @param object an object of class `"jcglasso"`.
#' @param GoF a function (`AIC` or `BIC`) or an object of class `"GoF"`.
#' @param print.all logical: print all the models or only the selected one?
#' @param digits number of significant digits.
#' @param ... further arguments passed to `GoF`.
#' @return Invisibly, a list with the table of models and the selected model
#'   (an object of class `"jcglasso"`).
#' @seealso [jcglasso()], [select_jcglasso()]
#' @examples
#' set.seed(1)
#' sim <- rjcggm(n = 50, p = 6, q = 2, K = 2, ncompPrecision = 2)
#' fit <- jcglasso(data = sim$ZX, nrho = 5, rho.min.ratio = 0.2)
#' summary(fit, GoF = BIC)
#' @export
summary.jcglasso <- function(object, GoF = AIC, print.all = TRUE, digits = 3L, ...) {
  if (is.function(GoF)) GoF <- GoF(object, ...)
  if (!inherits(GoF, "GoF")) stop(sQuote("GoF"), " must be AIC, BIC or an object of class ", sQuote("GoF"),
                                  call. = FALSE)
  tab <- data.frame(lambda.id = rep(seq_len(object$nlambda), times = object$nrho),
                    rho.id = rep(seq_len(object$nrho), each = object$nlambda),
                    lambda = rep(object$lambda, times = object$nrho), rho = rep(object$rho, each = object$nlambda),
                    df = as.vector(GoF$df), value = as.vector(GoF$value_gof))
  names(tab)[6L] <- GoF$type
  if (object$npred == 0L) tab$lambda.id <- tab$lambda <- NULL
  best <- which.min(tab[[GoF$type]])
  tab$best <- ""; tab$best[best] <- "<-"
  cat("\nModel:", object$model, "\n\n")
  print(format(if (print.all) tab else tab[best, ], digits = digits), row.names = FALSE)
  sel <- select_jcglasso(object, GoF = GoF)
  invisible(list(table = tab, best = sel))
}

#' Extract the Estimates of a jcglasso Model
#'
#' @param object an object of class `"jcglasso"`.
#' @param type which estimates: `"all"`, `"B"` (regression coefficients,
#'   including the intercepts), `"Theta"`, `"Omega"`, `"Sigma"` (the covariance
#'   matrix of the responses given the predictors) or `"Xi"` (the precision matrix
#'   of the responses and predictors).
#' @param class.id,lambda.id,rho.id indices of the classes and of the tuning
#'   parameters (by default, all).
#' @param drop logical: drop the dimensions of length one?
#' @param ... not used.
#' @return An array (or a list of arrays for `type = "all"`).
#' @seealso [jcglasso()], [fitted.jcglasso()]
#' @examples
#' set.seed(1)
#' sim <- rjcggm(n = 50, p = 6, q = 2, K = 2, ncompPrecision = 2)
#' fit <- jcglasso(data = sim$Z, nrho = 3, nlambda = 3, rho.min.ratio = 0.3,
#'                 lambda.min.ratio = 0.3)
#' coef(fit, type = "B", lambda.id = 3, rho.id = 3)
#' @export
coef.jcglasso <- function(object, type = c("all", "B", "Theta", "Omega", "Sigma", "Xi"), class.id, lambda.id,
                          rho.id, drop = TRUE, ...) {
  type <- match.arg(type)
  K <- length(object$nobs)
  if (missing(class.id)) class.id <- seq_len(K)
  if (missing(lambda.id)) lambda.id <- seq_len(object$nlambda)
  if (missing(rho.id)) rho.id <- seq_len(object$nrho)
  q <- object$npred
  get <- function(a) if (is.null(a)) NULL else a[, , class.id, lambda.id, rho.id, drop = drop]
  if (type == "all") return(list(B = get(object$B), Theta = get(object$Tht), Omega = get(object$Omega)))
  if (type %in% c("B", "Theta", "Omega"))
    return(get(switch(type, B = object$B, Theta = object$Tht, Omega = object$Omega)))
  # Sigma and Xi
  dims <- if (type == "Sigma") dim(object$Tht)[1:2] else rep(object$nresp + q, 2L)
  out <- array(0, c(dims, length(class.id), length(lambda.id), length(rho.id)))
  for (a in seq_along(class.id)) for (h in seq_along(lambda.id)) for (o in seq_along(rho.id)) {
    k <- class.id[a]; Th <- object$Tht[, , k, lambda.id[h], rho.id[o]]
    out[, , a, h, o] <- if (type == "Sigma") solve(Th) else
      build_Xi(Th, if (q) object$B[-1L, , k, lambda.id[h], rho.id[o]], if (q) object$Omega[, , k, lambda.id[h],
                                                                                           rho.id[o]])
  }
  nm <- if (type == "Sigma") dimnames(object$Tht)[1:2] else rep(list(c(dimnames(object$Tht)[[1L]],
                                                                      if (q) dimnames(object$Omega)[[1L]])), 2L)
  dimnames(out) <- c(nm, list(object$Z$classes[class.id], dimnames(object$Tht)[[4L]][lambda.id],
                              dimnames(object$Tht)[[5L]][rho.id]))
  if (drop) out <- drop(out)
  out
}

#' Imputation, Fitted Values and Residuals of jcglasso Models
#'
#' `impute()` returns the data with the censored and missing values replaced by
#' their conditional expected values under a fitted model (the E-step of the EM
#' algorithm). `fitted()` returns the fitted values of the responses
#' (\eqn{\beta_{0k} + B_k^\top x}, with the imputed predictors) and `residuals()`
#' the differences between the imputed responses and the fitted values.
#'
#' @param object an object of class `"jcglasso"`.
#' @param lambda.id,rho.id indices of the model (required if the path has more
#'   than one model).
#' @param ... not used.
#' @return A list of `K` matrices (one for each class).
#' @seealso [jcglasso()], [coef.jcglasso()]
#' @examples
#' set.seed(1)
#' sim <- rjcggm(n = 50, p = 6, q = 2, K = 2, ncompPrecision = 2)
#' fit <- jcglasso(data = sim$Z, nrho = 3, nlambda = 3, rho.min.ratio = 0.3,
#'                 lambda.min.ratio = 0.3)
#' Zhat <- impute(fit, lambda.id = 3, rho.id = 3)
#' head(Zhat$class1)
#' head(residuals(fit, lambda.id = 3, rho.id = 3)$class1)
#' @export
impute <- function(object, ...) UseMethod("impute")

#' @rdname impute
#' @export
impute.default <- function(object, ...) {
  if (requireNamespace("cglasso", quietly = TRUE)) return(cglasso::impute(object, ...))
  stop("no impute method for objects of class ", sQuote(class(object)[1L]), call. = FALSE)
}

#' @rdname impute
#' @export
impute.jcglasso <- function(object, lambda.id, rho.id, ...) {
  ids <- path_ids(object, lambda.id, rho.id)
  Z <- object$Z; q <- object$npred; p <- object$nresp
  out <- lapply(seq_along(Z$classes), function(k) {
    Zk <- if (q) cbind(Z$Y[[k]], Z$X[[k]]) else Z$Y[[k]]
    Th <- object$Tht[, , k, ids[1L], ids[2L]]
    Xi <- build_Xi(Th, if (q) object$B[-1L, , k, ids[1L], ids[2L]], if (q) object$Omega[, , k, ids[1L], ids[2L]])
    zlo <- Z$zm[, k] - 3 * sqrt(Z$zv[, k]); zup <- Z$zm[, k] + 3 * sqrt(Z$zv[, k])
    es <- k_estep(Zk, Z$R[[k]], object$psi[, k, ids[1L], ids[2L]], Xi, Z$lo, Z$up, zlo, zup)
    dimnames(es$Z) <- dimnames(Zk)
    es$Z
  })
  names(out) <- Z$classes
  out
}

path_ids <- function(object, lambda.id, rho.id) {
  if (missing(lambda.id)) { if (object$nlambda == 1L) lambda.id <- 1L else stop(sQuote("lambda.id"), " is missing",
                                                                              call. = FALSE) }
  if (missing(rho.id)) { if (object$nrho == 1L) rho.id <- 1L else stop(sQuote("rho.id"), " is missing",
                                                                     call. = FALSE) }
  check_index(lambda.id, object$nlambda, "lambda.id"); check_index(rho.id, object$nrho, "rho.id")
  c(lambda.id, rho.id)
}

#' @rdname impute
#' @export
fitted.jcglasso <- function(object, lambda.id, rho.id, ...) {
  ids <- path_ids(object, lambda.id, rho.id)
  Zi <- impute(object, lambda.id = ids[1L], rho.id = ids[2L])
  p <- object$nresp; q <- object$npred
  out <- lapply(seq_along(Zi), function(k) {
    B <- object$B[, , k, ids[1L], ids[2L]]
    B <- matrix(B, nrow = q + 1L, dimnames = dimnames(object$B)[1:2])
    X1 <- cbind(1, Zi[[k]][, p + seq_len(q), drop = FALSE])
    out <- X1 %*% B
    dimnames(out) <- list(rownames(Zi[[k]]), colnames(B))
    out
  })
  names(out) <- names(Zi)
  out
}

#' @rdname impute
#' @export
residuals.jcglasso <- function(object, lambda.id, rho.id, ...) {
  ids <- path_ids(object, lambda.id, rho.id)
  Zi <- impute(object, lambda.id = ids[1L], rho.id = ids[2L])
  Fi <- fitted(object, lambda.id = ids[1L], rho.id = ids[2L])
  p <- object$nresp
  out <- lapply(seq_along(Zi), function(k) Zi[[k]][, seq_len(p), drop = FALSE] - Fi[[k]])
  names(out) <- names(Zi)
  out
}

#' Predicted Values from a jcglasso Model
#'
#' @param object an object of class `"jcglasso"`.
#' @param newdata a list of `K` matrices (or data frames) with the values of the
#'   predictors for each class (complete data).
#' @param lambda.id,rho.id indices of the model.
#' @param ... not used.
#' @return A list of `K` matrices of predicted responses.
#' @seealso [fitted.jcglasso()]
#' @examples
#' set.seed(1)
#' sim <- rjcggm(n = 50, p = 6, q = 2, K = 2, ncompPrecision = 2)
#' fit <- jcglasso(data = sim$Z, nrho = 3, nlambda = 3, rho.min.ratio = 0.3,
#'                 lambda.min.ratio = 0.3)
#' nd <- list(matrix(0, 2, 2), matrix(1, 2, 2))
#' predict(fit, newdata = nd, lambda.id = 3, rho.id = 3)
#' @export
predict.jcglasso <- function(object, newdata, lambda.id, rho.id, ...) {
  ids <- path_ids(object, lambda.id, rho.id)
  q <- object$npred; K <- length(object$nobs)
  if (q == 0L) stop("the model has no predictors", call. = FALSE)
  if (!is.list(newdata) || is.data.frame(newdata) || length(newdata) != K)
    stop(sQuote("newdata"), " must be a list of ", K, " matrices", call. = FALSE)
  out <- lapply(seq_len(K), function(k) {
    X <- as.matrix(newdata[[k]])
    if (ncol(X) != q) stop("the matrices in ", sQuote("newdata"), " must have ", q, " columns", call. = FALSE)
    if (anyNA(X)) stop(sQuote("newdata"), " cannot contain missing values", call. = FALSE)
    cbind(1, X) %*% object$B[, , k, ids[1L], ids[2L]]
  })
  names(out) <- object$Z$classes
  out
}

#' Plot of the Coefficient Paths and of the Estimated Matrices
#'
#' With more than one model, plots the paths of the estimates against the tuning
#' parameter; with a single model (e.g., the output of [select_jcglasso()] or
#' [jcggm()]), plots the estimated matrices as images.
#'
#' @param x an object of class `"jcglasso"`.
#' @param what which estimates: `"Theta"`, `"B"` or `"Omega"`.
#' @param class.id the classes to plot (default all).
#' @param lambda.id,rho.id for paths, the index of the fixed tuning parameter
#'   (default 1).
#' @param ... further graphical parameters.
#' @return `x`, invisibly.
#' @examples
#' set.seed(1)
#' sim <- rjcggm(n = 50, p = 6, q = 2, K = 2, ncompPrecision = 2)
#' fit <- jcglasso(data = sim$ZX, nrho = 10, rho.min.ratio = 0.1)
#' plot(fit)
#' plot(select_jcglasso(fit, GoF = BIC))
#' @export
plot.jcglasso <- function(x, what = c("Theta", "B", "Omega"), class.id, lambda.id = 1L, rho.id = 1L, ...) {
  what <- match.arg(what)
  if (what != "Theta" && x$npred == 0L) stop("the model has no predictors", call. = FALSE)
  K <- length(x$nobs)
  if (missing(class.id)) class.id <- seq_len(K)
  A <- switch(what, Theta = x$Tht, B = x$B[-1L, , , , , drop = FALSE], Omega = x$Omega)
  op <- graphics::par(mfrow = grDevices::n2mfrow(length(class.id)))
  on.exit(graphics::par(op))
  single <- x$nlambda == 1L && x$nrho == 1L
  for (k in class.id) {
    if (single) {
      M <- A[, , k, 1L, 1L]
      if (what != "B") diag(M) <- 0
      M <- t(M[nrow(M):1L, , drop = FALSE])
      lim <- max(abs(M), 1e-10)
      graphics::image(seq_len(nrow(M)), seq_len(ncol(M)), M, zlim = c(-lim, lim), axes = FALSE, xlab = "", ylab = "",
                      col = grDevices::hcl.colors(51L, "Blue-Red 3"), main = paste(what, "-", x$Z$classes[k]), ...)
    } else {
      byrho <- x$nrho > 1L
      M <- if (byrho) A[, , k, lambda.id, ] else A[, , k, , rho.id]
      M <- matrix(M, ncol = if (byrho) x$nrho else x$nlambda)
      if (what != "B") M <- M[as.vector(upper.tri(A[, , k, 1L, 1L])), , drop = FALSE]
      tp <- if (byrho) x$rho else x$lambda
      graphics::matplot(tp, t(M), type = "l", lty = 1L, log = "x", xlab = if (byrho) expression(rho) else
                        expression(lambda), ylab = "estimates", main = paste(what, "-", x$Z$classes[k]), ...)
      graphics::abline(h = 0, lty = 2L)
    }
  }
  invisible(x)
}
