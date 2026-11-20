#' Q-Function of a Fitted jcglasso Model
#'
#' Computes the value of the Q-function (the expected complete-data
#' log-likelihood used in the M-step of the EM algorithm) of the models fitted
#' along the path, and the corresponding degrees of freedom. The Q-function is
#' used to define the goodness-of-fit measures [AIC.jcglasso()] and
#' [BIC.jcglasso()] (Sottile et al., 2025, Section 6).
#'
#' @param object an object of class `"jcglasso"` or `"jcggm"`.
#' @param mle logical: if `TRUE`, the Q-function is evaluated at the maximum
#'   likelihood estimates constrained to the sparsity pattern of each fitted
#'   model (computed by [jcggm()]); otherwise at the penalized estimates.
#' @param part `"conditional"` (default) for the model of the responses given the
#'   predictors (\eqn{Q_{y|x}}), `"joint"` for the sum of the conditional model and
#'   of the model of the predictors (\eqn{Q_{y|x} + Q_x}). Without predictors, the
#'   two are identical.
#' @param verbose logical: show a progress bar when `mle = TRUE`.
#' @param ... further arguments passed to [jcggm()].
#'
#' @details The degrees of freedom of each class are the numbers of non-zero
#'   regression coefficients and off-diagonal entries of \eqn{\Theta_k}, plus the
#'   `p` intercepts and the `p` diagonal entries (and, for `part = "joint"`, the
#'   non-zero off-diagonal entries of \eqn{\Omega_k} plus `2q`).
#' @return An object of class `"QFun"`, a list with the arrays `value` and `df`
#'   (`K x nlambda x nrho`), the sample sizes `n`, the tuning parameters and the
#'   argument `part`.
#' @seealso [AIC.jcglasso()], [BIC.jcglasso()], [jcglasso()]
#' @examples
#' set.seed(1)
#' sim <- rjcggm(n = 50, p = 6, q = 2, K = 2, ncompPrecision = 2)
#' fit <- jcglasso(data = sim$ZX, nrho = 5, rho.min.ratio = 0.2)
#' QFun(fit)
#' @export
QFun <- function(object, ...) UseMethod("QFun")

#' @rdname QFun
#' @export
QFun.default <- function(object, ...) {
  if (requireNamespace("cglasso", quietly = TRUE)) return(cglasso::QFun(object, ...))
  stop("no QFun method for objects of class ", sQuote(class(object)[1L]), call. = FALSE)
}

#' @rdname QFun
#' @export
QFun.jcglasso <- function(object, mle = FALSE, part = c("conditional", "joint"), verbose = FALSE, ...) {
  part <- match.arg(part)
  check_logical(mle, "mle")
  if (inherits(object, "jcggm")) mle <- FALSE
  q <- object$npred; p <- object$nresp
  joint <- part == "joint" && q > 0L
  if (!mle) {
    value <- object$Q$yx + if (joint) object$Q$x else 0
    df <- dfs(object, joint)
  } else {
    nl <- object$nlambda; nr <- object$nrho; K <- length(object$nobs)
    value <- df <- array(0, c(K, nl, nr), dimnames(object$Q$yx))
    seen <- list()
    if (verbose) pb <- utils::txtProgressBar(0, nl * nr, style = 3L)
    for (h in seq_len(nl)) for (o in seq_len(nr)) {
      key <- support_key(object, h, o, joint)
      if (!is.null(seen[[key]])) { value[, h, o] <- seen[[key]]$v; df[, h, o] <- seen[[key]]$d }
      else {
        m <- jcggm(object, lambda.id = h, rho.id = o, ...)
        v <- m$Q$yx[, 1L, 1L] + if (joint) m$Q$x[, 1L, 1L] else 0
        dd <- dfs(m, joint)[, 1L, 1L]
        value[, h, o] <- v; df[, h, o] <- dd
        seen[[key]] <- list(v = v, d = dd)
      }
      if (verbose) utils::setTxtProgressBar(pb, (h - 1L) * nr + o)
    }
    if (verbose) close(pb)
  }
  structure(list(value = value, df = df, n = object$nobs, p = p, q = q, lambda = object$lambda, rho = object$rho,
                 nu = object$nu, nlambda = object$nlambda, nrho = object$nrho, part = part, mle = mle,
                 model = object$model),
            class = "QFun")
}

dfs <- function(object, joint) {
  p <- object$nresp; q <- object$npred
  df <- object$dfTht + 2 * p
  if (q > 0L) df <- df + object$dfB
  if (joint) df <- df + object$dfOmg + 2 * q
  df
}

support_key <- function(object, h, o, joint) {
  s <- as.integer(object$Tht[, , , h, o, drop = FALSE] != 0)
  if (object$npred > 0L) s <- c(s, as.integer(object$B[-1L, , , h, o, drop = FALSE] != 0))
  if (joint) s <- c(s, as.integer(object$Omega[, , , h, o, drop = FALSE] != 0))
  paste(s, collapse = "")
}

#' @export
print.QFun <- function(x, digits = 3L, ...) {
  cat("Q-function (", x$part, " model", if (x$mle) ", maximum likelihood estimates", ")\n\n", sep = "")
  tab <- data.frame(lambda = rep(x$lambda, times = x$nrho), rho = rep(x$rho, each = x$nlambda),
                    df = as.vector(apply(x$df, c(2L, 3L), sum)), Q = as.vector(apply(x$value, c(2L, 3L), sum)))
  if (x$q == 0L) tab$lambda <- NULL
  print(format(tab, digits = digits), row.names = FALSE)
  invisible(x)
}

#' Information Criteria for jcglasso Models
#'
#' Computes the AIC, BIC and extended BIC of the models fitted along the path,
#' using the Q-function in place of the observed log-likelihood (Sottile et al.,
#' 2025, Section 6).
#'
#' @param object an object of class `"jcglasso"`.
#' @param k the penalty per degree of freedom (`2` for the AIC).
#' @param g the parameter of the extended BIC, in \eqn{[0, 1]}; `g = 0` gives the
#'   BIC.
#' @param type the extended BIC: `"FD"` (Foygel and Drton, 2010, penalty
#'   \eqn{\log n_k + 4 g \log p}) or `"CC"` (Chen and Chen, 2008, penalty
#'   \eqn{\log n_k + 2 g \log q}).
#' @param mle,part see [QFun()].
#' @param Qfun an optional object of class `"QFun"` (to avoid recomputing it).
#' @param ... further arguments passed to [QFun()].
#' @return An object of class `"GoF"`: a list with the matrix `value_gof`
#'   (`nlambda x nrho`) of the criterion, the degrees of freedom and the tuning
#'   parameters.
#' @seealso [select_jcglasso()], [QFun()]
#' @examples
#' set.seed(1)
#' sim <- rjcggm(n = 50, p = 6, q = 2, K = 2, ncompPrecision = 2)
#' fit <- jcglasso(data = sim$ZX, nrho = 5, rho.min.ratio = 0.2)
#' BIC(fit)
#' plot(BIC(fit))
#' @export
AIC.jcglasso <- function(object, k = 2, mle = FALSE, part = c("conditional", "joint"), Qfun = NULL, ...) {
  Qf <- if (is.null(Qfun)) QFun(object, mle = mle, part = match.arg(part), ...) else Qfun
  K <- length(Qf$n)
  if (length(k) == 1L) k <- rep(k, K)
  if (length(k) != K || any(k <= 0)) stop(sQuote("k"), " must be positive", call. = FALSE)
  gof_object(Qf, k, if (all(k == 2)) "AIC" else "GoF")
}

#' @rdname AIC.jcglasso
#' @export
BIC.jcglasso <- function(object, g = 0, type = c("FD", "CC"), mle = FALSE, part = c("conditional", "joint"),
                         Qfun = NULL, ...) {
  check_scalar(g, "g", 0, 1)
  type <- match.arg(type)
  Qf <- if (is.null(Qfun)) QFun(object, mle = mle, part = match.arg(part), ...) else Qfun
  k <- log(Qf$n)
  if (g > 0) {
    if (type == "CC" && Qf$q == 0L) stop("type = \"CC\" requires predictors", call. = FALSE)
    k <- k + if (type == "FD") 4 * g * log(Qf$p) else 2 * g * log(Qf$q)
  }
  gof_object(Qf, k, if (g == 0) "BIC" else paste0("eBIC_", type))
}

gof_object <- function(Qf, k, type) {
  K <- length(Qf$n)
  val <- -2 * Qf$value + array(k, dim(Qf$df)) * Qf$df
  value_gof <- apply(val, c(2L, 3L), sum)
  dim(value_gof) <- c(Qf$nlambda, Qf$nrho)
  dimnames(value_gof) <- dimnames(Qf$value)[2:3]
  structure(list(value_gof = value_gof, df = apply(Qf$df, c(2L, 3L), sum), value = Qf$value, n = Qf$n, p = Qf$p,
                 q = Qf$q, lambda = Qf$lambda, rho = Qf$rho, nlambda = Qf$nlambda, nrho = Qf$nrho, type = type,
                 part = Qf$part, mle = Qf$mle, model = Qf$model),
            class = "GoF")
}

#' @export
print.GoF <- function(x, digits = 3L, ...) {
  cat(x$type, " (", x$part, " model", if (x$mle) ", maximum likelihood estimates", ")\n\n", sep = "")
  tab <- data.frame(lambda = rep(x$lambda, times = x$nrho), rho = rep(x$rho, each = x$nlambda),
                    df = as.vector(x$df), value = as.vector(x$value_gof))
  if (x$q == 0L) tab$lambda <- NULL
  print(format(tab, digits = digits), row.names = FALSE)
  best <- which(x$value_gof == min(x$value_gof), arr.ind = TRUE)[1L, ]
  cat("\nminimum at lambda.id =", best[1L], "and rho.id =", best[2L], "\n")
  invisible(x)
}

#' @rdname AIC.jcglasso
#' @param x an object of class `"GoF"`.
#' @export
plot.GoF <- function(x, ...) {
  v <- x$value_gof
  if (x$nlambda == 1L || x$nrho == 1L) {
    tp <- if (x$nrho > 1L) x$rho else x$lambda
    graphics::plot(tp, as.vector(v), type = "b", log = "x", xlab = if (x$nrho > 1L) expression(rho) else
                   expression(lambda), ylab = x$type, main = x$model, ...)
    graphics::abline(v = tp[which.min(v)], lty = 2L, col = "red")
  } else {
    graphics::contour(log(x$lambda[order(x$lambda)]), log(x$rho[order(x$rho)]),
                      v[order(x$lambda), order(x$rho)], xlab = expression(log(lambda)), ylab = expression(log(rho)),
                      main = paste(x$type, "-", x$model), ...)
    best <- which(v == min(v), arr.ind = TRUE)[1L, ]
    graphics::points(log(x$lambda[best[1L]]), log(x$rho[best[2L]]), pch = 19L, col = "red")
  }
  invisible(x)
}

#' Select a Fitted Model from the Path
#'
#' Extracts from an object of class `"jcglasso"` the model that minimizes a
#' goodness-of-fit measure, or the model identified by `lambda.id` and `rho.id`.
#'
#' @param object an object of class `"jcglasso"`.
#' @param GoF a function (`AIC` or `BIC`) or an object of class `"GoF"`.
#' @param lambda.id,rho.id optional indices of the model.
#' @param ... further arguments passed to `GoF` when it is a function.
#' @return An object of class `"jcglasso"` with a single model.
#' @seealso [AIC.jcglasso()], [BIC.jcglasso()], [jcggm()]
#' @examples
#' set.seed(1)
#' sim <- rjcggm(n = 50, p = 6, q = 2, K = 2, ncompPrecision = 2)
#' fit <- jcglasso(data = sim$ZX, nrho = 5, rho.min.ratio = 0.2)
#' best <- select_jcglasso(fit, GoF = BIC)
#' best$rho
#' @export
select_jcglasso <- function(object, GoF = AIC, lambda.id, rho.id, ...) {
  if (!inherits(object, "jcglasso")) stop(sQuote("object"), " must be an object of class ", sQuote("jcglasso"),
                                          call. = FALSE)
  if (missing(lambda.id) && missing(rho.id)) {
    if (is.function(GoF)) GoF <- GoF(object, ...)
    if (!inherits(GoF, "GoF")) stop(sQuote("GoF"), " must be AIC, BIC or an object of class ", sQuote("GoF"),
                                    call. = FALSE)
    best <- which(GoF$value_gof == min(GoF$value_gof), arr.ind = TRUE)[1L, ]
    lambda.id <- best[1L]; rho.id <- best[2L]
    return(subset_path(object, lambda.id, rho.id))
  }
  ids <- path_ids(object, lambda.id, rho.id)
  subset_path(object, ids[1L], ids[2L])
}

check_index <- function(i, n, what) {
  if (!is.numeric(i) || length(i) != 1L || i < 1 || i > n || i != round(i))
    stop(sQuote(what), " must be an integer between 1 and ", n, call. = FALSE)
}

subset_path <- function(object, h, o) {
  sub5 <- function(a) if (is.null(a)) NULL else a[, , , h, o, drop = FALSE]
  sub3 <- function(a) if (is.null(a)) NULL else a[, h, o, drop = FALSE]
  object$B <- sub5(object$B); object$Tht <- sub5(object$Tht); object$Omega <- sub5(object$Omega)
  object$psi <- object$psi[, , h, o, drop = FALSE]
  object$Q <- list(yx = sub3(object$Q$yx), x = sub3(object$Q$x))
  object$dfB <- sub3(object$dfB); object$dfTht <- sub3(object$dfTht); object$dfOmg <- sub3(object$dfOmg)
  object$nit <- object$nit[, h, o, drop = FALSE]
  object$conv <- object$conv[h, o, drop = FALSE]
  object$lambda <- object$lambda[h]; object$rho <- object$rho[o]; object$nu <- object$nu[o]
  object$nlambda <- 1L; object$nrho <- 1L
  object
}
