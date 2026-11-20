#' Graphs of a Fitted jcglasso Model
#'
#' Creates, for each class, the graphs estimated by a `"jcglasso"` (or `"jcggm"`)
#' model, as objects of package 'igraph': the undirected graph of the responses
#' given the predictors (`Gyy`, from \eqn{\Theta_k}), the directed graph of the
#' effects of the predictors on the responses (`Gxy`, from \eqn{B_k}) and the
#' undirected graph of the predictors (`Gxx`, from \eqn{\Omega_k}).
#'
#' @param object an object of class `"jcglasso"`.
#' @param GoF,lambda.id,rho.id used to select the model, see [select_jcglasso()].
#' @param weighted logical: if `TRUE`, the edges have weights (partial
#'   correlations for `Gyy` and `Gxx`, regression coefficients for `Gxy`).
#' @param ... further arguments passed to [select_jcglasso()].
#' @return `to_graph()` returns an object of class `"jcglasso2igraph"`, a list
#'   with components `Gyy`, `Gxy` and `Gxx` (lists of `K` graphs; `Gxy` and `Gxx`
#'   are `NULL` without predictors). `getGraph()` returns a list of `K` graphs;
#'   the combined graphs (`"conditional"`, `"bipartite"`) are directed, with the
#'   edge attribute `directed` equal to `FALSE` for the undirected edges of
#'   \eqn{\Theta_k} and \eqn{\Omega_k}.
#' @details Package 'igraph' is required.
#' @seealso [jcglasso()]
#' @examplesIf requireNamespace("igraph", quietly = TRUE)
#' set.seed(1)
#' sim <- rjcggm(n = 50, p = 6, q = 2, K = 2, ncompPrecision = 2)
#' fit <- jcglasso(data = sim$Z, nrho = 3, nlambda = 3, rho.min.ratio = 0.3,
#'                 lambda.min.ratio = 0.3)
#' gr <- to_graph(fit, lambda.id = 3, rho.id = 3)
#' gr
#' plot(gr, type = "conditional")
#' @export
to_graph <- function(object, ...) UseMethod("to_graph")

#' @rdname to_graph
#' @export
to_graph.default <- function(object, ...) delegate("to_graph", object, ...)

#' @rdname to_graph
#' @export
to_graph.jcglasso <- function(object, GoF = AIC, lambda.id, rho.id, weighted = FALSE, ...) {
  if (!requireNamespace("igraph", quietly = TRUE)) stop("package 'igraph' is required", call. = FALSE)
  object <- select_jcglasso(object, GoF = GoF, lambda.id = lambda.id, rho.id = rho.id, ...)
  K <- length(object$nobs); q <- object$npred
  und <- function(M) {
    pc <- -stats::cov2cor(M); diag(pc) <- 0
    A <- (M != 0) * 1; diag(A) <- 0
    g <- igraph::graph_from_adjacency_matrix(if (weighted) pc * A else A, mode = "undirected",
                                             weighted = if (weighted) TRUE else NULL, diag = FALSE)
    igraph::V(g)$type <- "response"
    g
  }
  Gyy <- lapply(seq_len(K), function(k) und(object$Tht[, , k, 1L, 1L]))
  Gxy <- Gxx <- NULL
  if (q) {
    Gxx <- lapply(seq_len(K), function(k) und(object$Omega[, , k, 1L, 1L]))
    Gxy <- lapply(seq_len(K), function(k) {
      B <- object$B[-1L, , k, 1L, 1L]; B <- matrix(B, q, dimnames = list(dimnames(object$B)[[1L]][-1L], dimnames(object$B)[[2L]]))
      e <- which(B != 0, arr.ind = TRUE)
      d <- data.frame(from = rownames(B)[e[, 1L]], to = colnames(B)[e[, 2L]])
      if (weighted) d$weight <- B[e]
      v <- data.frame(name = c(colnames(B), rownames(B)), type = c(rep("response", ncol(B)), rep("predictor", nrow(B))))
      igraph::graph_from_data_frame(d, directed = TRUE, vertices = v)
    })
    for (k in seq_len(K)) igraph::V(Gxx[[k]])$type <- "predictor"
    names(Gxy) <- names(Gxx) <- object$Z$classes
  }
  names(Gyy) <- object$Z$classes
  structure(list(Gyy = Gyy, Gxy = Gxy, Gxx = Gxx), class = "jcglasso2igraph")
}

#' @rdname to_graph
#' @param x an object of class `"jcglasso2igraph"`.
#' @param type which graph: `"Gyy"`, `"Gxy"`, `"Gxx"`, `"conditional"` (`Gyy`
#'   and `Gxy` together) or `"bipartite"` (all of them).
#' @export
getGraph <- function(x, ...) UseMethod("getGraph")

#' @rdname to_graph
#' @export
getGraph.default <- function(x, ...) delegate("getGraph", x, ...)

#' @rdname to_graph
#' @export
getGraph.jcglasso2igraph <- function(x, type = c("Gyy", "Gxy", "Gxx", "conditional", "bipartite"), ...) {
  type <- match.arg(type)
  if (is.null(x$Gxy) && type != "Gyy") stop("the model has no predictors: use type = \"Gyy\"", call. = FALSE)
  if (type %in% c("Gyy", "Gxy", "Gxx")) return(x[[type]])
  out <- lapply(seq_along(x$Gyy), function(k) {
    gs <- if (type == "conditional") list(x$Gyy[[k]], x$Gxy[[k]]) else list(x$Gyy[[k]], x$Gxy[[k]], x$Gxx[[k]])
    e <- do.call(rbind, lapply(gs, function(g) {
      d <- igraph::as_data_frame(g, what = "edges")[, c("from", "to"), drop = FALSE]
      d$directed <- rep(igraph::is_directed(g), nrow(d))
      d
    }))
    v <- do.call(rbind, lapply(gs, function(g) igraph::as_data_frame(g, what = "vertices")[, c("name", "type"),
                                                                                         drop = FALSE]))
    v <- v[!duplicated(v$name), , drop = FALSE]
    igraph::graph_from_data_frame(e, directed = TRUE, vertices = v)
  })
  names(out) <- names(x$Gyy)
  out
}

#' @export
print.jcglasso2igraph <- function(x, ...) {
  cat("Graphs of a jcglasso model (", length(x$Gyy), " classes)\n", sep = "")
  for (k in names(x$Gyy)) {
    cat(" ", k, ": Gyy ", igraph::ecount(x$Gyy[[k]]), " edges",
        if (!is.null(x$Gxy)) paste0(", Gxy ", igraph::ecount(x$Gxy[[k]]), " edges, Gxx ",
                                    igraph::ecount(x$Gxx[[k]]), " edges"), "\n", sep = "")
  }
  invisible(x)
}

#' @rdname to_graph
#' @param which the classes to plot (default all).
#' @export
plot.jcglasso2igraph <- function(x, type, which, ...) {
  if (missing(type)) type <- if (is.null(x$Gxy)) "Gyy" else "conditional"
  gr <- getGraph(x, type = type)
  if (missing(which)) which <- seq_along(gr)
  op <- graphics::par(mfrow = grDevices::n2mfrow(length(which)))
  on.exit(graphics::par(op))
  for (k in which) {
    g <- gr[[k]]
    col <- ifelse(igraph::V(g)$type == "predictor", "lightsalmon", "lightblue")
    am <- if (!is.null(igraph::E(g)$directed)) ifelse(igraph::E(g)$directed, 2, 0) else
      if (igraph::is_directed(g)) 2 else 0
    graphics::plot(g, vertex.color = col, edge.arrow.size = 0.3, edge.arrow.mode = am, main = names(gr)[k], ...)
  }
  invisible(x)
}
