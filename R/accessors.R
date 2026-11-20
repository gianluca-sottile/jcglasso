# Accessor functions. They are S3 generics whose default methods call the
# homonymous functions of 'cglasso' (if installed), so that loading both packages
# does not break the functions of 'cglasso'.

delegate <- function(fun, object, ...) {
  if (requireNamespace("cglasso", quietly = TRUE) && exists(fun, envir = asNamespace("cglasso"), inherits = FALSE))
    return(get(fun, envir = asNamespace("cglasso"))(object, ...))
  stop("no ", fun, " method for objects of class ", sQuote(class(object)[1L]), call. = FALSE)
}

#' Accessors for datajcggm Objects
#'
#' Extract the components of a `"datajcggm"` object (or of the data used by a
#' `"jcglasso"` model).
#'
#' @param x,object an object of class `"datajcggm"` or `"jcglasso"`.
#' @param name which matrices: `"Y"` (responses), `"X"` (predictors) or `"both"`.
#' @param ... not used.
#' @return
#' * `getMatrix()`: a list of `K` matrices (for `"both"`, a list of lists with
#'   components `Y` and `X`);
#' * `event()`: a list of `K` matrices of codes (`0` observed, `-1` left censored,
#'   `1` right censored, `9` missing);
#' * `lower()`, `upper()`: the censoring values;
#' * `ColMeans()`, `ColVars()`: matrices of marginal means and variances (one
#'   column for each class), estimated under a censored Gaussian model;
#' * `nresp()`, `npred()`: the numbers of responses and predictors;
#' * `rowNames()`, `colNames()`: the names of the observations of each class and
#'   of the variables.
#' @seealso [datajcggm()]
#' @examples
#' set.seed(1)
#' sim <- rjcggm(n = 30, p = 4, q = 2, K = 2, ncompPrecision = 2)
#' str(getMatrix(sim$Z, "X"))
#' table(event(sim$Z)$class1)
#' ColMeans(sim$Z)
#' nresp(sim$Z); npred(sim$Z)
#' @name accessors
NULL

#' @rdname accessors
#' @export
getMatrix <- function(x, ...) UseMethod("getMatrix")
#' @rdname accessors
#' @export
getMatrix.default <- function(x, ...) delegate("getMatrix", x, ...)
#' @rdname accessors
#' @export
getMatrix.datajcggm <- function(x, name = c("Y", "X", "both"), ...) {
  name <- match.arg(name)
  switch(name, Y = x$Y, X = x$X,
         both = stats::setNames(lapply(seq_along(x$classes), function(k) list(Y = x$Y[[k]], X = x$X[[k]])),
                                x$classes))
}
#' @rdname accessors
#' @export
getMatrix.jcglasso <- function(x, name = c("Y", "X", "both"), ...) getMatrix(x$Z, name = name, ...)

#' @rdname accessors
#' @export
event <- function(x, ...) UseMethod("event")
#' @rdname accessors
#' @export
event.default <- function(x, ...) delegate("event", x, ...)
#' @rdname accessors
#' @export
event.datajcggm <- function(x, ...) x$R
#' @rdname accessors
#' @export
event.jcglasso <- function(x, ...) x$Z$R

#' @rdname accessors
#' @export
lower <- function(x, ...) UseMethod("lower")
#' @rdname accessors
#' @export
lower.default <- function(x, ...) delegate("lower", x, ...)
#' @rdname accessors
#' @export
lower.datajcggm <- function(x, ...) x$lo
#' @rdname accessors
#' @export
upper <- function(x, ...) UseMethod("upper")
#' @rdname accessors
#' @export
upper.default <- function(x, ...) delegate("upper", x, ...)
#' @rdname accessors
#' @export
upper.datajcggm <- function(x, ...) x$up

#' @rdname accessors
#' @export
ColMeans <- function(x, ...) UseMethod("ColMeans")
#' @rdname accessors
#' @export
ColMeans.default <- function(x, ...) delegate("ColMeans", x, ...)
#' @rdname accessors
#' @export
ColMeans.datajcggm <- function(x, ...) x$zm
#' @rdname accessors
#' @export
ColVars <- function(x, ...) UseMethod("ColVars")
#' @rdname accessors
#' @export
ColVars.default <- function(x, ...) delegate("ColVars", x, ...)
#' @rdname accessors
#' @export
ColVars.datajcggm <- function(x, ...) x$zv

#' @rdname accessors
#' @export
nresp <- function(object, ...) UseMethod("nresp")
#' @rdname accessors
#' @export
nresp.default <- function(object, ...) delegate("nresp", object, ...)
#' @rdname accessors
#' @export
nresp.datajcggm <- function(object, ...) ncol(object$Y[[1L]])
#' @rdname accessors
#' @export
nresp.jcglasso <- function(object, ...) object$nresp
#' @rdname accessors
#' @export
npred <- function(object, ...) UseMethod("npred")
#' @rdname accessors
#' @export
npred.default <- function(object, ...) delegate("npred", object, ...)
#' @rdname accessors
#' @export
npred.datajcggm <- function(object, ...) if (is.null(object$X)) 0L else ncol(object$X[[1L]])
#' @rdname accessors
#' @export
npred.jcglasso <- function(object, ...) object$npred

#' @rdname accessors
#' @export
rowNames <- function(x, ...) UseMethod("rowNames")
#' @rdname accessors
#' @export
rowNames.default <- function(x, ...) delegate("rowNames", x, ...)
#' @rdname accessors
#' @export
rowNames.datajcggm <- function(x, ...) lapply(x$Y, rownames)
#' @rdname accessors
#' @export
colNames <- function(x, ...) UseMethod("colNames")
#' @rdname accessors
#' @export
colNames.default <- function(x, ...) delegate("colNames", x, ...)
#' @rdname accessors
#' @export
colNames.datajcggm <- function(x, ...) list(Y = colnames(x$Y[[1L]]), X = if (!is.null(x$X)) colnames(x$X[[1L]]))

#' @export
nobs.jcglasso <- function(object, ...) object$nobs
