# Without compiled code (e.g., when testing in WebAssembly) use the R engine.
has_C <- exists("C_estep", envir = asNamespace("jcglasso"), inherits = FALSE)
if (!has_C) options(jcglasso.engine = "R")
ns <- asNamespace("jcglasso")

small_sim <- function(seed = 1, n = 60, p = 6, q = 3, K = 2) {
  set.seed(seed)
  rjcggm(n = n, p = p, q = q, K = K, ncompPrecision = 2, ncompB = 2)
}
