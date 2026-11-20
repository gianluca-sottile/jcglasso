# Dispatch of the numerical kernels to the compiled code (default) or to the
# reference R implementations (option jcglasso.engine = "R", used for testing).

use_C <- function() !identical(getOption("jcglasso.engine"), "R")

# Rows grouped by pattern of unobserved values (0-based, for the C code): rows
# row[start[i] + 1], ..., row[start[i + 1]] share the same codes.
make_patterns <- function(R) {
  un <- which(rowSums(R != 0L) > 0L)
  if (!length(un)) return(list(start = 0L, row = integer(0)))
  key <- apply(R[un, , drop = FALSE], 1L, paste, collapse = ",")
  grp <- split(un, factor(key, levels = unique(key)))
  list(start = as.integer(c(0L, cumsum(lengths(grp)))), row = as.integer(unlist(grp, use.names = FALSE) - 1L))
}

k_estep <- function(Z, R, psi, Xi, lo, up, zlo, zup, pat = make_patterns(R)) {
  if (use_C()) {
    Z0 <- Z; Z0[R != 0L] <- 0
    storage.mode(Z0) <- "double"; storage.mode(R) <- "integer"
    out <- .Call(C_estep, Z0, R, as.double(psi), Xi, as.double(lo), as.double(up), as.double(zlo), as.double(zup),
                 pat$start, pat$row)
    names(out) <- c("Z", "T1", "T2", "conv")
    if (out$conv != 0L) stop("the E-step failed (singular conditional covariance matrix)", call. = FALSE)
    out
  } else estep_R(Z, R, psi, Xi, lo, up, zlo, zup, pat)
}

k_apg <- function(B, Sxx, Sxy, Tht, f, lambda, alpha, w, gw, maxit, thr) {
  if (use_C()) {
    storage.mode(B) <- storage.mode(w) <- storage.mode(gw) <- "double"
    out <- .Call(C_apg_B, B, Sxx, Sxy, Tht, as.double(f), as.double(lambda), as.double(alpha), w, gw,
                 as.integer(maxit), as.double(thr))
    names(out) <- c("B", "nit", "conv")
    dim(out$B) <- dim(B)
    out
  } else apg_B_R(B, Sxx, Sxy, Tht, f, lambda, alpha, w, maxit, thr, gw)
}

k_jgl <- function(S, f, rho, alpha, w, gw, Tht, pendiag, maxit, thr) {
  if (use_C()) {
    storage.mode(S) <- storage.mode(w) <- storage.mode(gw) <- storage.mode(Tht) <- "double"
    out <- .Call(C_jgl, S, as.double(f), as.double(rho), as.double(alpha), w, gw, Tht, as.integer(pendiag),
                 as.integer(maxit), as.double(thr))
    names(out) <- c("Tht", "comp", "nit", "conv")
    dim(out$Tht) <- dim(S)
    out
  } else jgl_R(S, f, rho, alpha, w, Tht, pendiag, maxit, thr, gw)
}
