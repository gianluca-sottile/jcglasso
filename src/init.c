/* .Call interfaces and registration of the native routines. */
#include "jcglasso.h"
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

SEXP C_estep(SEXP Z, SEXP R, SEXP psi, SEXP Xi, SEXP lo, SEXP up, SEXP zlo, SEXP zup, SEXP pstart, SEXP prow)
{
    int n = nrows(Z), d = ncols(Z), npat = length(pstart) - 1;
    SEXP out = PROTECT(allocVector(VECSXP, 4));
    SEXP Zout = PROTECT(allocMatrix(REALSXP, n, d)), T1 = PROTECT(allocVector(REALSXP, d));
    SEXP T2 = PROTECT(allocMatrix(REALSXP, d, d));
    int conv = jcg_estep(n, d, REAL(Z), INTEGER(R), REAL(psi), REAL(Xi), REAL(lo), REAL(up), REAL(zlo), REAL(zup),
                         INTEGER(pstart), INTEGER(prow), npat, REAL(Zout), REAL(T1), REAL(T2));
    SET_VECTOR_ELT(out, 0, Zout); SET_VECTOR_ELT(out, 1, T1); SET_VECTOR_ELT(out, 2, T2);
    SET_VECTOR_ELT(out, 3, ScalarInteger(conv));
    UNPROTECT(4);
    return out;
}

SEXP C_apg_B(SEXP B, SEXP Sxx, SEXP Sxy, SEXP Tht, SEXP f, SEXP lambda, SEXP alpha, SEXP w, SEXP maxit, SEXP thr)
{
    const int *dim = INTEGER(getAttrib(B, R_DimSymbol));
    int q = dim[0], p = dim[1], K = dim[2], nit = 0;
    SEXP out = PROTECT(allocVector(VECSXP, 3)), Bn = PROTECT(duplicate(B));
    int conv = jcg_apg_B(q, p, K, REAL(Bn), REAL(Sxx), REAL(Sxy), REAL(Tht), REAL(f), asReal(lambda), asReal(alpha),
                         REAL(w), asInteger(maxit), asReal(thr), &nit);
    SET_VECTOR_ELT(out, 0, Bn); SET_VECTOR_ELT(out, 1, ScalarInteger(nit)); SET_VECTOR_ELT(out, 2, ScalarInteger(conv));
    UNPROTECT(2);
    return out;
}

SEXP C_jgl(SEXP S, SEXP f, SEXP rho, SEXP alpha, SEXP w, SEXP Tht, SEXP pendiag, SEXP maxit, SEXP thr)
{
    const int *dim = INTEGER(getAttrib(S, R_DimSymbol));
    int p = dim[0], K = dim[2], nit = 0;
    SEXP out = PROTECT(allocVector(VECSXP, 4)), Tn = PROTECT(duplicate(Tht)), comp = PROTECT(allocVector(INTSXP, p));
    int conv = jcg_jgl(p, K, REAL(S), REAL(f), asReal(rho), asReal(alpha), REAL(w), REAL(Tn), asInteger(pendiag),
                       asInteger(maxit), asReal(thr), INTEGER(comp), &nit);
    SET_VECTOR_ELT(out, 0, Tn); SET_VECTOR_ELT(out, 1, comp); SET_VECTOR_ELT(out, 2, ScalarInteger(nit));
    SET_VECTOR_ELT(out, 3, ScalarInteger(conv));
    UNPROTECT(3);
    return out;
}

static const R_CallMethodDef CallEntries[] = {
    {"C_estep", (DL_FUNC) &C_estep, 10},
    {"C_apg_B", (DL_FUNC) &C_apg_B, 10},
    {"C_jgl", (DL_FUNC) &C_jgl, 9},
    {NULL, NULL, 0}
};

void attribute_visible R_init_jcglasso(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
