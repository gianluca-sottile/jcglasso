/* jcglasso: numerical kernels. The core functions use plain C arrays (column
 * major) and BLAS/LAPACK; the .Call interfaces are in init.c. */
#ifndef JCGLASSO_H
#define JCGLASSO_H

#ifdef JCG_STANDALONE
/* stand-alone build for testing outside R */
#include <stdlib.h>
#include <math.h>
#define F77_CALL(x) x ## _
#define FCONE , (size_t) 1
#define jcg_interrupt() ((void) 0)
#define jcg_alloc(n, type) ((type *) calloc((size_t) (n), sizeof(type)))
#define jcg_free(p) free(p)
#define R_FINITE_(x) isfinite(x)
double jcg_lmills_lower(double z); /* log(dnorm(z) / pnorm(z)) */
double jcg_lmills_upper(double z); /* log(dnorm(z) / (1 - pnorm(z))) */
void dgemm_(const char *, const char *, const int *, const int *, const int *, const double *, const double *,
            const int *, const double *, const int *, const double *, double *, const int *, size_t, size_t);
void dsyrk_(const char *, const char *, const int *, const int *, const double *, const double *, const int *,
            const double *, double *, const int *, size_t, size_t);
void dpotrf_(const char *, const int *, double *, const int *, int *, size_t);
void dpotri_(const char *, const int *, double *, const int *, int *, size_t);
void dsyev_(const char *, const char *, const int *, double *, const int *, double *, double *, const int *, int *,
            size_t, size_t);
#else
#define USE_FC_LEN_T
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#ifndef FCONE
#define FCONE
#endif
#define jcg_interrupt() R_CheckUserInterrupt()
#define jcg_alloc(n, type) ((type *) R_alloc((size_t) (n), sizeof(type)))
#define jcg_free(p) ((void) 0)
#define R_FINITE_(x) R_FINITE(x)
static inline double jcg_lmills_lower(double z) { return dnorm(z, 0.0, 1.0, 1) - pnorm(z, 0.0, 1.0, 1, 1); }
static inline double jcg_lmills_upper(double z) { return dnorm(z, 0.0, 1.0, 1) - pnorm(z, 0.0, 1.0, 0, 1); }
#endif

int jcg_estep(int n, int d, const double *Z, const int *R, const double *psi, const double *Xi, const double *lo,
              const double *up, const double *zlo, const double *zup, const int *pstart, const int *prow, int npat,
              double *Zout, double *T1, double *T2);

int jcg_apg_B(int q, int p, int K, double *B, const double *Sxx, const double *Sxy, const double *Tht,
              const double *f, double lambda, double alpha, const double *w, const double *gw, int maxit, double thr,
              int *nit);

int jcg_jgl(int p, int K, const double *S, const double *f, double rho, double alpha, const double *w,
            const double *gw, double *Tht, int pendiag, int maxit, double thr, int *comp, int *nit);

#endif
