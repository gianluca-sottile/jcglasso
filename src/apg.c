/* Accelerated proximal gradient algorithm for the B sub-problem of the M-step:
 *
 *   min_B  (1/2) sum_k f_k tr(Tht_k (B_k' Sxx_k B_k - 2 B_k' Sxy_k)) + lambda * SGL(B)
 *
 * SGL(B) = alpha sum_k sum_{hj} w_hjk |b_hjk| + (1 - alpha) sum_{hj} gw_hj (sum_k b_hjk^2)^(1/2)
 * (sparse group lasso across the K classes). Barzilai-Borwein initial step, adaptive
 * restart and TFOCS-style step-size update (as in the Fortran subroutine 'apg' of the
 * development version of 'cglasso'). The algorithm stops when the norm of the gradient
 * mapping is smaller than thr times the norm of the initial gradient (or lambda). B, Sxy, w: q x p x K; Sxx: q x q x K; Tht: p x p x K.
 * Returns 0 (converged) or 1 (maximum number of iterations). */
#include "jcglasso.h"

static void grad_B(int q, int p, int K, const double *B, const double *Sxx, const double *Sxy, const double *Tht,
                   const double *f, double *G, double *tmp)
{
    int k, i, qp = q * p;
    double one = 1.0, zero = 0.0;
    for (k = 0; k < K; k++) {
        const double *Bk = B + (size_t) qp * k, *Sxk = Sxx + (size_t) q * q * k, *Syk = Sxy + (size_t) qp * k;
        const double *Tk = Tht + (size_t) p * p * k;
        /* tmp = Sxx_k B_k - Sxy_k */
        for (i = 0; i < qp; i++) tmp[i] = -Syk[i];
        F77_CALL(dgemm)("N", "N", &q, &p, &q, &one, Sxk, &q, Bk, &q, &one, tmp, &q FCONE FCONE);
        /* G_k = f_k tmp Tht_k */
        double fk = f[k];
        F77_CALL(dgemm)("N", "N", &q, &p, &p, &fk, tmp, &q, Tk, &p, &zero, G + (size_t) qp * k, &q FCONE FCONE);
    }
}

static void prox_sgl(int q, int p, int K, double *x, double a, double g, const double *w, const double *gw)
{
    int i, k, qp = q * p;
    for (i = 0; i < qp * K; i++) {
        double t = a * w[i], v = x[i];
        x[i] = v > t ? v - t : (v < -t ? v + t : 0.0);
    }
    if (K > 1 && g > 0.0) {
        for (i = 0; i < qp; i++) {
            double nrm = 0.0;
            for (k = 0; k < K; k++) nrm += x[i + (size_t) qp * k] * x[i + (size_t) qp * k];
            nrm = sqrt(nrm);
            for (k = 0; k < K; k++) {
                double fac = nrm > 0.0 ? 1.0 - g * gw[i] / nrm : 0.0;
                x[i + (size_t) qp * k] = fac > 0.0 ? x[i + (size_t) qp * k] * fac : 0.0;
            }
        }
    }
}

int jcg_apg_B(int q, int p, int K, double *B, const double *Sxx, const double *Sxy, const double *Tht,
              const double *f, double lambda, double alpha, const double *w, const double *gw, int maxit, double thr,
              int *nit)
{
    int i, it, N = q * p * K, conv = 1;
    double *x = B, *y = jcg_alloc(N, double), *g = jcg_alloc(N, double), *xold = jcg_alloc(N, double);
    double *yold = jcg_alloc(N, double), *gold = jcg_alloc(N, double), *tmp = jcg_alloc(q * p, double);
    double t, theta = 1.0, nrmg = 0.0, a1, a2, scale;

    *nit = 0;
    for (i = 0; i < N; i++) y[i] = x[i];
    grad_B(q, p, K, y, Sxx, Sxy, Tht, f, g, tmp);
    for (i = 0; i < N; i++) nrmg += g[i] * g[i];
    nrmg = sqrt(nrmg);
    if (nrmg == 0.0) { conv = 0; goto done; }
    /* scale of the stopping rule: norm of the initial gradient (or lambda) */
    scale = nrmg > lambda ? nrmg : lambda;
    /* Barzilai-Borwein initial step */
    t = 1.0 / nrmg;
    for (i = 0; i < N; i++) xold[i] = x[i] - t * g[i];
    grad_B(q, p, K, xold, Sxx, Sxy, Tht, f, gold, tmp);
    a1 = 0.0; a2 = 0.0;
    for (i = 0; i < N; i++) { double dg = g[i] - gold[i]; a1 += (x[i] - xold[i]) * dg; a2 += dg * dg; }
    if (a2 > 0.0) t = fabs(a1 / a2);

    for (it = 1; it <= maxit; it++) {
        jcg_interrupt();
        double num = 0.0, den = 0.0, nx = 0.0, prod = 0.0;
        for (i = 0; i < N; i++) { xold[i] = x[i]; yold[i] = y[i]; x[i] = y[i] - t * g[i]; }
        prox_sgl(q, p, K, x, t * lambda * alpha, t * lambda * (1.0 - alpha), w, gw);
        for (i = 0; i < N; i++) { double dd = y[i] - x[i]; num += dd * dd; nx += x[i] * x[i]; }
        *nit = it;
        /* stopping rule on the gradient mapping (y - x) / t, relative to the scale of the
         * problem, or no progress at machine precision */
        if (sqrt(num) / t <= thr * scale || sqrt(num) <= 1e-14 * (sqrt(nx) > 1.0 ? sqrt(nx) : 1.0)) {
            conv = 0; break;
        }
        theta = 2.0 / (1.0 + sqrt(1.0 + 4.0 / (theta * theta)));
        for (i = 0; i < N; i++) prod += (y[i] - x[i]) * (x[i] - xold[i]);
        if (prod > 0.0) {                        /* adaptive restart */
            for (i = 0; i < N; i++) { x[i] = xold[i]; y[i] = x[i]; }
            theta = 1.0;
        } else {
            for (i = 0; i < N; i++) y[i] = x[i] + (1.0 - theta) * (x[i] - xold[i]);
        }
        for (i = 0; i < N; i++) gold[i] = g[i];
        grad_B(q, p, K, y, Sxx, Sxy, Tht, f, g, tmp);
        num = 0.0; den = 0.0;
        for (i = 0; i < N; i++) { double dy = y[i] - yold[i]; num += dy * dy; den += dy * (gold[i] - g[i]); }
        den = fabs(den);
        double that = den > 0.0 ? 0.5 * num / den : t;
        double lo = 0.5 * t > that ? 0.5 * t : that;
        t = 1.01 * t < lo ? 1.01 * t : lo;
    }
done:
    jcg_free(y); jcg_free(g); jcg_free(xold); jcg_free(yold); jcg_free(gold); jcg_free(tmp);
    return conv;
}
