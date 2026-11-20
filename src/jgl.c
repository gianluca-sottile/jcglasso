/* Joint graphical lasso with group penalty (Danaher et al., 2014):
 *
 *   min sum_k f_k (tr(S_k Tht_k) - log det Tht_k)
 *       + rho [alpha sum_k sum_{i!=j} w_ijk |tht_ijk| + (1 - alpha) sum_{i!=j} (sum_k tht_ijk^2)^(1/2)]
 *
 * The variables are first split into connected components (Danaher et al., 2014,
 * Theorem 2); isolated variables have closed-form solutions, the other components
 * are solved by ADMM with adaptive penalty parameter (Boyd et al., 2011).
 * S, w, Tht: p x p x K (Tht is used as starting value and overwritten).
 * comp: p integers (component labels, 1-based). Returns 0 or 1 (maximum number of
 * iterations reached in some component). */
#include "jcglasso.h"

static int components(int p, int K, const double *S, const double *f, double rho, double alpha, const double *w,
                      int *comp)
{
    int i, j, k, nc = 0;
    int *stack = jcg_alloc(p, int);
    char *adj = jcg_alloc((size_t) p * p, char);
    for (j = 0; j < p; j++)
        for (i = 0; i < p; i++) {
            double s = 0.0;
            if (i == j) { adj[i + (size_t) p * j] = 0; continue; }
            for (k = 0; k < K; k++) {
                size_t id = i + (size_t) p * j + (size_t) p * p * k;
                double v = f[k] * fabs(S[id]) - alpha * rho * w[id];
                if (v > 0.0) s += v * v;
            }
            adj[i + (size_t) p * j] = sqrt(s) > (1.0 - alpha) * rho * w[i + (size_t) p * j];
        }
    for (i = 0; i < p; i++) comp[i] = 0;
    for (i = 0; i < p; i++) {
        if (comp[i]) continue;
        int top = 0;
        nc++;
        comp[i] = nc; stack[top++] = i;
        while (top > 0) {
            int u = stack[--top];
            for (j = 0; j < p; j++)
                if (!comp[j] && adj[u + (size_t) p * j]) { comp[j] = nc; stack[top++] = j; }
        }
    }
    jcg_free(stack); jcg_free(adj);
    return nc;
}

/* ADMM for one block of size m (arrays m x m x K) */
static int admm_block(int m, int K, const double *S, const double *f, double rho, double alpha, const double *w,
                      double *Tht, int pendiag, int maxit, double thr, int *nit)
{
    size_t mm = (size_t) m * m, N = mm * K;
    int i, j, k, it, info, lwork = -1, conv = 1;
    double tau = 2.0, wopt;
    double *Z = jcg_alloc(N, double), *U = jcg_alloc(N, double), *Zold = jcg_alloc(N, double);
    double *A = jcg_alloc(mm, double), *ev = jcg_alloc(m, double);
    double *work;
    F77_CALL(dsyev)("V", "U", &m, A, &m, ev, &wopt, &lwork, &info FCONE FCONE);
    lwork = (int) wopt;
    if (lwork < 3 * m) lwork = 3 * m;
    work = jcg_alloc(lwork, double);
    for (i = 0; i < (int) N; i++) { Z[i] = Tht[i]; U[i] = 0.0; }
    *nit = 0;
    for (it = 1; it <= maxit; it++) {
        jcg_interrupt();
        double r = 0.0, s = 0.0, e1 = 0.0, e2 = 0.0;
        for (i = 0; i < (int) N; i++) Zold[i] = Z[i];
        /* Tht-update: eigen-decomposition of S_k - (tau Z_k - U_k) / f_k */
        for (k = 0; k < K; k++) {
            const double *Sk = S + mm * k;
            double *Zk = Z + mm * k, *Uk = U + mm * k, *Tk = Tht + mm * k;
            for (i = 0; i < (int) mm; i++) A[i] = Sk[i] - (tau * Zk[i] - Uk[i]) / f[k];
            F77_CALL(dsyev)("V", "U", &m, A, &m, ev, work, &lwork, &info FCONE FCONE);
            if (info != 0) { conv = 2; goto done; }
            for (j = 0; j < m; j++) ev[j] = f[k] * (-ev[j] + sqrt(ev[j] * ev[j] + 4.0 * tau / f[k])) / (2.0 * tau);
            for (j = 0; j < m; j++)
                for (i = 0; i <= j; i++) {
                    double v = 0.0;
                    for (int l = 0; l < m; l++) v += A[i + (size_t) m * l] * ev[l] * A[j + (size_t) m * l];
                    Tk[i + (size_t) m * j] = v; Tk[j + (size_t) m * i] = v;
                }
        }
        /* Z-update: sparse group lasso proximal operator (off-diagonal entries) */
        for (j = 0; j < m; j++)
            for (i = 0; i < m; i++) {
                double nrm = 0.0;
                for (k = 0; k < K; k++) {
                    size_t id = i + (size_t) m * j + mm * k;
                    double v = Tht[id] + U[id] / tau, t = (i == j && !pendiag) ? 0.0 : alpha * rho * w[id] / tau;
                    v = v > t ? v - t : (v < -t ? v + t : 0.0);
                    Z[id] = v; nrm += v * v;
                }
                if (K > 1 && !(i == j && !pendiag) && alpha < 1.0) {
                    nrm = sqrt(nrm);
                    for (k = 0; k < K; k++) {
                        size_t id = i + (size_t) m * j + mm * k;
                        double g = (1.0 - alpha) * rho * w[id] / tau, fac = nrm > 0.0 ? 1.0 - g / nrm : 0.0;
                        Z[id] = fac > 0.0 ? Z[id] * fac : 0.0;
                    }
                }
            }
        /* dual update and stopping rule */
        for (k = 0; k < K; k++) {
            double rk = 0.0, sk = 0.0, nt = 0.0, nz = 0.0, nu = 0.0;
            for (i = 0; i < (int) mm; i++) {
                size_t id = i + mm * k;
                U[id] += tau * (Tht[id] - Z[id]);
                rk += (Tht[id] - Z[id]) * (Tht[id] - Z[id]);
                sk += (Z[id] - Zold[id]) * (Z[id] - Zold[id]);
                nt += Tht[id] * Tht[id]; nz += Z[id] * Z[id]; nu += U[id] * U[id];
            }
            r += f[k] * sqrt(rk); s += f[k] * tau * sqrt(sk);
            e1 += f[k] * (m * thr + thr * (sqrt(nt) > sqrt(nz) ? sqrt(nt) : sqrt(nz)));
            e2 += f[k] * (m * thr + thr * sqrt(nu));
        }
        *nit = it;
        if (r < e1 && s < e2) { conv = 0; break; }
        if (r > 10.0 * s) tau *= 2.0; else if (s > 10.0 * r) tau /= 2.0;
    }
done:
    for (i = 0; i < (int) N; i++) Tht[i] = Z[i];
    jcg_free(Z); jcg_free(U); jcg_free(Zold); jcg_free(A); jcg_free(ev); jcg_free(work);
    return conv;
}

int jcg_jgl(int p, int K, const double *S, const double *f, double rho, double alpha, const double *w,
            double *Tht, int pendiag, int maxit, double thr, int *comp, int *nit)
{
    size_t pp = (size_t) p * p;
    int c, i, j, k, nc, conv = 0, nnit;
    int *idx = jcg_alloc(p, int);
    nc = components(p, K, S, f, rho, alpha, w, comp);
    *nit = 0;
    double *out = jcg_alloc(pp * K, double);
    for (i = 0; i < (int) (pp * K); i++) out[i] = 0.0;
    for (c = 1; c <= nc; c++) {
        int m = 0;
        for (i = 0; i < p; i++) if (comp[i] == c) idx[m++] = i;
        if (m == 1) {
            int v = idx[0];
            for (k = 0; k < K; k++) {
                size_t id = v + (size_t) p * v + pp * k;
                double s = S[id] + (pendiag ? alpha * rho * w[id] / f[k] : 0.0);
                out[id] = 1.0 / s;
            }
            continue;
        }
        size_t mm = (size_t) m * m;
        double *Sb = jcg_alloc(mm * K, double), *wb = jcg_alloc(mm * K, double), *Tb = jcg_alloc(mm * K, double);
        for (k = 0; k < K; k++)
            for (j = 0; j < m; j++)
                for (i = 0; i < m; i++) {
                    size_t ib = i + (size_t) m * j + mm * k, id = idx[i] + (size_t) p * idx[j] + pp * k;
                    Sb[ib] = S[id]; wb[ib] = w[id]; Tb[ib] = Tht[id];
                }
        int cc = admm_block(m, K, Sb, f, rho, alpha, wb, Tb, pendiag, maxit, thr, &nnit);
        if (cc > conv) conv = cc;
        *nit += nnit;
        for (k = 0; k < K; k++)
            for (j = 0; j < m; j++)
                for (i = 0; i < m; i++) out[idx[i] + (size_t) p * idx[j] + pp * k] = Tb[i + (size_t) m * j + mm * k];
        jcg_free(Sb); jcg_free(wb); jcg_free(Tb);
    }
    for (i = 0; i < (int) (pp * K); i++) Tht[i] = out[i];
    jcg_free(idx); jcg_free(out);
    return conv;
}
