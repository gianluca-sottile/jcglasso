/* E-step of the EM algorithm for one class.
 *
 * For each pattern of unobserved values (rows prow[pstart[i]], ..., prow[pstart[i+1]-1]
 * share the same codes R), the conditional distribution of the unobserved entries m
 * given the observed ones o is N(psi_m - Xi_mm^{-1} Xi_mo (z_o - psi_o), Xi_mm^{-1}).
 * Missing-at-random entries are replaced by their conditional means (clipped to
 * [zlo, zup]); censored entries by the means of the univariate truncated normal
 * distributions, with the approximation of Guo et al. (2015) for the second moments.
 *
 * Input:  Z (n x d, values of unobserved entries are ignored), R (codes: 0 observed,
 *         -1 left censored, 1 right censored, 9 missing at random), psi (d), Xi (d x d).
 * Output: Zout (imputed data), T1 = colSums(Zout), T2 = Zout' Zout + diag(conditional
 *         variances summed over the rows).
 * Returns 0, or 1 if a conditional covariance matrix is not positive definite. */
#include "jcglasso.h"

int jcg_estep(int n, int d, const double *Z, const int *R, const double *psi, const double *Xi, const double *lo,
              const double *up, const double *zlo, const double *zup, const int *pstart, const int *prow, int npat,
              double *Zout, double *T1, double *T2)
{
    int i, j, jj, h, r, info = 0;
    double *addv = jcg_alloc(d, double);
    int *mi = jcg_alloc(d, int), *oi = jcg_alloc(d, int);
    double *Smm = jcg_alloc((size_t) d * d, double), *Bmo = jcg_alloc((size_t) d * d, double);
    double *mu = jcg_alloc(d, double), *res = jcg_alloc(d, double);

    for (i = 0; i < n * d; i++) Zout[i] = Z[i];
    for (j = 0; j < d; j++) addv[j] = 0.0;

    for (h = 0; h < npat; h++) {
        jcg_interrupt();
        int r0 = prow[pstart[h]], nm = 0, no = 0;
        for (j = 0; j < d; j++) {
            if (R[r0 + (size_t) n * j] != 0) mi[nm++] = j; else oi[no++] = j;
        }
        if (nm == 0) continue;
        /* Smm = (Xi_mm)^{-1} */
        for (jj = 0; jj < nm; jj++)
            for (i = 0; i < nm; i++) Smm[i + (size_t) nm * jj] = Xi[mi[i] + (size_t) d * mi[jj]];
        F77_CALL(dpotrf)("L", &nm, Smm, &nm, &info FCONE);
        if (info != 0) { info = 1; goto done; }
        F77_CALL(dpotri)("L", &nm, Smm, &nm, &info FCONE);
        if (info != 0) { info = 1; goto done; }
        for (jj = 0; jj < nm; jj++)
            for (i = jj + 1; i < nm; i++) Smm[jj + (size_t) nm * i] = Smm[i + (size_t) nm * jj];
        /* Bmo = Smm Xi_mo (nm x no) */
        for (jj = 0; jj < no; jj++)
            for (i = 0; i < nm; i++) {
                double s = 0.0;
                for (r = 0; r < nm; r++) s += Smm[i + (size_t) nm * r] * Xi[mi[r] + (size_t) d * oi[jj]];
                Bmo[i + (size_t) nm * jj] = s;
            }
        for (r = pstart[h]; r < pstart[h + 1]; r++) {
            int row = prow[r];
            for (jj = 0; jj < no; jj++) res[jj] = Z[row + (size_t) n * oi[jj]] - psi[oi[jj]];
            for (i = 0; i < nm; i++) {
                double s = psi[mi[i]];
                for (jj = 0; jj < no; jj++) s -= Bmo[i + (size_t) nm * jj] * res[jj];
                mu[i] = s;
            }
            for (i = 0; i < nm; i++) {
                int jv = mi[i], code = R[row + (size_t) n * jv];
                double m = mu[i], v = Smm[i + (size_t) nm * i], s = sqrt(v), tm, tv;
                if (code == 9) {
                    tm = m < zlo[jv] ? zlo[jv] : (m > zup[jv] ? zup[jv] : m);
                    tv = v;
                } else {
                    double z, ratio, bnd;
                    if (code == -1) {
                        z = (lo[jv] - m) / s;
                        ratio = exp(jcg_lmills_lower(z));
                        tm = m - s * ratio;
                        tv = v * (1.0 - z * ratio - ratio * ratio);
                        bnd = zlo[jv] < lo[jv] ? zlo[jv] : lo[jv];
                        if (!R_FINITE_(tm) || !R_FINITE_(tv) || tv <= 0.0 || tm < bnd) { tm = bnd; tv = v; }
                    } else {
                        z = (up[jv] - m) / s;
                        ratio = exp(jcg_lmills_upper(z));
                        tm = m + s * ratio;
                        tv = v * (1.0 + z * ratio - ratio * ratio);
                        bnd = zup[jv] > up[jv] ? zup[jv] : up[jv];
                        if (!R_FINITE_(tm) || !R_FINITE_(tv) || tv <= 0.0 || tm > bnd) { tm = bnd; tv = v; }
                    }
                }
                Zout[row + (size_t) n * jv] = tm;
                addv[jv] += tv;
            }
        }
    }
    /* sufficient statistics */
    for (j = 0; j < d; j++) {
        double s = 0.0;
        for (i = 0; i < n; i++) s += Zout[i + (size_t) n * j];
        T1[j] = s;
    }
    {
        double one = 1.0, zero = 0.0;
        F77_CALL(dsyrk)("U", "T", &d, &n, &one, Zout, &n, &zero, T2, &d FCONE FCONE);
        for (j = 0; j < d; j++) {
            for (i = j + 1; i < d; i++) T2[i + (size_t) d * j] = T2[j + (size_t) d * i];
            T2[j + (size_t) d * j] += addv[j];
        }
    }
done:
    jcg_free(addv); jcg_free(mi); jcg_free(oi); jcg_free(Smm); jcg_free(Bmo); jcg_free(mu); jcg_free(res);
    return info;
}
