#include "mathtools.h"
#include "global.h"
#include <math.h>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_cdf.h>

/* See .h file for documentation. */
double *getgammaquantiles(double alpha, double beta, int k) {
    double *statevec = malloc(k * sizeof *statevec);
    for (int i = 0; i < k; i++)
        statevec[i] =
            gsl_cdf_gamma_Pinv((i + 1) / (double)(k + 1), alpha, 1.0 / beta);

    return statevec;
}

void initlogfactorial(void) {
    int i;

    Logfactorial[0] = 0.0;
    Logfactorial[1] = 0.0;
    for (i = 2; i < MAXNUMBEROFINDINSPECIES; i++)
        Logfactorial[i] = Logfactorial[i - 1] + log((double)i);
}

/* See .h file for documentation. */
int eigendecompose_gtr(const gsl_matrix *Q, const gsl_vector *pi,
                       gsl_vector *eigenvalues, gsl_matrix *V,
                       gsl_matrix *Vinv) {

    int n = (int)Q->size1;

#if 0
    /* Use for debugging. */
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (fabs(gsl_vector_get(pi, i) * gsl_matrix_get(Q, i, j) -
                     gsl_vector_get(pi, j) * gsl_matrix_get(Q, j, i)) > 1e-8)
                printf("WARNING: Q does not satisfy detailed balance condition "
                       "for given pi (%s:%d)\n",
                       __FILE__, __LINE__);
#endif

    double *sqrt_pi = malloc(n * sizeof(double));
    for (int i = 0; i < n; ++i) sqrt_pi[i] = sqrt(gsl_vector_get(pi, i));

    /* Symmetrize: S = D^{1/2} Q D^{-1/2}, S_ij = sqrt(pi_i/pi_j) * Q_ij */
    gsl_matrix *S = gsl_matrix_alloc(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            gsl_matrix_set(S, i, j,
                           gsl_matrix_get(Q, i, j) * sqrt_pi[i] / sqrt_pi[j]);

    gsl_matrix *Vs = gsl_matrix_alloc(n, n); /* orthogonal eigenvectors of S */
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(n);
    int status = gsl_eigen_symmv(S, eigenvalues, Vs, w);
    gsl_eigen_symmv_free(w);
    gsl_matrix_free(S);
    if (status) {
        gsl_matrix_free(Vs);
        free(sqrt_pi);
        return status;
    }

    gsl_eigen_symmv_sort(eigenvalues, Vs, GSL_EIGEN_SORT_ABS_DESC);

    /* Recover V = D^{-1/2} Vs  and  Vinv = Vs^T D^{1/2}  (D = diag(pi)) */
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) {
            double vs_ik = gsl_matrix_get(Vs, i, k);
            gsl_matrix_set(V, i, k, vs_ik / sqrt_pi[i]);
            gsl_matrix_set(Vinv, k, i,
                           vs_ik * sqrt_pi[i]); /* = Vs^T scaled by columns */
        }
    }

    gsl_matrix_free(Vs);
    free(sqrt_pi);
    return 0;
}
