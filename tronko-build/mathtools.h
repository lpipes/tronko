#ifndef _MATHTOOLS_H
#define _MATHTOOLS_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>

/* Returns a vector of the (k+1)-quantiles of a Gamma distribution with
 * parameters alpha, beta. The caller is responsible for freeing the returned
 * vector! */
double *getgammaquantiles(double alpha, double beta, int k);

void initlogfactorial();

/* Eigendecomposes a reversible CTMC rate matrix Q = V D V^{-1}, given its
 * stationary distribution pi (with pi_i * Q_ij = pi_j * Q_ji).
 *
 * Q and pi must be n x n / length n, for some n >= 1.
 *
 * Outputs (caller-allocated, must match Q's dimension n):
 *   eigenvalues = gsl_vector_alloc(n);        // real eigenvalues
 *   V           = gsl_matrix_alloc(n, n);     // eigenvectors of Q, as columns
 *   Vinv        = gsl_matrix_alloc(n, n);     // inverse of V
 *
 * Returns 0 on success, nonzero GSL error code otherwise.
 */
int eigendecompose_gtr(const gsl_matrix *Q, const gsl_vector *pi,
                       gsl_vector *eigenvalues, gsl_matrix *V,
                       gsl_matrix *Vinv);

#endif /* _MATHTOOLS_H_ */
