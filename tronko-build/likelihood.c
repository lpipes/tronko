#include <string.h>
#include <assert.h>
#include "likelihood.h"
#include "mathtools.h"

/* Two transition matrices for nucleotide substitutions over a given time
 * interval (will be used in calculations below). */
static double PMATnc[2][4][5];

/* Eigendecomposition of a GTR rate matrix. Only used in
 * maketransitionmatrixnc_Arr and inittransitionmatrixnc. */
static double V[4][4], eigenvalues[4], Vinv[4][4];

/* Compute the matrix e^{Qt}, where Q is the rate matrix constructed in
 * inittransitionmatrixnc. This computation is done using the eigendecomposition
 * of A, which is precomputed by inittransitionmatrixnc. The result is stored in
 * PMATnc[n]. Note that the matrix PMATnc[n] is of dimension 4x5; the last
 * column has all entries 1.0 (for missing data). */
void maketransitionmatrixnc_Arr(int n, double t) {
    double expos[4];
    for (int k = 0; k < 4; k++) expos[k] = exp(t * eigenvalues[k]);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            PMATnc[n][i][j] = 0.0;
            for (int k = 0; k < 4; k++)
                PMATnc[n][i][j] += V[i][k] * expos[k] * Vinv[k][j];
        }
    }
}

/* See comment for maketransitionmatrixnc_Arr */
void inittransitionmatrixnc(double pi[4], double par[8]) {
    gsl_matrix *Q = gsl_matrix_alloc(4, 4);

    // Q is the rate matrix of the continuous-time Markov chain that governs
    // nucleotide substitutions. It satisfies the property that e^{Qt} is the
    // transition matrix for nucleotide substitutions over a time t.
    gsl_matrix_set(Q, 0, 1, pi[1] * par[3]);
    gsl_matrix_set(Q, 0, 2, pi[2] * par[4]);
    gsl_matrix_set(Q, 0, 3, pi[3] * par[5]);
    gsl_matrix_set(Q, 1, 0, pi[0] * par[3]);
    gsl_matrix_set(Q, 1, 2, pi[2] * par[6]);
    gsl_matrix_set(Q, 1, 3, pi[3] * par[7]);
    gsl_matrix_set(Q, 2, 0, pi[0] * par[4]);
    gsl_matrix_set(Q, 2, 1, pi[1] * par[6]);
    gsl_matrix_set(Q, 2, 3, pi[3]); // unscaled rate of GT = 1.0
    gsl_matrix_set(Q, 3, 0, pi[0] * par[5]);
    gsl_matrix_set(Q, 3, 1, pi[1] * par[7]);
    gsl_matrix_set(Q, 3, 2, pi[2]); // unscaled rate of GT = 1.0
    for (int i = 0; i < 4; i++) {
        gsl_matrix_set(Q, i, i, 0.0);
        double sum = 0.0;
        for (int j = 0; j < 4; j++) sum = sum + gsl_matrix_get(Q, i, j);
        gsl_matrix_set(Q, i, i, -sum);
    }

    gsl_vector_view pi_gsl = gsl_vector_view_array(pi, 4);

    // compute eigendecomposition of Q
    gsl_vector *eig_gsl = gsl_vector_alloc(4);
    gsl_matrix *V_gsl = gsl_matrix_alloc(4, 4);
    gsl_matrix *Vinv_gsl = gsl_matrix_alloc(4, 4);
    eigendecompose_gtr(Q, &pi_gsl.vector, eig_gsl, V_gsl, Vinv_gsl);

    // store eigendecomposition in global variables V, eigenvalues, Vinv
    for (int i = 0; i < 4; i++) {
        eigenvalues[i] = gsl_vector_get(eig_gsl, i);
        for (int j = 0; j < 4; j++) {
            V[i][j] = gsl_matrix_get(V_gsl, i, j);
            Vinv[i][j] = gsl_matrix_get(Vinv_gsl, i, j);
        }
    }

    for (int i = 0; i < 4; i++) {
        PMATnc[0][i][4] = 1.0; // missing data
        PMATnc[1][i][4] = 1.0; // missing data
    }

    gsl_matrix_free(Q);
    gsl_matrix_free(V_gsl);
    gsl_matrix_free(Vinv_gsl);
    gsl_vector_free(eig_gsl);
}

void _makeconnc_Arr(int node, int whichRoot, double *UFCnc) {
    double lambda = 1.0; /* might use in future... */

    int child = treeArr[whichRoot][node].up[0];
    if (treeArr[whichRoot][child].up[0] == -1) {
        maketransitionmatrixnc_Arr(0, lambda * treeArr[whichRoot][child].bl);
        int seqn = child - numspecArr[whichRoot] + 1;
        for (int site = 0; site < numbaseArr[whichRoot]; site++) {
            for (int i = 0; i < 4; i++)
                treeArr[whichRoot][node].likenc[site][i] =
                    PMATnc[0][i][seqArr[whichRoot][seqn][site]];
        }
    } else {
        _makeconnc_Arr(child, whichRoot, UFCnc);
        maketransitionmatrixnc_Arr(0, lambda * treeArr[whichRoot][child].bl);
        for (int site = 0; site < numbaseArr[whichRoot]; site++) {
            for (int i = 0; i < 4; i++) {
                treeArr[whichRoot][node].likenc[site][i] = 0.0;
                for (int j = 0; j < 4; j++)
                    treeArr[whichRoot][node].likenc[site][i] +=
                        PMATnc[0][i][j] *
                        treeArr[whichRoot][child].likenc[site][j];
            }
        }
    }
    child = treeArr[whichRoot][node].up[1];
    if (treeArr[whichRoot][child].up[1] == -1) {
        maketransitionmatrixnc_Arr(0, lambda * treeArr[whichRoot][child].bl);
        int seqn = child - numspecArr[whichRoot] + 1;
        for (int site = 0; site < numbaseArr[whichRoot]; site++)
            for (int i = 0; i < 4; i++)
                treeArr[whichRoot][node].likenc[site][i] =
                    treeArr[whichRoot][node].likenc[site][i] *
                    PMATnc[0][i][seqArr[whichRoot][seqn][site]];
    } else {
        _makeconnc_Arr(child, whichRoot, UFCnc);
        maketransitionmatrixnc_Arr(0, lambda * treeArr[whichRoot][child].bl);
        for (int site = 0; site < numbaseArr[whichRoot]; site++) {
            double max = 0.0;
            for (int i = 0; i < 4; i++) {
                double L = 0.0;
                for (int j = 0; j < 4; j++)
                    L += PMATnc[0][i][j] *
                         treeArr[whichRoot][child].likenc[site][j];

                // underflow control
                if ((treeArr[whichRoot][node].likenc[site][i] =
                         treeArr[whichRoot][node].likenc[site][i] * L) > max)
                    max = treeArr[whichRoot][node].likenc[site][i];
            }

            if (max < 1e-300)
                printf("WARNING: max = %lg (%s:%d)\n", max, __FILE__, __LINE__);

            for (int i = 0; i < 4; i++)
                treeArr[whichRoot][node].likenc[site][i] =
                    treeArr[whichRoot][node].likenc[site][i] / max;
            UFCnc[site] = UFCnc[site] + log(max);
        }
    }
}

/* Compute likelihood of tree given by `whichRoot` under the following
 * parameters:
 *
 *      GTR model
 *        par[0]: inverse function for piA;
 *        par[1]: inverse function for piC;
 *        par[2]: inverse function for piG;
 *        par[3]: AC;
 *        par[4]: AG;
 *        par[5]: AT;
 *        par[6]: CG;
 *        par[7]: CT;
 *
 *        (the unscaled rate of G<->T is defined to be 1.0)
 **/
double getlike_Arr(double par[8], int whichRoot) {

    double stand = 1.0 + par[0] + par[1] + par[2];

    double pi[4];
    pi[0] = par[0] / stand;
    pi[1] = par[1] / stand;
    pi[2] = par[2] / stand;
    pi[3] = 1.0 - pi[0] - pi[1] - pi[2];

    double *UFCnc = malloc((numbaseArr[whichRoot]) * (sizeof(double)));
    double *locloglike = malloc(numbaseArr[whichRoot] * sizeof *locloglike);

    inittransitionmatrixnc(pi, par);

    for (int i = 0; i < numbaseArr[whichRoot]; i++) UFCnc[i] = 0.0;
    _makeconnc_Arr(rootArr[whichRoot], whichRoot, UFCnc);
    for (int i = 0; i < numbaseArr[whichRoot]; i++) {
        double L = 0.0;
        for (int k = 0; k < 4; k++)
            L += treeArr[whichRoot][rootArr[whichRoot]].likenc[i][k] * pi[k];
        if (L > 0.0) locloglike[i] = log(L) + UFCnc[i];
    }

    double like = 0.0;
    for (int i = 0; i < numbaseArr[whichRoot]; i++) {
        double loclike = 0.0;
        double max = -1e300;
        if (locloglike[i] > max) // underflow protection
            max = locloglike[i];
        double d = locloglike[i] - max;
        if (d > -100) loclike += exp(d);
        like = like + log(loclike) + max;
    }
    free(locloglike);
    free(UFCnc);
    return -like;
}

// OPTIMIZATION STUFF:
/* The global tree on which optimization is currently being performed */
static int _local_whichRoot = -2;

/* Set tree to be optimized over */
void _set_local_whichRoot(int new_local_whichRoot) {
    _local_whichRoot = new_local_whichRoot;
}

/* A wrapper function for the likelihood under a fixed tree (_local_whichRoot).
 * This function is used for optimization. */
double _getlike_Arr_wrapper(double parameters[8]) {
    double ret = -getlike_Arr(parameters, _local_whichRoot);
    if (!isfinite(ret)) {
        printf("Warning: non-finite likelihood encountered: L=%lg (%s:%d)\n",
               ret, __FILE__, __LINE__);
        ret = -1e300; // "negative infinity"
    }
    return ret;
}

/* Call wrapper function to maximize the likelihood under the 8 GTR parameters
 * described elsewhere. THIS FUNCTION ASSUMES _local_whichRoot IS ALREADY SET!!!
 */
double _maximizelikelihoodnc_globals_Arr(double parameters[8],
                                         double *step_size, int num_iter) {

    double params_min[8] = {0.005, 0.005, 0.005, 0.005, // reasonable minima
                            0.005, 0.005, 0.005, 0.005};

    findmax_nm(_getlike_Arr_wrapper, 8, parameters, params_min, step_size,
               num_iter);
    return _getlike_Arr_wrapper(parameters);
}

/* The pi and node which _like_bl_Arr is currently computing the likelihood
 * over. */
static double *_local_pi;
static int _local_node;

// Calculate log-likelihood of the branch length of `localnode` toward its
// parent. This branch length should be in *bl. The function maxbl_nc_Arr
// then finds the maximum likelihood branch length. THIS FUNCTION REQUIRES THAT
// _local_whichRoot, _localpi and _localnode ARE SET!!!!! That is the tree, pi
// and node on which it operates, respectively.
double _like_bl_Arr(double *bl) {
    int whichRoot = _local_whichRoot;

    int i, j, s, base;
    double b, p, L = 0.0;

    maketransitionmatrixnc_Arr(0, *bl);

    for (s = 0; s < numbaseArr[whichRoot]; s++) {
        p = 0.0;
        if (treeArr[whichRoot][_local_node].up[0] == -1) {
            base =
                seqArr[whichRoot][_local_node - numspecArr[whichRoot] + 1][s];
            assert(base >= 0 && base <= 4);
            if (base < 4) {
                for (j = 0; j < 4; j++) {
                    p += _local_pi[base] * PMATnc[0][base][j] *
                         templike_nc[s][j];
                }
            } else p = 1.0;
        } else {
            for (i = 0; i < 4; i++) {
                b = 0;
                for (j = 0; j < 4; j++)
                    b += PMATnc[0][i][j] * templike_nc[s][j];
                p += b * _local_pi[i] *
                     treeArr[whichRoot][_local_node].likenc[s][i];
            }
        }
        L += log(p);
    }


    if (!isfinite(L)) {
        printf("Warning: non-finite likelihood encountered: L=%lg p=%lg "
               "*bl=%lg (%s:%d)\n",
               L, p, *bl, __FILE__, __LINE__);
        L = -1e300; // "negative infinity"
    }
    return L;
}

// Find the branch length of node `node` toward its parent whose likelihood is
// maximal, and update that branch length in the tree.
void _maxbl_nc_Arr(int node, double pi[4], int whichRoot) {
    double bl = treeArr[whichRoot][node].bl, par_min = MINBL;

    _local_pi = pi;
    _local_node = node;

    // Other hyperparameters would also be possible.
    double step_size = 0.001;
    findmax_nm(_like_bl_Arr, 1, &bl, &par_min, &step_size, 10 /*num_iter*/);

    // printf("%lg => %lg\n", treeArr[whichRoot][node].bl, *bl);
    treeArr[whichRoot][node].bl = bl;
}

// Recursive helper function for the estimatebranchlengths_Arr() (see below).
void _recurse_estimatebranchlengths_Arr(int node, double pi[4], int whichRoot) {
    int i, j, s, parent, otherb, child1, child2;
    double max;
    child1 = treeArr[whichRoot][node].up[0];
    parent = treeArr[whichRoot][node].down;
    if ((otherb = treeArr[whichRoot][parent].up[0]) == node)
        otherb = treeArr[whichRoot][parent].up[1];
    maketransitionmatrixnc_Arr(1, treeArr[whichRoot][otherb].bl);
    for (s = 0; s < numbaseArr[whichRoot]; s++) {
        for (i = 0; i < 4; i++) {
            max = 0.0;
            if (treeArr[whichRoot][otherb].up[0] == -1) {
                templike_nc[s][i] =
                    PMATnc[1][i][seqArr[whichRoot]
                                       [otherb - numspecArr[whichRoot] + 1][s]];
            } else {
                templike_nc[s][i] = 0.0;
                for (j = 0; j < 4; j++)
                    templike_nc[s][i] =
                        templike_nc[s][i] +
                        treeArr[whichRoot][otherb].likenc[s][j] *
                            PMATnc[1][i][j];
            }
            if ((templike_nc[s][i] =
                     templike_nc[s][i] *
                     treeArr[whichRoot][parent].posteriornc[s][i]) > max)
                max = templike_nc[s][i];
        }
        // for (i=0; i<4; i++)//underflow protection
        // templike_nc[s][i]=templike_nc[s][i]/max;
    }

    _maxbl_nc_Arr(node, pi, whichRoot);

    maketransitionmatrixnc_Arr(0, treeArr[whichRoot][node].bl);
    for (s = 0; s < numbaseArr[whichRoot]; s++) {
        max = 0.0;
        for (i = 0; i < 4; i++) {
            treeArr[whichRoot][node].posteriornc[s][i] = 0.0;
            for (j = 0; j < 4; j++)
                treeArr[whichRoot][node].posteriornc[s][i] =
                    treeArr[whichRoot][node].posteriornc[s][i] +
                    PMATnc[0][i][j] * templike_nc[s][j];
            if (treeArr[whichRoot][node].posteriornc[s][i] >
                max) // more hysterical underflow protection
                max = treeArr[whichRoot][node].posteriornc[s][i];
        }
        for (i = 0; i < 4; i++)
            treeArr[whichRoot][node].posteriornc[s][i] =
                treeArr[whichRoot][node].posteriornc[s][i] / max;
    }
    if (child1 > -1) {
        child2 = treeArr[whichRoot][node].up[1];
        _recurse_estimatebranchlengths_Arr(child1, pi, whichRoot);
        _recurse_estimatebranchlengths_Arr(child2, pi, whichRoot);
    }
}

/* This function does one pass at estimating the ML branch lengths. The new
 * branch lengths (with higher likelihood) are written to the tree. */
void estimatebranchlengths_Arr(double par[8], int whichRoot) {
    int i, j, s, child1, child2;
    double stand, pi[4];
    stand = 1.0 + par[0] + par[1] + par[2];
    pi[0] = par[0] / stand;
    pi[1] = par[1] / stand;
    pi[2] = par[2] / stand;
    pi[3] = 1.0 - pi[0] - pi[1] - pi[2];
    child1 = treeArr[whichRoot][rootArr[whichRoot]].up[0];
    child2 = treeArr[whichRoot][rootArr[whichRoot]].up[1];
    templike_nc = malloc(numbaseArr[whichRoot] * (sizeof(double *)));
    for (i = 0; i < numbaseArr[whichRoot]; i++)
        templike_nc[i] = malloc(4 * (sizeof(double)));

    /*Root is different because only one parameter can be estimated for the two
     * child nodes. We arbitrarily fix one of the branchlengths to equal MINBL
     * and only estimate the other one*/

    if ((treeArr[whichRoot][child2].bl = treeArr[whichRoot][child2].bl +
                                         treeArr[whichRoot][child1].bl -
                                         MINBL) < MINBL)
        treeArr[whichRoot][child2].bl = MINBL;
    treeArr[whichRoot][child1].bl = MINBL;
    getlike_Arr(par, whichRoot); /*this is not necessary if the likelihood
                                          fucntion has already been called*/
    maketransitionmatrixnc_Arr(0, treeArr[whichRoot][child2].bl +
                                      treeArr[whichRoot][child1].bl);
    for (s = 0; s < numbaseArr[whichRoot]; s++) {
        for (i = 0; i < 4; i++) {
            treeArr[whichRoot][rootArr[whichRoot]].posteriornc[s][i] =
                1.0; /*tree[child1].likenc[s][i];*/
            if (treeArr[whichRoot][child2].up[0] > -1) {
                treeArr[whichRoot][child1].posteriornc[s][i] = 0.0;
                for (j = 0; j < 4; j++)
                    treeArr[whichRoot][child1].posteriornc[s][i] +=
                        treeArr[whichRoot][child2].likenc[s][j] *
                        PMATnc[0][i][j];
            } else
                treeArr[whichRoot][child1].posteriornc[s][i] =
                    PMATnc[0][i][seqArr[whichRoot]
                                       [child2 - numspecArr[whichRoot] + 1][s]];
        }
    }
    if (treeArr[whichRoot][child1].up[0] > -1)
        _recurse_estimatebranchlengths_Arr(treeArr[whichRoot][child1].up[0], pi,
                                           whichRoot);
    if (treeArr[whichRoot][child1].up[1] > -1)
        _recurse_estimatebranchlengths_Arr(treeArr[whichRoot][child1].up[1], pi,
                                           whichRoot);
    _recurse_estimatebranchlengths_Arr(child2, pi, whichRoot);
    for (i = 0; i < numbaseArr[whichRoot]; i++) free(templike_nc[i]);
    free(templike_nc);
}

// Given a tree, this function optimizes two things at the same time (by doing
// multiple passes): THE BRANCH LENGTHS, and the following 8 global parameters:
//
//    par[0]: inverse function for piA;
//    par[1]: inverse function for piC;
//    par[2]: inverse function for piG;
//    par[3]: AC;
//    par[4]: AG;
//    par[5]: AT;
//    par[6]: CG;
//    par[7]: CT;
//
//    (the unscaled rate of G<->T is defined to be 1.0)
//
// The input array `parameters` should be initialized with initial parameters.
// When this function returns, the argument will hold the updated parameters.
void estimatenucparameters_Arr(double parameters[8], int whichRoot) {

    clearGlobals();

    _set_local_whichRoot(whichRoot); // very important!!!!!

    double L = _getlike_Arr_wrapper(parameters);
    printf("============================================\nTREE %d\n"
           "Initial ML value = %lf (with params",
           whichRoot, L);
    for (int j = 0; j < 8; j++) printf(" %.2lf", parameters[j]);
    printf(")\n");

    // step sizes will be automatically adjusted below
    double step_size[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    // XXX: Other hyperparameters are also possible, but these work well.
    for (int i = 0; i < 15; i++) {
        double oldparameters[8];
        memcpy(oldparameters, parameters, sizeof oldparameters);

        estimatebranchlengths_Arr(parameters, whichRoot);

        L = _maximizelikelihoodnc_globals_Arr(parameters, step_size, 6);

        // some educated guessing to try to find initial step sizes for the
        // next iteration that match the change in that parameter during the
        // current iteration. NM works in logspace so we use d(log x) = x/dx.
        for (int j = 0; j < 8; j++)
            step_size[j] = 0.6 * step_size[j] + 0.001 +
                           0.9 * fabs((parameters[j] - oldparameters[j]) /
                                      (parameters[j]));

        printf("Current ML value = %lf (with params", L);
        for (int j = 0; j < 8; j++) printf(" %.2lf", parameters[j]);
        printf(")\n");
    }
}

void clearGlobals(void) {
    int i, j;
    for (i = 0; i < STATESPACE; i++) {
        RRVAL[i] = 0;
        for (j = 0; j < STATESPACE; j++) {
            LRVEC[i][j] = 0;
            RRVEC[i][j] = 0;
            PMAT1[i][j] = 0;
            PMAT2[i][j] = 0;
        }
    }
    for (i = 4; i < 4; i++) {
        eigenvalues[i] = 0;
        for (i = 4; i < 4; i++) {
            V[i][j] = 0;
            Vinv[i][j] = 0;
        }
    }
    for (i = 0; i < 4; i++) {
        PMATnc[0][i][4] = 1.0; // missing data
        PMATnc[1][i][4] = 1.0; // missing data
    }
}

void makeposterior_nc_Arr(int node, int whichRoot) {
    int i, j, s, parent, otherb, child1, child2, b;
    double bl, max;


    child1 = treeArr[whichRoot][node].up[0];
    child2 = treeArr[whichRoot][node].up[1];
    parent = treeArr[whichRoot][node].down;
    bl = treeArr[whichRoot][node].bl;
    maketransitionmatrixnc_Arr(0, bl);
    if ((otherb = treeArr[whichRoot][parent].up[0]) == node)
        otherb = treeArr[whichRoot][parent].up[1];
    maketransitionmatrixnc_Arr(1, treeArr[whichRoot][otherb].bl);
    for (s = 0; s < numbaseArr[whichRoot]; s++) {
        if (treeArr[whichRoot][otherb].up[0] > -1) {
            for (i = 0; i < 4; i++) {
                templike_nc[s][i] = 0;
                for (j = 0; j < 4; j++)
                    templike_nc[s][i] +=
                        treeArr[whichRoot][otherb].likenc[s][j] *
                        PMATnc[1][i][j];
                templike_nc[s][i] =
                    templike_nc[s][i] *
                    treeArr[whichRoot][parent].posteriornc[s][i];
            }
        } else { // if leaf node (check otherb leaf node)
            b = seqArr[whichRoot][otherb - numspecArr[whichRoot] + 1][s];
            for (i = 0; i < 4; i++)
                templike_nc[s][i] =
                    PMATnc[1][i][b] *
                    treeArr[whichRoot][parent].posteriornc[s][i];
        }
        for (i = 0; i < 4; i++) {
            treeArr[whichRoot][node].posteriornc[s][i] = 0.0;
            max = 0.0;
            for (j = 0; j < 4; j++) {
                if ((treeArr[whichRoot][node].posteriornc[s][i] =
                         treeArr[whichRoot][node].posteriornc[s][i] +
                         PMATnc[0][i][j] * templike_nc[s][j]) > max)
                    max = treeArr[whichRoot][node]
                              .posteriornc[s][i]; // more underflow protection
            }
        }
        for (i = 0; i < 4; i++)
            treeArr[whichRoot][node].posteriornc[s][i] =
                treeArr[whichRoot][node].posteriornc[s][i] / max;
    }
    if (treeArr[whichRoot][child1].up[0] > -1)
        makeposterior_nc_Arr(child1, whichRoot);
    if (treeArr[whichRoot][child2].up[0] > -1)
        makeposterior_nc_Arr(child2, whichRoot);
}

void set_posteriors(int whichRoot) {
    int i, j, k;
    for (i = 0; i < 2 * numspecArr[whichRoot] - 1; i++) {
        for (j = 0; j < numbaseArr[whichRoot]; j++) {
            for (k = 0; k < 4; k++) {
                treeArr[whichRoot][i].posteriornc[j][k] =
                    log(0.01 / 4 +
                        (1 - 0.01) * treeArr[whichRoot][i].posteriornc[j][k]);
            }
        }
    }
}

// THIS FUNCTION IS IMPLEMENTED WITHOUT GAMMA CATEGORIES
void getposterior_nc_Arr(double par[8], int whichRoot) {
    int i, j, s, k, parent, b;
    double sum, pi[4], stand;
    getlike_Arr(par, whichRoot); /*need to call likelihood again to get
                                    .likenc[.][.] for the current parameters*/
    templike_nc = malloc(numbaseArr[whichRoot] * (sizeof(double *)));
    for (i = 0; i < numbaseArr[whichRoot]; i++)
        templike_nc[i] = malloc(4 * (sizeof(double)));
    stand = 1.0 + par[0] + par[1] + par[2];
    pi[0] = par[0] / stand;
    pi[1] = par[1] / stand;
    pi[2] = par[2] / stand;
    pi[3] = 1.0 - pi[0] - pi[1] - pi[2];
    for (s = 0; s < numbaseArr[whichRoot]; s++)
        for (i = 0; i < 4; i++)
            treeArr[whichRoot][rootArr[whichRoot]].posteriornc[s][i] = 1.0;
    if (treeArr[whichRoot][treeArr[whichRoot][rootArr[whichRoot]].up[0]].up[0] >
        -1)
        makeposterior_nc_Arr(treeArr[whichRoot][rootArr[whichRoot]].up[0],
                             whichRoot);
    if (treeArr[whichRoot][treeArr[whichRoot][rootArr[whichRoot]].up[1]].up[0] >
        -1)
        makeposterior_nc_Arr(treeArr[whichRoot][rootArr[whichRoot]].up[1],
                             whichRoot);
    for (j = 0; j < 2 * numspecArr[whichRoot] - 1; j++) {
        if (treeArr[whichRoot][j].up[0] > -1) {
            for (s = 0; s < numbaseArr[whichRoot]; s++) {
                sum = 0.0;
                for (i = 0; i < 4; i++)
                    sum = sum +
                          (treeArr[whichRoot][j].posteriornc[s][i] =
                               treeArr[whichRoot][j].likenc[s][i] *
                               treeArr[whichRoot][j].posteriornc[s][i] * pi[i]);
                for (i = 0; i < 4; i++)
                    treeArr[whichRoot][j].posteriornc[s][i] =
                        treeArr[whichRoot][j].posteriornc[s][i] / sum;
            }
        } else {
            for (s = 0; s < numbaseArr[whichRoot]; s++) {
                b = seqArr[whichRoot][j - numspecArr[whichRoot] + 1][s];
                if (b == 4) {
                    maketransitionmatrixnc_Arr(0, treeArr[whichRoot][j].bl);
                    parent = treeArr[whichRoot][j].down;
                    sum = 0.0;
                    for (i = 0; i < 4; i++) {
                        treeArr[whichRoot][j].posteriornc[s][i] = 0.0;
                        for (k = 0; k < 4; k++)
                            treeArr[whichRoot][j].posteriornc[s][i] +=
                                treeArr[whichRoot][parent].posteriornc[s][k] *
                                PMATnc[0][i][k];
                        sum =
                            sum + (treeArr[whichRoot][j].posteriornc[s][i] =
                                       treeArr[whichRoot][j].posteriornc[s][i] *
                                       pi[i]);
                    }
                    for (i = 0; i < 4; i++)
                        treeArr[whichRoot][j].posteriornc[s][i] =
                            treeArr[whichRoot][j].posteriornc[s][i] / sum;
                } else {
                    for (i = 0; i < 4; i++) {
                        if (i == b)
                            treeArr[whichRoot][j].posteriornc[s][i] = 1.0;
                        else treeArr[whichRoot][j].posteriornc[s][i] = 0.0;
                    }
                }
            }
        }
    }
    for (i = 0; i < numbaseArr[whichRoot]; i++) free(templike_nc[i]);
    free(templike_nc);
}
