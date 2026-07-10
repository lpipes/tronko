#include "opt.h"
#include <gsl/gsl_multimin.h>

typedef struct {
    double (*f)(double *);
    int n;
    double *offset;
} Ctx;

static double wrapper(const gsl_vector *x, void *data) {
    Ctx *ctx = data;
    double p[ctx->n];
    for (int i = 0; i < ctx->n; i++)
        p[i] = ctx->offset[i] + exp(gsl_vector_get(x, i));
    return -ctx->f(p);
}

/*
 * findmax_nm: Maximize f(params) over parameters bounded below by params_min
 * using GSL's Nelder-Mead simplex. Parameters are reparameterized as
 * p = params_min[i] + exp(x) so the unconstrained search space maps to
 * (params_min[i], +inf).
 *
 * params[] holds initial values on input and optimal values on output.
 * Stops after num_iter iterations; step_size is an array of length n containing
 * initial step sizes.
 */
void findmax_nm(double (*f)(double *), int n, double *params,
                double *params_min, double *step_size, int num_iter) {
    Ctx ctx = {f, n, params_min};
    gsl_multimin_function fn = {wrapper, n, &ctx};

    gsl_vector *x = gsl_vector_alloc(n);
    gsl_vector *step = gsl_vector_alloc(n);
    for (int i = 0; i < n; i++) {
        gsl_vector_set(x, i,
                       params[i] > params_min[i]
                           ? log(params[i] - params_min[i])
                           : -1e300);
        gsl_vector_set(step, i, step_size[i]);
    }

    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T, n);
    gsl_multimin_fminimizer_set(s, &fn, x, step);

    for (int iter = 0; iter < num_iter; iter++) {
        if (gsl_multimin_fminimizer_iterate(s) != GSL_SUCCESS) break;
    }

    for (int i = 0; i < n; i++)
        params[i] = params_min[i] + exp(gsl_vector_get(s->x, i));

    gsl_multimin_fminimizer_free(s);
    gsl_vector_free(x);
    gsl_vector_free(step);
}
