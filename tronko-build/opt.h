#ifndef OPT_H
#define OPT_H

/* See file opt.c for documentation. */
void findmax_nm(double (*f)(double *), int n, double *params,
                double *params_min, double *step_size, int max_iter);

#endif /* OPT_H */
