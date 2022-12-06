#ifndef OPT_H
#define OPT_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define NMAX 10000000
#define NR_END 1
#define FREE
#define FREE_ARG char*
#define ALF 1.0e-4
#define EPS 3.0e-8

extern double ITMAX;
extern double TOLX;
extern double TOLX2;
extern double STPMX;
extern double *dg,*g,*hdg,*pnew,*xi,**hessin; 
extern int npar, CENTRALMODE, PRECISIONLEVEL;
extern FILE *tempfile;
extern double oldf0;

void nrerror(char error_text[]);
double *dvector(long nl, long nh);
void free_dvector(double *v, long nl);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long ncl);
//int lnsrch(int n, double xold[], double fold, double g[], double p[], double x[], double *f, double stpmax, int *check, double (*func)(double []), double lowbound[], double upbound[]);
int lnsrch_Arr(int n, double xold[], double fold, double g[], double p[], double x[], double *f, double stpmax, int *check, double (*func)(double [], int), double lowbound[], double upbound[], int whichRoot);
void doNRinits(int n);
void freeNRinits(int n);
//void dfpmin(double p[], int n, double gtol, int *iter, double *fret, double(*func)(double []), void (*dfunc)(double [], double [],double [], double [], double(*fu)(double [])), double lowbound[], double upbound[]);
void dfpmin_Arr(double p[], int n, double gtol, int *iter, double *fret, double(*func)(double [], int), void (*dfunc)(double [], double [],double [], double [], double(*fu)(double [], int), int), double lowbound[], double upbound[], int whichRoot);
//void Yanggradient (int n, double x[], double f0, double g[], double (*fun)(double x[]), double space[], int central, double lowbound[], double upbound[]);
void Yanggradient_Arr (int n, double x[], double f0, double g[], double (*fun)(double x[], int), double space[], int central, double lowbound[], double upbound[], int whichRoot);
//void getgradient(double invec[], double outvec[], double lowbound[], double upbound[], double(*func)(double []));
void getgradient_Arr(double invec[], double outvec[], double lowbound[], double upbound[], double(*func)(double [], int), int whichRoot);
//double findmax(double newinvecter[], double lowbound[], double upbound[], int n, double (*fun)(double x[]), int precisionlevel);
double findmax_Arr(double newinvecter[], double lowbound[], double upbound[], int n, double (*fun)(double x[], int), int precisionlevel, int whichRoot);

#endif /* OPT_H */
