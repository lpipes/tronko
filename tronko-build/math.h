#ifndef _MATH_H
#define _MATH_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "tools.h"
#include "opt.h"
#include "global.h"

#define PointGamma(prob,alpha,beta) PointChi2(prob,2.0*(alpha))/(2.0*(beta))
#define BASE        2    /* base of floating point arithmetic */
#define DIGITS     40    /* no. of digits to the base BASE in the fraction */
#define MAXITER    30    /* max. no. of iterations to converge */
#define pos(i,j,n)      ((i)*(n)+(j))
#define csize(a) (fabs(a.re)+fabs(a.im))

double LnGamma (double alpha);
double IncompleteGamma (double x, double alpha, double ln_gamma_alpha);
double PointNormal (double prob);
double PointChi2 (double prob, double v);
double CDFfunGamma(double x, double par[2]);
void definegammaquantiles(int k, double par[2]);
void initlogfactorial();
int eigen(int job, double A[], int n, double rr[], double ri[], double vr[], double vi[], double w[]);
void balance(double mat[], int n, int *low, int *hi, double scale[]);
void unbalance(int n, double vr[], double vi[], int low, int hi, double scale[]);
int realeig(int job, double mat[], int n,int low, int hi, double valr[], double vali[], double vr[], double vi[]);
void elemhess(int job, double mat[], int n, int low, int hi, double vr[], double vi[], int work[]);
int eigen(int job, double A[], int n, double rr[], double ri[],double vr[], double vi[], double work[]);
complex compl (double re,double im);
complex conj (complex a);
complex cplus (complex a, complex b);
complex cminus (complex a, complex b);
complex cby (complex a, complex b);
complex cdiv (complex a,complex b);
complex cexp (complex a);
complex cfactor (complex x, double a);
int cxtoy (complex x[], complex y[], int n);
int cmatby (complex a[], complex b[], complex c[], int n,int m,int k);
int cmatout (FILE * fout, complex x[], int n, int m);
int cmatinv( complex x[], int n, int m, double space[]);
void balance(double mat[], int n,int *low, int *hi, double scale[]);
void unbalance(int n,double vr[],double vi[], int low, int hi, double scale[]);
void elemhess(int job,double mat[],int n,int low,int hi, double vr[], double vi[], int work[]);
int realeig(int job,double mat[],int n,int low, int hi, double valr[], double vali[], double vr[],double vi[]);
int matinv( double x[], int n, int m, double space[]);
void inittransitionmatrix(double lambda);
void maketransitionmatrix(int matnum, double t);
double pdata(double theta, int s, int nt, int nb);

#endif /* _MATH_H_ */
