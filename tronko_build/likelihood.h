#ifndef _LIKELIHOOD_
#define _LIKELIHOOD_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "tools.h"
#include "opt.h"
#include "global.h"

//void makecon(int node);
void makecon_Arr(int node, int whichRoot);
//double tlike();
double tlike_Arr(int whichRoot);
//double getlikelihood(double par[4]);
double getlikelihood_Arr(double par[4], int whichRoot);
//double maximizelikelihood(double parameters[4], int precision);
double maximizelikelihood_Arr(double parameters[4], int precision, int whichRoot);
void makeposterior(int node);
void getposterior(double par[4]);
//void maketransitionmatrixnc(int n, double t);
void maketransitionmatrixnc_Arr(int n, double t, int whichRoot);
//void makeconnc(int node, double lambda);
void makeconnc_Arr(int node, double lambda, int whichRoot);
void inittransitionmatrixnc(double pi[4], double par[10]);
//double getlike_gamma(double par[10]);
double getlike_gamma_Arr(double par[10], int whichRoot);
//double maximizelikelihoodnc_globals(double parameters[10], int precision);
double maximizelikelihoodnc_globals_Arr(double parameters[10], int precision, int whichRoot);
//double like_bl(double par[2]);
double like_bl_Arr(double par[2], int whichRoot);
//void maxbl_nc(int node, int parent, double pi[4], int precision);
void maxbl_nc_Arr(int node, int parent, double pi[4], int precision, int whichRoot);
//void recurse_estimatebracnhlengths(int node, double pi[4], int precision);
void recurse_estimatebracnhlengths_Arr(int node, double pi[4], int precision, int whichRoot);
//void estimatebracnhlengths(double par[10], int precision);
void estimatebracnhlengths_Arr(double par[10], int precision, int whichRoot);
//void estimatenucparameters(double parameters[10]);
void estimatenucparameters_Arr(double parameters[10], int whichRoot);
void clearGlobals();
//void makeposterior_nc(int node);
void makeposterior_nc_Arr(int node, int whichRoot);
//void getposterior_nc(double par[10]);
void getposterior_nc_Arr(double par[10], int whichRoot);

#endif
