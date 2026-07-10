#ifndef _LIKELIHOOD_
#define _LIKELIHOOD_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "opt.h"
#include "global.h"
#include "math.h"

void maketransitionmatrixnc_Arr(int n, double t);
void inittransitionmatrixnc(double pi[4], double par[8]);
double getlike_Arr(double par[8], int whichRoot);

void estimatenucparameters_Arr(double parameters[8], int whichRoot);

void clearGlobals(void);
void makeposterior_nc_Arr(int node, int whichRoot);
void getposterior_nc_Arr(double par[8], int whichRoot);

#endif
