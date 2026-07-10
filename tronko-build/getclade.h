/*
 *	getclade.h
 */
#ifndef _GET_CLADE_
#define _GET_CLADE_

#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include "opt.h"
#include "global.h"

int getcladeArr(FILE *tre, struct masterArr *m, int max_nodename);
void linknodesArr(int i, int j, int node, struct masterArr *m);
int getnodenumbArr(FILE *tre);
int specsearchArr(FILE *tre, struct masterArr *m, int max_node_name);

#endif /* _GET_CLADE_ */
