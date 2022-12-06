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
#include "tools.h"

void linknodes(int i,int j,int node);
int specsearch();
int getnodenumb();
int getclade();
int getcladeArr(FILE *tre, struct masterArr *m, int max_nodename);
int getcladeArr_UsePartitions(FILE *tre, int whichPartitions, char*** nodeIDsArr_heap);
void linknodesArr(int i, int j, int node, struct masterArr *m);
int getnodenumbArr(FILE *tre);
int specsearchArr(FILE *tre, struct masterArr *m, int max_node_name);
int specsearchArr_UsePartitions(FILE *tre, int whichPartitions, char*** nodeIDsArr_heap);
#endif /* _GET_CLADE_ */
