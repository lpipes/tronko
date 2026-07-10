#ifndef _ALLOC_TREE_
#define _ALLOC_TREE_

#include <stdio.h>
#include <stdlib.h>
// #include <math.h>
#include <ctype.h>
#include "global.h"
#include "opt.h"

void allocatetreememory_for_nucleotide_Arr(int numberOfTrees);
void allocateTreeArrMemory(struct masterArr *m, int max_nodename);
void freeTreeMemory(int whichPartition);
void *calloc_check(size_t nmemb, size_t size);
#endif /* _ALLOC_TREE */
