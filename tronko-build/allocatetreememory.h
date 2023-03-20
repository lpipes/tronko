#ifndef _ALLOC_TREE_
#define _ALLOC_TREE_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "tools.h"
#include "global.h"
#include "opt.h"

void allocatetreememmory();
void allocatetreememmory_for_nucleotide();
void allocatetreememory_for_nucleotide_Arr(int numberOfTrees);
void allocateTreeArrMemory(struct masterArr *m, int max_nodename);
void freeTreeMemory(int whichPartition);
void *calloc_check(size_t nmemb, size_t size);
void allocateMemoryForTaxArr(int whichPartitions, int max_tax_name);
#endif /* _ALLOC_TREE */
