#ifndef _ALLOC_TREE_
#define _ALLOC_TREE_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "global.h"

void allocatetreememory_for_nucleotide_Arr(int numberOfTrees);
//void printNewick(FILE* file, int whichRoot, int idx);
void allocateTreeArrMemory(int whichPartition, int max_nodename);
void *calloc_check(size_t nmemb, size_t size);
void allocateMemoryForTaxArr(int whichPartitions, int max_tax_name_len);
void getReverseComplement(char *read, char *reverseComplement, int max_query_length);
#endif /* _ALLOC_TREE */
