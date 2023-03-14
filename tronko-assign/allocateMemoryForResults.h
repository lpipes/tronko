#ifndef _ALLOCFORPTHREADS_
#define _ALLOCFORPTHREADS_

#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "needleman_wunsch.h"
#include "alignment.h"
#include "alignment_scoring.h"

void allocateMemForResults( resultsStruct *results, int sizeOfChunk, int num_threads, int numberOfTrees, int print_alignments, int maxNumSpec, int paired, int use_nw, int max_lineTaxonomy, int max_name_length, int max_query_length, int max_numbase, int use_portion, int padding_size, int number_of_total_nodes);
void freeMemForResults ( resultsStruct *results, int sizeOfChunk, int numberOfTrees, int num_threads, int paired, int use_nw, int use_portion, int maxNumSpec, int number_of_total_nodes);
#endif /* _ALLOCFORPTHREADS_ */
