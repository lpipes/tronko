#ifndef _ASSIGNMEMFORTHREADS_
#define _ASSIGNMEMFORTHREADS_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global.h"

type_of_PP **assignScores_Arr(int rootNum, int node, char *locQuery, int *positions, type_of_PP **scores, int alength);
void assignScores_Arr_paired( int rootNum, int node, char *locQuery, int *positions, type_of_PP ***scores, int alength, int search_number, int print_all_nodes, FILE* site_scores_file);
type_of_PP getscore_Arr(int alength, int node, int rootNum, char *locQuery, int *positions, int print_all_nodes, FILE* site_scores_file);
#endif /* _ASSIGNMEMFORTHREADS_ */
