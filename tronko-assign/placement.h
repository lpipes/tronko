#ifndef _PLACEMENT_
#define _PLACEMENT_

#include <stdio.h>
#include "global.h"
#include <unistd.h>
#include <time.h>
#include "WFA2/wavefront_align.h"

void place_paired_withBLAST( char *query_1, char *query_2, char **rootSeqs, int numberOfTotalRoots, int *positions, char *locQuery, type_of_PP ***nodeScores, int **voteRoot, int number_of_matches , int **leaf_coordinates, int paired, type_of_PP* minimum_score, char *alignments_dir, char *forward_name, char *reverse_name, int print_alignments, char *leaf_sequence, int *positionsInRoot, int maxNumSpec, affine_penalties_t affine_penalties, int* starts_forward, char** cigars_forward, int* starts_reverse, char** cigars_reverse , int print_alignments_to_file, int use_leaf_portion, int padding);
#endif /* _PLACEMENT_ */
