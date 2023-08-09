#ifndef _PLACEMENT_
#define _PLACEMENT_

#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include "global.h"
#include "WFA2/wavefront_align.h"
#include "printAlignments.h"
#include "assignment.h"

int perform_WFA_alignment(cigar_t* const cigar, mm_allocator_t* mm_allocator,char* seq1, char* seq2,char* const pattern_alg,char* const text_alg, char* const ops_alg, int begin_offset, int end_offset);
void place_paired( char *query_1, char *query_2, char **rootSeqs, int numberOfTotalRoots, int *positions, char *locQuery, type_of_PP ***nodeScores, int **voteRoot, int number_of_matches , int **leaf_coordinates, int paired, type_of_PP* minimum_score, char *alignments_dir, char *forward_name, char *reverse_name, int print_alignments, char *leaf_sequence, int *positionsInRoot, int maxNumSpec, int* starts_forward, char** cigars_forward, int* starts_reverse, char** cigars_reverse, int print_alignments_to_file, int use_leaf_portion, int padding, int max_query_length, int max_numbase, int print_all_nodes);
void place_paired_with_nw( char *query_1, char *query_2, char **rootSeqs, int numberOfTotalRoots, int *positions, char *locQuery, nw_aligner_t *nw, alignment_t *aln, scoring_t *scoring, type_of_PP ***nodeScores, int **voteRoot, int number_of_matches , int **leaf_coordinates, int paired, type_of_PP* minimum_score, char *alignments_dir, char *forward_name, char *reverse_name, int print_alignments, char *leaf_sequence, int *positionsInRoot, int maxNumSpec, int* starts_forward, char** cigars_forward, int* starts_reverse, char** cigars_reverse, int print_alignments_to_file, int use_leaf_portion, int padding, int max_query_length, int max_numbase, int print_all_nodes);
#endif /* _PLACEMENT_ */
