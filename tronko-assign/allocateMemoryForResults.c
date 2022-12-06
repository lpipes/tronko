#include "allocateMemoryForResults.h"

void allocateMemForResults( resultsStruct *results, int sizeOfChunk, int num_threads, int numberOfTrees, int print_alignments, int maxNumSpec, int paired, int use_nw, int max_lineTaxonomy, int max_name_length, int max_query_length, int max_numbase, int use_portion, int padding_size){
	int i,j, k;
	if (use_portion==1){
		results->positions = malloc((max_query_length+max_query_length+2*padding_size+1)*(sizeof(int)));
		results->locQuery = malloc((max_query_length+max_query_length+2*padding_size+1)*(sizeof(char)));
	}else{
		results->positions = malloc((max_query_length+max_numbase+1)*(sizeof(int)));
		results->locQuery = malloc((max_query_length+max_numbase+1)*(sizeof(char)));
	}
	results->nodeScores = (type_of_PP ***)malloc(MAX_NUM_BWA_MATCHES*(sizeof(type_of_PP **)));
	for (i=0; i<MAX_NUM_BWA_MATCHES; i++){
		results->nodeScores[i] = (type_of_PP **)malloc(numberOfTrees*(sizeof(type_of_PP *)));
		for (j=0; j<numberOfTrees; j++){
			results->nodeScores[i][j] = (type_of_PP *)malloc((2*numspecArr[j]-1)*(sizeof(type_of_PP)));
			for(k=0; k<2*numspecArr[j]-1; k++){
				results->nodeScores[i][j][k] = 0;
			}
		}
	}
	results->voteRoot = (int **)malloc(numberOfTrees*sizeof(int *));
	for (i=0; i<numberOfTrees; i++){
		results->voteRoot[i]=(int *)malloc((2*numspecArr[i]-1)*sizeof(int));
		for (j=0; j<2*numspecArr[i]-1; j++){
			results->voteRoot[i][j]=0;
		}
	}
	if ( use_portion==1){
		results->starts_forward = (int *)malloc(MAX_NUM_BWA_MATCHES*sizeof(int));
		if (paired == 1 ){
			results->starts_reverse = (int *)malloc(MAX_NUM_BWA_MATCHES*sizeof(int));
		}
		results->cigars_forward = (char **)malloc(MAX_NUM_BWA_MATCHES*sizeof(char *));
		if (paired == 1){
			results->cigars_reverse = (char **)malloc(MAX_NUM_BWA_MATCHES*sizeof(char *));
		}
		for(i=0; i<MAX_NUM_BWA_MATCHES; i++){
			results->starts_forward[i] = -1;
			results->cigars_forward[i] = (char *)malloc(MAX_CIGAR*sizeof(char));
			memset(results->cigars_forward[i],'\0',MAX_CIGAR);
			if ( paired==1){
				results->starts_reverse[i] = -1;
				results->cigars_reverse[i] = (char *)malloc(MAX_CIGAR*sizeof(char));
				memset(results->cigars_reverse[i],'\0',MAX_CIGAR);
			}
		}
	}
	if ( use_nw == 1 ){
		results->nw=needleman_wunsch_new();
		results->aln=alignment_create(max_query_length+max_numbase+1);
		results->scoring=malloc(sizeof(scoring_t));
	}
	results->minimum=(type_of_PP *)malloc(3*sizeof(type_of_PP));
	if ( print_alignments == 1){
		results->print_alignments=1;
	}else{
		results->print_alignments=0;
	}
	if (use_nw==1){
		int match=2;
		int mismatch=-1;
		int gap_open=-3;
		int gap_extend=-1;
		int num_of_mismatches=0;
		int num_of_indels = 0;
		bool no_start_gap_penalty=true;
		bool no_end_gap_penalty=true;
		bool no_gaps_in_a = false;
		bool no_gaps_in_b = false;
		bool no_mismatches = false;
		bool case_sensitive=false;
		scoring_init(results->scoring, match, mismatch, gap_open, gap_extend, no_start_gap_penalty, no_end_gap_penalty, no_gaps_in_a, no_gaps_in_b, no_mismatches, case_sensitive);
	}
	results->minNodes = (int *)malloc((2*maxNumSpec-1)*sizeof(int));
	results->LCAnames = (char **)malloc((2*maxNumSpec-1)*sizeof(char *));
	for(i=0;i<2*maxNumSpec-1;i++){
		results->LCAnames[i]=(char *)malloc(max_lineTaxonomy*sizeof(char));
	}
	results->leaf_coordinates = (int **)malloc(numberOfTrees*sizeof(int *));
	for(i=0; i<numberOfTrees; i++){
		results->leaf_coordinates[i] = (int *)malloc(2*sizeof(int));
		results->leaf_coordinates[i][0]=-1;
		results->leaf_coordinates[i][1]=-1;
	}
}
void freeMemForResults ( resultsStruct *results, int sizeOfChunk, int num_threads, int numberOfTrees, int paired, int use_nw, int use_portion, int maxNumSpec){
	int i, j, k;
	free(results->positions);
	free(results->locQuery);
	for(i=0; i<MAX_NUM_BWA_MATCHES; i++){
		for(j=0; j<numberOfTrees; j++){
			free(results->nodeScores[i][j]);
		}
		free(results->nodeScores[i]);
	}
	free(results->nodeScores);
	for(i=0; i<numberOfTrees; i++){
		free(results->voteRoot[i]);
		free(results->leaf_coordinates[i]);
	}
	free(results->voteRoot);
	for(i=0; i<2*maxNumSpec-1; i++){
		free(results->LCAnames[i]);
	}
	if (use_portion == 1){
		for(i=0; i<MAX_NUM_BWA_MATCHES; i++){
			free(results->cigars_forward[i]);
			if (paired==1){
				free(results->cigars_reverse[i]);
			}
		}
		free(results->starts_forward);
		if (paired==1){
			free(results->starts_reverse);
		}
		free(results->cigars_forward);
		if (paired==1){
			free(results->cigars_reverse);
		}
	}
	free(results->minimum);
	free(results->minNodes);
	free(results->LCAnames);
	free(results->leaf_coordinates);
	if (use_nw==1){
		needleman_wunsch_free(results->nw);
		alignment_free(results->aln);
		free(results->scoring);
	}
	free(results);
}
