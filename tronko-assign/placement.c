#include "placement.h"

int perform_WFA_alignment(cigar_t* const cigar, mm_allocator_t* mm_allocator,char* seq1, char* seq2,char* const pattern_alg,char* const text_alg, char* const ops_alg, int begin_offset, int end_offset){
	char* const operations = cigar->operations;
	int k, alg_pos =0, pattern_pos= 0, text_pos =0;
	for (k=begin_offset;k<end_offset;++k) {
		switch (operations[k]) {
			case 'M':
				if (seq1[pattern_pos] != seq2[text_pos]) {
					pattern_alg[alg_pos] = seq1[pattern_pos];
					ops_alg[alg_pos] = 'X';
					text_alg[alg_pos++] = seq2[text_pos];
				}else{
					pattern_alg[alg_pos] = seq1[pattern_pos];
					ops_alg[alg_pos] = '|';
					text_alg[alg_pos++] = seq2[text_pos];
				}
				pattern_pos++; text_pos++;
				break;
			case 'X':
				if (seq1[pattern_pos] != seq2[text_pos]) {
					pattern_alg[alg_pos] = seq1[pattern_pos++];
					ops_alg[alg_pos] = ' ';
					text_alg[alg_pos++] = seq2[text_pos++];
				}else{
					pattern_alg[alg_pos] = seq1[pattern_pos++];
					ops_alg[alg_pos] = 'X';
					text_alg[alg_pos++] = seq2[text_pos++];
				}
				break;
			case 'I':
				pattern_alg[alg_pos] = '-';
				ops_alg[alg_pos] = ' ';
				text_alg[alg_pos++] = seq2[text_pos++];
				break;
			case 'D':
				pattern_alg[alg_pos] = seq1[pattern_pos++];
				ops_alg[alg_pos] = ' ';
				text_alg[alg_pos++] = '-';
				break;
			default:
				break;
		}
	}
	k=0;
	while (pattern_pos < strlen(seq1)) {
		pattern_alg[alg_pos+k] = seq1[pattern_pos++];
		ops_alg[alg_pos+k] = '?';
		++k;
	}
	while (text_pos < strlen(seq2)) {
		text_alg[alg_pos+k] = seq2[text_pos++];
		ops_alg[alg_pos+k] = '?';
		++k;
	}
	int alignment_length = strlen(pattern_alg);
	return alignment_length;
}

void place_paired( char *query_1, char *query_2, char **rootSeqs, int numberOfTotalRoots, int *positions, char *locQuery, type_of_PP ***nodeScores, int **voteRoot, int number_of_matches , int **leaf_coordinates, int paired, type_of_PP* minimum_score, char *alignments_dir, char *forward_name, char *reverse_name, int print_alignments, char *leaf_sequence, int *positionsInRoot, int maxNumSpec, int* starts_forward, char** cigars_forward, int* starts_reverse, char** cigars_reverse, int print_alignments_to_file, int use_leaf_portion, int padding, int max_query_length, int max_numbase, int print_all_nodes){
	int i, j, k, node, match;
	type_of_PP forward_mismatch, reverse_mismatch;
	forward_mismatch=0;
	reverse_mismatch=0;
	/*affine_penalties.match = 0;
	affine_penalties.mismatch =4;
	affine_penalties.gap_opening=6;
	affine_penalties.gap_extension = 2;
	#define BUFFER_SIZE_8M   (1ul<<23)*/
	wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
	attributes.distance_metric = gap_affine;
	attributes.affine_penalties.mismatch =4;
	attributes.affine_penalties.gap_opening = 6;
	attributes.affine_penalties.gap_extension = 2;
	attributes.alignment_form.span = alignment_endsfree;
	//attributes.alignment_form.pattern_begin_free = 0;
	//attributes.alignment_form.pattern_end_free = 300;
	//attributes.alignment_form.text_begin_free = 50;
	//attributes.alignment_form.text_end_free = 50;
	for(match=0; match<number_of_matches; match++){
		if (use_leaf_portion==1){
			for(i=0; i<max_query_length+max_query_length+2*padding+1;i++){
				positions[i]=-1;
				locQuery[i]='\0';
				leaf_sequence[i] = '\0';
				positionsInRoot[i]=-1;
			}
		}else{
			for(i=0; i<max_query_length+max_numbase+1;i++){
				positions[i]=-1;
				locQuery[i]='\0';
				leaf_sequence[i] = '\0';
				positionsInRoot[i]=-1;
			}
		}
		int query_length = strlen(query_1);
		if (use_leaf_portion==1){
			if ( cigars_forward[match][0] == '*'){ break; }
			if ( cigars_forward[match][0] == '\0'){ break; }
		}
		if (use_leaf_portion == 1 && starts_forward[match] != -1){
			int start_position = getStartPosition(starts_forward[match],leaf_coordinates[match][0],leaf_coordinates[match][1],padding);
			int end_position = getEndPosition(cigars_forward[match],leaf_coordinates[match][0],leaf_coordinates[match][1],start_position+padding,padding);
			getSequenceInNodeWithoutNs(leaf_coordinates[match][0],leaf_coordinates[match][1],leaf_sequence,positionsInRoot,start_position,end_position);
		}else{
			getSequenceInNodeWithoutNs(leaf_coordinates[match][0],leaf_coordinates[match][1],leaf_sequence,positionsInRoot,0,numbaseArr[leaf_coordinates[match][0]]);
		}
		int leaf_length = strlen(leaf_sequence);
		if (leaf_length > 0){
		/*mm_allocator_t* mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
		affine_wavefronts_t* affine_wavefronts;
		char* const pattern_alg=mm_allocator_calloc(mm_allocator,leaf_length+query_length+1,char,true);
		char* const text_alg=mm_allocator_calloc(mm_allocator,leaf_length+query_length+1,char,true);*/
		// Initialize Wavefront Aligner
		wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
		wavefront_align(wf_aligner,leaf_sequence,leaf_length,query_1,query_length);
		char* const pattern_alg = mm_allocator_calloc(wf_aligner->mm_allocator,leaf_length+query_length+1,char,true);
		char* const ops_alg = mm_allocator_calloc(wf_aligner->mm_allocator,leaf_length+query_length+1,char,true);
		char* const text_alg = mm_allocator_calloc(wf_aligner->mm_allocator,leaf_length+query_length+1,char,true);
		int alignment_length = 0;
		/*affine_wavefronts = affine_wavefronts_new_complete(leaf_length,query_length,&affine_penalties,NULL,mm_allocator);
		affine_wavefronts_align(affine_wavefronts,leaf_sequence,leaf_length,query_1,query_length);*/
		//pattern_alg = mm_allocator_calloc(mm_allocator,leaf_length+query_length+1,char,true);
		//text_alg = mm_allocator_calloc(mm_allocator,leaf_length+query_length+1,char,true);
		alignment_length = perform_WFA_alignment(wf_aligner->cigar,wf_aligner->mm_allocator,leaf_sequence,query_1,pattern_alg,text_alg,ops_alg,wf_aligner->cigar->begin_offset,wf_aligner->cigar->end_offset);
		//clock_gettime(CLOCK_MONOTONIC, &tend);
		//printf("finished %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
		int alength=0;
		int keepTrackOfPosInRoot=0;
		int loop=0;
		int lengthOfResultA=0;
		/*for(i=0; aln->result_a[i] != '\0'; i++){
			if (aln->result_a[i]=='-'){*/
				/*positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[loop];
				loop++;
				alength++;*/
				//keepTrackOfPosInRoot++;
			/*}else if (aln->result_a[i]=='N'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[loop];
				loop++;
				alength++;
				keepTrackOfPosInRoot++;
			}else{
				if (alength != 0 && aln->result_b[loop]=='-' && aln->result_a[loop]!='-' && lengthOfResultA<query_length){
				//if ( aln->result_b[loop]=='-' && aln->result_a[loop]!='-' && lengthOfResultA<query_length){
					locQuery[alength]=aln->result_b[loop];
					positions[alength]=keepTrackOfPosInRoot;
					alength++;
				}
				if ( aln->result_b[loop] == '-' && aln->result_a[loop] != '-'){
					keepTrackOfPosInRoot++;
					loop++;
				}else{
					if (aln->result_a[i]!='-'){
						locQuery[alength]=aln->result_b[loop];
						positions[alength]=keepTrackOfPosInRoot;
						keepTrackOfPosInRoot++;
						alength++;
						loop++;
						lengthOfResultA++;
					}
				}
			}
		}*/
		/*for(i=0; aln->result_a[i] != '\0'; i++){
			if (keepTrackOfPosInRoot > leaf_length){
				positions[0]=-1;
				break;
			}
			//if (aln->result_a[i] != '-' && aln->result_b[i] != '-'){
			if (aln->result_a[i] != 'N' && aln->result_b[i] != '-' && aln->result_a[i] != '-'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[i];
				alength++;
				lengthOfResultA++;
				if (aln->result_a[i] != aln->result_b[i] && match==0){
					forward_mismatch++;
				}
			}
			if (alength != 0 && lengthOfResultA<query_length && aln->result_b[i] =='-'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[i];
				alength++;
				lengthOfResultA++;
				keepTrackOfPosInRoot++;
				if (aln->result_a[i] != aln->result_b[i] && match==0 && aln->result_a[i] != 'N' && aln->result_b[i] != '-' && aln->result_a[i] != '-'){
					forward_mismatch++;
				}
			}else if (aln->result_a[i] != '-'){
			//}else if (aln->result_a[i] != 'N'){
				keepTrackOfPosInRoot++;
			}
		}
		locQuery[alength]='\0';*/
		int counterInRoot=0;
		int counterinB=0;
			for (i=0; pattern_alg[i] != '\0'; i++){
				if (pattern_alg[i] != '-' && text_alg[i] != '-' ){
					positions[alength] = positionsInRoot[counterInRoot];
					locQuery[alength] = text_alg[i];
					alength++;
					if (pattern_alg[i] != text_alg[i] && match==0){
						forward_mismatch++;
					}
				}
				if ( pattern_alg[i] == '-' && text_alg[i] != '-' && alength > 0 && counterinB<query_length){
					positions[alength]=-1;
					locQuery[alength] = text_alg[i];
					alength++;
				}
				if (pattern_alg[i] != '-' && text_alg[i] == '-' && alength > 0 && counterinB<query_length){
					positions[alength] = positionsInRoot[counterInRoot];
					locQuery[alength] = text_alg[i];
					alength++;
				}
				if ( pattern_alg[i] != '-'){
					counterInRoot++;
				}
				if ( text_alg[i] != '-'){
					counterinB++;
				}
			}
		locQuery[alength]='\0';
		//if (access(alignmentFileName, F_OK ) != -1 && match==0){
		//	printToFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, forward_name);
		//}else if (match==0){
		//	createNewFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, forward_name);
		//}
		if (match==0 && print_alignments_to_file==1){
			char alignmentFileName[BUFFER_SIZE];
			snprintf(alignmentFileName,BUFFER_SIZE,"%s/%s.fasta",alignments_dir,treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			snprintf(alignmentFileName,BUFFER_SIZE,"%s",alignments_dir);
			printToFile_WFA(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, pattern_alg, text_alg, forward_name, query_length, positionsInRoot);
		}
		if (print_alignments==1){
			int breaks = 100;
			i=0;
			while(pattern_alg[i] != '\0'){
			printf("\n");
			printf("%s\n",forward_name);
			  printf("query_1\t\t\t\t\t");
			  for(i=breaks-100; i<breaks; i++){
				  if (pattern_alg[i]=='\0'){ break; }
					printf("%c",text_alg[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=breaks-100; i<breaks; i++){
				  if (pattern_alg[i]=='\0'){ break;}
					  printf("%c",pattern_alg[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=breaks-100; i<breaks; i++){
				  if (pattern_alg[i]=='\0'){break;}
					if ( pattern_alg[i] == text_alg[i] && pattern_alg[i] != '-' && text_alg[i] != '-' ){
						printf("|");
					}else if ( pattern_alg[i] != text_alg[i] && pattern_alg[i] != 'N' && text_alg[i] != '-' && pattern_alg[i] != '-'){
						printf("*");
					}else{
						printf(" ");
					}
			  }
			  printf("\n\n");
			  breaks += 100;
			}
			  /*printf("query_1\t\t\t\t\t");
			  for(i=101;i<200;i++){
					printf("%c",text_alg[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=101;i<200;i++){
					printf("%c",pattern_alg[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=101;i<200;i++){
					  if ( pattern_alg[i] == text_alg[i] && pattern_alg[i] != '-' && text_alg[i] != '-' ){
						  printf("|");
					  }else if ( pattern_alg[i] != text_alg[i] && pattern_alg[i] != 'N' && text_alg[i] != '-'){
						  printf("*");
					  }else{
						  printf(" ");
					  }
			  }
			  printf("\n\n");
			  printf("query_1\t\t\t\t\t");
			  for(i=201; i<300; i++){
			  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=201;i<300;i++){
			  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=201;i<300;i++){
			  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  printf("|");
			  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
			  printf("*");
			  }else{
			  printf(" ");
			  }
			  }
			  printf("\n\n");
			  printf("query_1\t\t\t\t\t");
			  for(i=301; i<400 ; i++){
			  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=301;i<400; i++){
			  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=301;i<400;i++){
			  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  printf("|");
			  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
				printf("*");
			  }else{
				  printf(" ");
			  }
			  }
			  printf("\n\n");
			  printf("query_1\t\t\t\t\t");
			  for(i=401; i<500 ; i++){
				printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=401; i<500 ; i++){
				  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=401; i<500 ; i++){
				  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					  printf("|");
					  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						  printf("*");
						  }else{
							  printf(" ");
						  }
			  }
			  printf("\n\n");
			  printf("query_1\t\t\t\t\t");
			  for(i=501; i<600; i++){
			  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=501;i<600; i++){
			  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=501;i<600;i++){
			  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  printf("|");
			  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
				  printf("*");
			  }else{
				  printf(" ");
			  }
			  }
			  printf("\n\n");
			  printf("query_1\t\t\t\t\t");
			  for(i=601; i<700; i++){
				  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=601;i<700; i++){
				  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=601;i<700;i++){
				  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					  printf("|");
					  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						  printf("*");
						  }else{
							  printf(" ");
						  }
			  }
			  printf("\n\n");
			printf("query_1\t\t\t\t\t");
			for(i=701; i<800; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=701;i<800; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=701;i<800;i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						printf("*");
					}else{
						printf(" ");
					}
			}
			printf("\n\n");
			printf("query_1\t\t\t\t\t");
			for(i=801;i<900; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=801;i<900; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=801;i<900;i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					 printf("|");
					  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						   printf("*");
					  }else{
						  printf(" ");
					  }
			}
			 printf("\n\n");
			printf("query_1\t\t\t\t\t");
			for(i=901; aln->result_b[i] != '\0'; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=901;aln->result_b[i]!='\0'; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=901;aln->result_b[i]!='\0';i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						printf("*");
					}else{
						printf(" ");
					}
			}
			printf("\n\n");*/
		}
		//assignScores_Arr_paired(topRoots[rootNum],rootArr[topRoots[rootNum]],locQuery, positions, nodeScores, alength);
		//clock_gettime(CLOCK_MONOTONIC, &tstart);
		//printf("assignScores_Arr_paired...\n");
		FILE* site_scores_file;
		if ( print_all_nodes == 1 && match > 0 && access("site_scores.txt", F_OK ) != -1 ){
			if (( site_scores_file = fopen("site_scores.txt","a")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
			fprintf(site_scores_file,"%s\t",forward_name);
		}
		if ( print_all_nodes == 1 && match > 0 ){
			if (( site_scores_file = fopen("site_scores.txt","w")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
			fprintf(site_scores_file,"Readname\tSite\tScore\n");
			fprintf(site_scores_file,"%s\t",forward_name);
		}
		assignScores_Arr_paired(leaf_coordinates[match][0],rootArr[leaf_coordinates[match][0]],locQuery, positions, nodeScores, alength, match,print_all_nodes,site_scores_file);
		//clock_gettime(CLOCK_MONOTONIC, &tend);
		if ( print_all_nodes == 1){
			fclose(site_scores_file);
		}
		//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
		mm_allocator_free(wf_aligner->mm_allocator,pattern_alg);
		mm_allocator_free(wf_aligner->mm_allocator,ops_alg);
		mm_allocator_free(wf_aligner->mm_allocator,text_alg);
		wavefront_aligner_delete(wf_aligner);
	}
	}
	if (paired==1){
	for(match=0; match<number_of_matches; match++){
		if (use_leaf_portion==1){
			for(i=0; i<max_query_length+max_query_length+2*padding+1;i++){
				positions[i]=-1;
				locQuery[i]='\0';
				leaf_sequence[i] = '\0';
				positionsInRoot[i]=-1;
			}
		}else{
			for(i=0; i<max_query_length+max_numbase+1;i++){
				positions[i]=-1;
				locQuery[i]='\0';
				leaf_sequence[i] = '\0';
				positionsInRoot[i]=-1;
			}
		}
		int query_length = strlen(query_2);
		if (use_leaf_portion==1){
			if ( cigars_reverse[match][0] == '*'){
				break;
			}
			if ( cigars_reverse[match][0] == '\0'){ break; }
		}
		if (use_leaf_portion==1){
			if (starts_reverse[match]!=1){
		//needleman_wunsch_align(rootSeqs[topRoots[rootNum]], query_2, scoring, nw, aln);
		//char *leaf_sequence = (char *)malloc(numbaseArr[leaf_coordinates[rootNum][0]]*sizeof(char));
		//getSequenceInNode(leaf_coordinates[match][0],leaf_coordinates[match][1],leaf_sequence);
				int start_position = getStartPosition(starts_reverse[match],leaf_coordinates[match][0],leaf_coordinates[match][1],padding);
				int end_position = getEndPosition(cigars_reverse[match],leaf_coordinates[match][0],leaf_coordinates[match][1],start_position+padding,padding);
				getSequenceInNodeWithoutNs(leaf_coordinates[match][0],leaf_coordinates[match][1],leaf_sequence,positionsInRoot,start_position,end_position);
			}
		}else{
			getSequenceInNodeWithoutNs(leaf_coordinates[match][0],leaf_coordinates[match][1],leaf_sequence,positionsInRoot,0,numbaseArr[leaf_coordinates[match][0]]);
		}
		int leaf_length = strlen(leaf_sequence);
		if (leaf_length > 0){
		//printf("leaf sequence: ");
		//for(i=0; i<strlen(leaf_sequence); i++){
		//	printf("%c",leaf_sequence[i]);
		//}
		//printf("\n");
		//clock_gettime(CLOCK_MONOTONIC, &tstart);
		//printf("aligning 2nd pair...\n");
		/*mm_allocator_t* mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
		affine_wavefronts_t* affine_wavefronts;
		char* const pattern_alg=mm_allocator_calloc(mm_allocator,leaf_length+query_length+1,char,true);
		char* const text_alg=mm_allocator_calloc(mm_allocator,leaf_length+query_length+1,char,true);*/
		wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
		wavefront_align(wf_aligner,leaf_sequence,leaf_length,query_2,query_length);
		char* const pattern_alg = mm_allocator_calloc(wf_aligner->mm_allocator,leaf_length+query_length+1,char,true);
		char* const ops_alg = mm_allocator_calloc(wf_aligner->mm_allocator,leaf_length+query_length+1,char,true);
		char* const text_alg = mm_allocator_calloc(wf_aligner->mm_allocator,leaf_length+query_length+1,char,true);
		int alignment_length=0;
		alignment_length = perform_WFA_alignment(wf_aligner->cigar,wf_aligner->mm_allocator,leaf_sequence,query_2,pattern_alg,text_alg,ops_alg,wf_aligner->cigar->begin_offset,wf_aligner->cigar->end_offset);
			/*affine_wavefronts = affine_wavefronts_new_complete(leaf_length,query_length,&affine_penalties,NULL,mm_allocator);
			affine_wavefronts_align(affine_wavefronts,leaf_sequence,leaf_length,query_2,query_length);
			alignment_length = perform_WFA_alignment(affine_wavefronts,mm_allocator,leaf_sequence,query_2,pattern_alg,text_alg);*/
		//clock_gettime(CLOCK_MONOTONIC, &tend);
		//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
		/*printf("leaf sequence: ");
		for(i=0; i<strlen(leaf_sequence); i++){
			printf("%c",leaf_sequence[i]);
		}*/
		int alength=0;
		int keepTrackOfPosInRoot=0;
		int loop=0;
		int lengthOfResultA=0;
		/*for(i=0; aln->result_a[i] != '\0'; i++){
			if (aln->result_a[i]=='-'){*/
				/*positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[loop];
				loop++;
				alength++;*/
				//keepTrackOfPosInRoot++;
			/*}else if (aln->result_a[i]=='N'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[loop];
				loop++;
				alength++;
				keepTrackOfPosInRoot++;
			}else{
				if (alength != 0 && aln->result_b[loop]=='-' && aln->result_a[loop]!='-' && lengthOfResultA<query_length){
				//if ( aln->result_b[loop]=='-' && aln->result_a[loop]!='-' && lengthOfResultA<query_length){
					locQuery[alength]=aln->result_b[loop];
					positions[alength]=keepTrackOfPosInRoot;
					alength++;
				}
				if ( aln->result_b[loop] == '-' && aln->result_a[loop] != '-'){
					keepTrackOfPosInRoot++;
					loop++;
				}else{
					if (aln->result_a[i]!='-'){
						locQuery[alength]=aln->result_b[loop];
						positions[alength]=keepTrackOfPosInRoot;
						keepTrackOfPosInRoot++;
						alength++;
						loop++;
						lengthOfResultA++;
					}
				}
			}
		}*/
		/*for(i=0; aln->result_a[i] != '\0'; i++){
			if (keepTrackOfPosInRoot > leaf_length){
				positions[0]=-1;
				break;
			}
			//if (aln->result_a[i] != '-' && aln->result_b[i] != '-'){
			if (aln->result_a[i] != 'N' && aln->result_b[i] != '-' && aln->result_a[i] != '-'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[i];
				alength++;
				lengthOfResultA++;
				if (aln->result_a[i] != aln->result_b[i] && match==0){
					reverse_mismatch++;
				}
			}
			if (alength != 0 && lengthOfResultA<query_length && aln->result_b[i]=='-'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[i];
				alength++;
				lengthOfResultA++;
				keepTrackOfPosInRoot++;
				if (aln->result_a[i] != aln->result_b[i] && match==0 && aln->result_a[i] != 'N' && aln->result_b[i] != '-' && aln->result_a[i] != '-'){
					reverse_mismatch++;
				}
			}else if (aln->result_a[i] != '-'){
			//}else if (aln->result_a[i] != 'N'){
				keepTrackOfPosInRoot++;
			}
		}
		locQuery[alength]='\0';*/
		int counterInRoot=0;
		int counterinB=0;
			for (i=0; pattern_alg[i] != '\0'; i++){
				if (pattern_alg[i] != '-' && text_alg[i] != '-' ){
					positions[alength] = positionsInRoot[counterInRoot];
					locQuery[alength] = text_alg[i];
					alength++;
					if (pattern_alg[i] != text_alg[i] && match==0){
						reverse_mismatch++;
					}
				}
				if ( pattern_alg[i] == '-' && text_alg[i] != '-' && alength > 0 && counterinB<query_length){
					positions[alength]=-1;
					locQuery[alength] = text_alg[i];
					alength++;
				}
				if (pattern_alg[i] != '-' && text_alg[i] == '-' && alength > 0 && counterinB<query_length){
					positions[alength] = positionsInRoot[counterInRoot];
					locQuery[alength] = text_alg[i];
					alength++;
				}
				if ( pattern_alg[i] != '-'){
					counterInRoot++;
				}
				if ( text_alg[i] != '-'){
					counterinB++;
				}
			}
		locQuery[alength]='\0';
		//if (access(alignmentFileName, F_OK ) != -1 && match==0){
		//	printToFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, reverse);
		//}else if (match==0){
		//	createNewFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, readname);
		//}
		if (match==0 && print_alignments_to_file==1){
			char alignmentFileName[1000];
			snprintf(alignmentFileName,1000,"%s/%s.fasta",alignments_dir,treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			printToFile_WFA(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, pattern_alg, text_alg, reverse_name, query_length, positionsInRoot);
		}
		if (print_alignments==1){
			int breaks = 100;
			i=0;
			while(pattern_alg[i] != '\0'){
			printf("\n");
			printf("%s\n",reverse_name);
			printf("query_2\t\t\t\t\t");
			for(i=breaks-100; i<breaks; i++){
				if (text_alg[i]=='\0'){ break; }
					printf("%c",text_alg[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=breaks-100; i<breaks; i++){
				if (pattern_alg[i]=='\0'){ break; }
					printf("%c",pattern_alg[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=breaks-100; i<breaks; i++){
				if (pattern_alg[i]=='\0'){ break; }
					if ( pattern_alg[i] == text_alg[i] && pattern_alg[i] != '-' && text_alg[i] != '-' ){
						printf("|");
					}else if ( pattern_alg[i] != text_alg[i] && pattern_alg[i] != 'N' && text_alg[i] != '-' && pattern_alg[i] != '-'){
						printf("*");
					}else{
						printf(" ");
					}
			}
			printf("\n\n");
			breaks += 100;
			}
			/*
			printf("query_2\t\t\t\t\t");
			for(i=101;i<200;i++){
			  printf("%c",text_alg[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=101;i<200;i++){
			  printf("%c",pattern_alg[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=101;i<200;i++){
			  if(pattern_alg[i] == text_alg[i] && pattern_alg[i] !='-' && text_alg[i] != '-'){
			  printf("|");
			  }else if ( pattern_alg[i] != text_alg[i] && pattern_alg[i] != 'N' && text_alg[i] != '-' ){
			  printf("*");
			  }else{
			  printf(" ");
			  }
			}
			printf("\n\n");
			  printf("query_2\t\t\t\t\t");
			  for(i=201; i<300; i++){
			  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=201;i<300;i++){
			  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=201;i<300;i++){
			  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  printf("|");
			  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
			  printf("*");
			  }else{
			  printf(" ");
			  }
			  }
			  printf("\n\n");
			  printf("\nquery_2\t\t\t\t\t");
			  for(i=301; i<400; i++){
			  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=301;i<400; i++){
			  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=301;i<400;i++){
			  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  printf("|");
			  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
				  printf("*");
			  }else{
				  printf(" ");
			  }
			  }
			  printf("\n\n");
			  printf("\nquery_2\t\t\t\t\t");
			  for(i=401; i<500; i++){
				  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=401; i<500; i++){
				  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=401; i<500; i++){
				  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					  printf("|");
					  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						  printf("*");
						  }else{
							  printf(" ");
						  }
			  }
			  printf("\n\n");
			  printf("\nquery_2\t\t\t\t\t");
			for(i=501;i<600; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=501;i<600; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=501;i<600; i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						printf("*");
						}else{
							printf(" ");
						}
			}
			printf("\n\n");
			printf("\nquery_2\t\t\t\t\t");
			for(i=601;i<700; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=601;i<700; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=601;i<700; i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						printf("*");
						}else{
							printf(" ");
						}
			}
			printf("\n\n");
			printf("\nquery_2\t\t\t\t\t");
			for(i=701;i<800; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=701;i<800; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=701;i<800; i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						printf("*");
						}else{
							printf(" ");
						}
			}
			printf("\n\n");
			printf("\nquery_2\t\t\t\t\t");
			for(i=801;i<900; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=801;i<900; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=801;i<900; i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='-' && aln->result_b[i]!='-' && aln->result_a[i]!='N'){
						printf("*");
						}else{
							printf(" ");
						}
			}
			printf("\n\n");
			printf("\nquery_2\t\t\t\t\t");
			for(i=901; aln->result_b[i]!='\0'; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=901;aln->result_b[i]!='\0'; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=901;aln->result_b[i]!='\0'; i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='-' && aln->result_b[i]!='-' && aln->result_a[i]!='N'){
						printf("*");
						}else{
							printf(" ");
						}
			}
			printf("\n\n");*/
		//}
		//clock_gettime(CLOCK_MONOTONIC, &tstart);
		//printf("assigningscores to 2nd pair...\n");
		}
		FILE* site_scores_file;
		if ( print_all_nodes == 1 ){
			if (( site_scores_file = fopen("site_scores.txt","a")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
			fprintf(site_scores_file,"%s\t",reverse_name);
		}
		assignScores_Arr_paired(leaf_coordinates[match][0],rootArr[leaf_coordinates[match][0]],locQuery,positions,nodeScores,alength,match,print_all_nodes,site_scores_file);
		//clock_gettime(CLOCK_MONOTONIC, &tend);
		//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
		//free(leaf_sequence);
		if ( print_all_nodes == 1 ){
			fclose(site_scores_file);
		}
		mm_allocator_free(wf_aligner->mm_allocator,pattern_alg);
		mm_allocator_free(wf_aligner->mm_allocator,ops_alg);
		mm_allocator_free(wf_aligner->mm_allocator,text_alg);
		wavefront_aligner_delete(wf_aligner);
		}
		}
	}
	type_of_PP maximum=-9999999999999999;
	int minRoot=0;
	int minNode=0;
	int match_number=0;
	//clock_gettime(CLOCK_MONOTONIC, &tstart);
	//printf("finding minimum score...\n");
	FILE* node_scores_file;
	if ( print_all_nodes == 1 ){
		if (( node_scores_file = fopen("scores_all_nodes.txt","w")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
		fprintf(node_scores_file,"Tree_Number\tNode_Number\tScore\n");
	}
	for (i=0; i<number_of_matches;i++){
		for (j=leaf_coordinates[i][0]; j<leaf_coordinates[i][0]+1; j++){
			for(k=0; k<2*numspecArr[j]-1; k++){
				//printf("Tree %d Node %d Taxonomy (%s) Score: %lf\n",j,k,taxonomyArr[j][treeArr[j][k].taxIndex[0]][treeArr[j][k].taxIndex[1]],nodeScores[i][j][k]);
				if ( maximum < nodeScores[i][j][k]){
					maximum=nodeScores[i][j][k];
					match_number=i;
					minRoot=j;
					minNode=k;
				}
				if ( print_all_nodes == 1){
					fprintf(node_scores_file,"%d\t%d\t%lf\n",j,k,nodeScores[i][j][k]);
				}
			}
		}
	}
	if (print_all_nodes == 1){
		fclose(node_scores_file);
	}
	//clock_gettime(CLOCK_MONOTONIC, &tend);
	//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	//printf("match number: %d\n",match_number);
	//printf("minimum score: %lf\n",maximum);
	//printf("minimum root: %d\n",minRoot);
	//printf("minimum node: %d\n",minNode);
	//printf("C interval: %lf\n",Cinterval);
	//clock_gettime(CLOCK_MONOTONIC, &tstart);
	//printf("clearing voteroot...\n");
	/*for(i=0; i<numberOfTotalRoots; i++){
		for(j=0;j<2*numspecArr[i]-1;j++){
			voteRoot[i][j]=0;
		}
	}*/
	//clock_gettime(CLOCK_MONOTONIC, &tend);
	//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	int index = 0;
	//clock_gettime(CLOCK_MONOTONIC, &tstart);
	//printf("filling out voteroot...\n");
	for(i=0; i<number_of_matches; i++){
		//for(j=leaf_coordinates[i][0]; j<leaf_coordinates[i][0]+1; j++){
			for(k=0; k<2*numspecArr[leaf_coordinates[i][0]]-1; k++){
				if ( nodeScores[i][leaf_coordinates[i][0]][k] >= (maximum-Cinterval) && nodeScores[i][leaf_coordinates[i][0]][k] <= (maximum+Cinterval) ){
					//printf("Match : %d Min Root: %d Min node: %d, score: %lf\n",i,j,k,nodeScores[i][leaf_coordinates[i][0]][k]);
					voteRoot[leaf_coordinates[i][0]][k]=1;
					index++;
				}
			}
		//}
	}
	//clock_gettime(CLOCK_MONOTONIC, &tend);
	//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	minimum_score[0] = maximum;
	minimum_score[1] = forward_mismatch;
	minimum_score[2] = reverse_mismatch;
	//printf("%lf\t",minimum);
	//free(leaf_sequence);
	//free(positionsInRoot);
}
void place_paired_with_nw( char *query_1, char *query_2, char **rootSeqs, int numberOfTotalRoots, int *positions, char *locQuery, nw_aligner_t *nw, alignment_t *aln, scoring_t *scoring, type_of_PP ***nodeScores, int **voteRoot, int number_of_matches , int **leaf_coordinates, int paired, type_of_PP* minimum_score, char *alignments_dir, char *forward_name, char *reverse_name, int print_alignments, char *leaf_sequence, int *positionsInRoot, int maxNumSpec, int* starts_forward, char** cigars_forward, int* starts_reverse, char** cigars_reverse, int print_alignments_to_file, int use_leaf_portion, int padding, int max_query_length, int max_numbase, int print_all_nodes){
	int i, j, k, node, match;
	type_of_PP forward_mismatch, reverse_mismatch;
	forward_mismatch=0;
	reverse_mismatch=0;
	for(match=0; match<number_of_matches; match++){
		if (use_leaf_portion==1){
			for(i=0; i<max_query_length+max_query_length+2*padding+1;i++){
				positions[i]=-1;
				locQuery[i]='\0';
				leaf_sequence[i] = '\0';
				positionsInRoot[i]=-1;
			}
		}else{
			for(i=0; i<max_query_length+max_numbase+1;i++){
				positions[i]=-1;
				locQuery[i]='\0';
				leaf_sequence[i] = '\0';
				positionsInRoot[i]=-1;
			}
		}
		int query_length = strlen(query_1);
		if (use_leaf_portion==1){
			if ( cigars_forward[match][0] == '*'){ break; }
			if ( cigars_forward[match][0] == '\0'){ break; }
		}
		if (use_leaf_portion == 1 && starts_forward[match] != -1){
			int start_position = getStartPosition(starts_forward[match],leaf_coordinates[match][0],leaf_coordinates[match][1],padding);
			int end_position = getEndPosition(cigars_forward[match],leaf_coordinates[match][0],leaf_coordinates[match][1],start_position+padding,padding);
			getSequenceInNodeWithoutNs(leaf_coordinates[match][0],leaf_coordinates[match][1],leaf_sequence,positionsInRoot,start_position,end_position);
		}else{
			getSequenceInNodeWithoutNs(leaf_coordinates[match][0],leaf_coordinates[match][1],leaf_sequence,positionsInRoot,0,numbaseArr[leaf_coordinates[match][0]]);
		}
		int leaf_length = strlen(leaf_sequence);
		if (leaf_length > 0){
		/*printf("leaf sequence: ");
		for(i=0; i<strlen(leaf_sequence); i++){
			printf("%c",leaf_sequence[i]);
		}*/
		//printf("\n");
		//printf("aligning...\n");
		//clock_gettime(CLOCK_MONOTONIC, &tstart);
		needleman_wunsch_align(leaf_sequence, query_1, scoring, nw, aln);
		//clock_gettime(CLOCK_MONOTONIC, &tend);
		//printf("finished %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
		int alength=0;
		int keepTrackOfPosInRoot=0;
		int loop=0;
		int lengthOfResultA=0;
		/*for(i=0; aln->result_a[i] != '\0'; i++){
			if (aln->result_a[i]=='-'){*/
				/*positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[loop];
				loop++;
				alength++;*/
				//keepTrackOfPosInRoot++;
			/*}else if (aln->result_a[i]=='N'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[loop];
				loop++;
				alength++;
				keepTrackOfPosInRoot++;
			}else{
				if (alength != 0 && aln->result_b[loop]=='-' && aln->result_a[loop]!='-' && lengthOfResultA<query_length){
				//if ( aln->result_b[loop]=='-' && aln->result_a[loop]!='-' && lengthOfResultA<query_length){
					locQuery[alength]=aln->result_b[loop];
					positions[alength]=keepTrackOfPosInRoot;
					alength++;
				}
				if ( aln->result_b[loop] == '-' && aln->result_a[loop] != '-'){
					keepTrackOfPosInRoot++;
					loop++;
				}else{
					if (aln->result_a[i]!='-'){
						locQuery[alength]=aln->result_b[loop];
						positions[alength]=keepTrackOfPosInRoot;
						keepTrackOfPosInRoot++;
						alength++;
						loop++;
						lengthOfResultA++;
					}
				}
			}
		}*/
		/*for(i=0; aln->result_a[i] != '\0'; i++){
			if (keepTrackOfPosInRoot > leaf_length){
				positions[0]=-1;
				break;
			}
			//if (aln->result_a[i] != '-' && aln->result_b[i] != '-'){
			if (aln->result_a[i] != 'N' && aln->result_b[i] != '-' && aln->result_a[i] != '-'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[i];
				alength++;
				lengthOfResultA++;
				if (aln->result_a[i] != aln->result_b[i] && match==0){
					forward_mismatch++;
				}
			}
			if (alength != 0 && lengthOfResultA<query_length && aln->result_b[i] =='-'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[i];
				alength++;
				lengthOfResultA++;
				keepTrackOfPosInRoot++;
				if (aln->result_a[i] != aln->result_b[i] && match==0 && aln->result_a[i] != 'N' && aln->result_b[i] != '-' && aln->result_a[i] != '-'){
					forward_mismatch++;
				}
			}else if (aln->result_a[i] != '-'){
			//}else if (aln->result_a[i] != 'N'){
				keepTrackOfPosInRoot++;
			}
		}
		locQuery[alength]='\0';*/
		int counterInRoot=0;
		int counterinB=0;
			for(i=0; aln->result_a[i] != '\0'; i++){
				if (aln->result_a[i] != '-' && aln->result_b[i] != '-'){
					positions[alength] = positionsInRoot[counterInRoot];
					locQuery[alength] = aln->result_b[i];
					alength++;
					if (aln->result_a[i] != aln->result_b[i] && match==0){
						forward_mismatch++;
					}
				}
				if (aln->result_a[i] == '-' && aln->result_b[i] != '-' && alength>0 && counterinB<query_length){
					positions[alength]=-1;
					locQuery[alength] = aln->result_b[i];
					alength++;	
				}
				if (aln->result_a[i] != '-' && aln->result_b[i] == '-' && alength > 0 && counterinB<query_length){
					positions[alength]= positionsInRoot[counterInRoot];
					locQuery[alength] = aln->result_b[i];
					alength++;
				}
				if (aln->result_a[i] != '-'){
					counterInRoot++;
				}
				if (aln->result_b[i] != '-'){
					counterinB++;
				}
			}
		locQuery[alength]='\0';
		//if (access(alignmentFileName, F_OK ) != -1 && match==0){
		//	printToFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, forward_name);
		//}else if (match==0){
		//	createNewFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, forward_name);
		//}
		if (match==0 && print_alignments_to_file==1){
			char alignmentFileName[1000];
			snprintf(alignmentFileName,1000,"%s/%s.fasta",alignments_dir,treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			//snprintf(alignmentFileName,1000,"%s",alignments_dir);
			if ( access(alignmentFileName, F_OK ) != -1 ){
				printToFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, forward_name, query_length, positionsInRoot);
			}else{
				createNewFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, forward_name, query_length);
			}
		}
		if (print_alignments==1){
			printf("Using Needleman-Wunsch\n");
			printf("\n");
			printf("%s\n",forward_name);
			  printf("query_1\t\t\t\t\t");
			  for(i=0; i<100; i++){
			  		printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=0; i<100; i++){
					printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=0; i<100; i++){
			  		if ( aln->result_a[i] == aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  			printf("|");
			  		}else if (aln->result_a[i] != aln->result_b[i] && aln->result_a[i] !='N' && aln->result_b[i] != '-'){
			  			printf("*");
			  		}else{
			  			printf(" ");
			  		}
			  }
			  printf("\n\n");
			printf("%s\n",forward_name);
			  printf("query_1\t\t\t\t\t");
			  for(i=101;i<200;i++){
				  	printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=101;i<200;i++){
					printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=101;i<200;i++){
					  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
						  printf("|");
					  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						  printf("*");
					  }else{
						  printf(" ");
					  }
			  }
			  printf("\n\n");
			printf("%s\n",forward_name);
			  printf("query_1\t\t\t\t\t");
			  for(i=201; i<300; i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=201;i<300;i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=201;i<300;i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  printf("|");
			  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
			  printf("*");
			  }else{
			  printf(" ");
			  }
			  }
			  printf("\n\n");
			printf("%s\n",forward_name);
			  printf("query_1\t\t\t\t\t");
			  for(i=301; i<400 ; i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=301;i<400; i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=301;i<400;i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  printf("|");
			  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
				printf("*");
			  }else{
				  printf(" ");
			  }
			  }
			  printf("\n\n");
			  /*
			  printf("query_1\t\t\t\t\t");
			  for(i=401; i<500 ; i++){
				printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=401; i<500 ; i++){
				  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=401; i<500 ; i++){
				  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					  printf("|");
					  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						  printf("*");
						  }else{
							  printf(" ");
						  }
			  }
			  printf("\n\n");
			  printf("query_1\t\t\t\t\t");
			  for(i=501; i<600; i++){
			  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=501;i<600; i++){
			  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=501;i<600;i++){
			  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  printf("|");
			  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
				  printf("*");
			  }else{
				  printf(" ");
			  }
			  }
			  printf("\n\n");
			  printf("query_1\t\t\t\t\t");
			  for(i=601; i<700; i++){
				  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=601;i<700; i++){
				  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=601;i<700;i++){
				  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					  printf("|");
					  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						  printf("*");
						  }else{
							  printf(" ");
						  }
			  }
			  printf("\n\n");
			printf("query_1\t\t\t\t\t");
			for(i=701; i<800; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=701;i<800; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=701;i<800;i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						printf("*");
					}else{
						printf(" ");
					}
			}
			printf("\n\n");
			printf("query_1\t\t\t\t\t");
			for(i=801;i<900; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=801;i<900; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=801;i<900;i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					 printf("|");
					  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						   printf("*");
					  }else{
						  printf(" ");
					  }
			}
			 printf("\n\n");
			printf("query_1\t\t\t\t\t");
			for(i=901; aln->result_b[i] != '\0'; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=901;aln->result_b[i]!='\0'; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=901;aln->result_b[i]!='\0';i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						printf("*");
					}else{
						printf(" ");
					}
			}
			printf("\n\n");*/
		}
		//assignScores_Arr_paired(topRoots[rootNum],rootArr[topRoots[rootNum]],locQuery, positions, nodeScores, alength);
		//clock_gettime(CLOCK_MONOTONIC, &tstart);
		//printf("assignScores_Arr_paired...\n");
		FILE* site_scores_file;
		if ( print_all_nodes == 1 && match > 0 && access("site_scores.txt", F_OK ) != -1 ){
			if (( site_scores_file = fopen("site_scores.txt","a")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
			fprintf(site_scores_file,"%s\t",forward_name);
		} else if ( print_all_nodes == 1 && match ==0 ){
			if (( site_scores_file = fopen("site_scores.txt","w")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
			fprintf(site_scores_file,"Readname\tSite\tScore\n");
			fprintf(site_scores_file,"%s\t",forward_name);
		}	
		assignScores_Arr_paired(leaf_coordinates[match][0],rootArr[leaf_coordinates[match][0]],locQuery, positions, nodeScores, alength, match, print_all_nodes, site_scores_file);
		if ( print_all_nodes == 1){
			fclose(site_scores_file);
		}
		//clock_gettime(CLOCK_MONOTONIC, &tend);
		//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	}
	}
	if (paired==1){
	for(match=0; match<number_of_matches; match++){
		if (use_leaf_portion==1){
			for(i=0; i<max_query_length+2*padding+1;i++){
				positions[i]=-1;
				locQuery[i]='\0';
				leaf_sequence[i] = '\0';
				positionsInRoot[i]=-1;
			}
		}else{
			for(i=0; i<max_query_length+max_numbase+1;i++){
				positions[i]=-1;
				locQuery[i]='\0';
				leaf_sequence[i] = '\0';
				positionsInRoot[i]=-1;
			}
		}
		int query_length = strlen(query_2);
		//needleman_wunsch_align(rootSeqs[topRoots[rootNum]], query_2, scoring, nw, aln);
		//char *leaf_sequence = (char *)malloc(numbaseArr[leaf_coordinates[rootNum][0]]*sizeof(char));
		//getSequenceInNode(leaf_coordinates[match][0],leaf_coordinates[match][1],leaf_sequence);
		if (use_leaf_portion==1){
			if ( cigars_reverse[match][0] == '*'){
				break;
			}
			if ( cigars_forward[match][0] == '\0'){ break; }
		}
		if (use_leaf_portion == 1 ){
			if ( starts_reverse[match]!=1){
				int start_position = getStartPosition(starts_reverse[match],leaf_coordinates[match][0],leaf_coordinates[match][1],padding);
				int end_position = getEndPosition(cigars_reverse[match],leaf_coordinates[match][0],leaf_coordinates[match][1],start_position+padding,padding);
				getSequenceInNodeWithoutNs(leaf_coordinates[match][0],leaf_coordinates[match][1],leaf_sequence,positionsInRoot,start_position,end_position);
			}
		}else{
			getSequenceInNodeWithoutNs(leaf_coordinates[match][0],leaf_coordinates[match][1],leaf_sequence,positionsInRoot,0,numbaseArr[leaf_coordinates[match][0]]);
		}
		int leaf_length = strlen(leaf_sequence);
		if(leaf_length > 0){
		//printf("leaf sequence: ");
		//for(i=0; i<strlen(leaf_sequence); i++){
		//	printf("%c",leaf_sequence[i]);
		//}
		//printf("\n");
		//clock_gettime(CLOCK_MONOTONIC, &tstart);
		//printf("aligning 2nd pair...\n");
		needleman_wunsch_align(leaf_sequence, query_2, scoring, nw, aln);
		//clock_gettime(CLOCK_MONOTONIC, &tend);
		//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
		/*printf("leaf sequence: ");
		for(i=0; i<strlen(leaf_sequence); i++){
			printf("%c",leaf_sequence[i]);
		}*/
		int alength=0;
		int keepTrackOfPosInRoot=0;
		int loop=0;
		int lengthOfResultA=0;
		/*for(i=0; aln->result_a[i] != '\0'; i++){
			if (aln->result_a[i]=='-'){*/
				/*positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[loop];
				loop++;
				alength++;*/
				//keepTrackOfPosInRoot++;
			/*}else if (aln->result_a[i]=='N'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[loop];
				loop++;
				alength++;
				keepTrackOfPosInRoot++;
			}else{
				if (alength != 0 && aln->result_b[loop]=='-' && aln->result_a[loop]!='-' && lengthOfResultA<query_length){
				//if ( aln->result_b[loop]=='-' && aln->result_a[loop]!='-' && lengthOfResultA<query_length){
					locQuery[alength]=aln->result_b[loop];
					positions[alength]=keepTrackOfPosInRoot;
					alength++;
				}
				if ( aln->result_b[loop] == '-' && aln->result_a[loop] != '-'){
					keepTrackOfPosInRoot++;
					loop++;
				}else{
					if (aln->result_a[i]!='-'){
						locQuery[alength]=aln->result_b[loop];
						positions[alength]=keepTrackOfPosInRoot;
						keepTrackOfPosInRoot++;
						alength++;
						loop++;
						lengthOfResultA++;
					}
				}
			}
		}*/
		/*for(i=0; aln->result_a[i] != '\0'; i++){
			if (keepTrackOfPosInRoot > leaf_length){
				positions[0]=-1;
				break;
			}
			//if (aln->result_a[i] != '-' && aln->result_b[i] != '-'){
			if (aln->result_a[i] != 'N' && aln->result_b[i] != '-' && aln->result_a[i] != '-'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[i];
				alength++;
				lengthOfResultA++;
				if (aln->result_a[i] != aln->result_b[i] && match==0){
					reverse_mismatch++;
				}
			}
			if (alength != 0 && lengthOfResultA<query_length && aln->result_b[i]=='-'){
				positions[alength]=keepTrackOfPosInRoot;
				locQuery[alength]=aln->result_b[i];
				alength++;
				lengthOfResultA++;
				keepTrackOfPosInRoot++;
				if (aln->result_a[i] != aln->result_b[i] && match==0 && aln->result_a[i] != 'N' && aln->result_b[i] != '-' && aln->result_a[i] != '-'){
					reverse_mismatch++;
				}
			}else if (aln->result_a[i] != '-'){
			//}else if (aln->result_a[i] != 'N'){
				keepTrackOfPosInRoot++;
			}
		}
		locQuery[alength]='\0';*/
		int counterInRoot=0;
		int counterinB=0;
			for(i=0; aln->result_a[i] != '\0'; i++){
				if (aln->result_a[i] != '-' && aln->result_b[i] != '-'){
					positions[alength] = positionsInRoot[counterInRoot];
					locQuery[alength] = aln->result_b[i];
					alength++;
					if (aln->result_a[i] != aln->result_b[i] && match==0){
						reverse_mismatch++;
					}
				}
				if (aln->result_a[i] == '-' && aln->result_b[i] != '-' && alength>0 && counterinB<query_length){
					positions[alength]=-1;
					locQuery[alength] = aln->result_b[i];
					alength++;
				}
				if (aln->result_a[i] != '-' && aln->result_b[i] == '-' && alength > 0 && counterinB <query_length){
					positions[alength]= positionsInRoot[counterInRoot];
					locQuery[alength] = aln->result_b[i];
					alength++;
				}
				if (aln->result_a[i] != '-'){
					counterInRoot++;
				}
				if (aln->result_b[i] != '-'){
					counterinB++;
				}
			}
		locQuery[alength]='\0';
		//if (access(alignmentFileName, F_OK ) != -1 && match==0){
		//	printToFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, reverse);
		//}else if (match==0){
		//	createNewFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, readname);
		//}
		if (match==0 && print_alignments_to_file==1){
			char alignmentFileName[1000];
			snprintf(alignmentFileName,1000,"%s/%s.fasta",alignments_dir,treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			if ( access(alignmentFileName, F_OK ) != -1 ){
				printToFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, reverse_name, query_length, positionsInRoot);
			}else{
				createNewFile2(treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name, alignments_dir, aln, reverse_name, query_length);
			}
		}
		if (print_alignments==1){
			printf("\n");
			printf("%s\n",reverse_name);
			printf("query_2\t\t\t\t\t");
			for(i=0; i<100; i++){
					printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=0; i<100; i++){
					printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t");
			for(i=0; i<100; i++){
					if ( aln->result_a[i] == aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
						printf("|");
					}else if (aln->result_a[i] !='N' && aln->result_b[i] != '-'){
						printf("*");
					}else{
						printf(" ");
					}
			}
			printf("\n\n");
			printf("%s\n",reverse_name);
			printf("query_2\t\t\t\t\t");
			  for(i=101;i<200;i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=101;i<200;i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=101;i<200;i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  printf("|");
			  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
			  printf("*");
			  }else{
			  printf(" ");
			  }
			  }
			  printf("\n\n");
			printf("%s\n",reverse_name);
			  printf("query_2\t\t\t\t\t");
			  for(i=201; i<300; i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=201;i<300;i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=201;i<300;i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  printf("|");
			  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
			  printf("*");
			  }else{
			  printf(" ");
			  }
			  }
			  printf("\n\n");
			printf("%s\n",reverse_name);
			  printf("\nquery_2\t\t\t\t\t");
			  for(i=301; i<400; i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=301;i<400; i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=301;i<400;i++){
				  if (aln->result_a[i]=='\0'){ break; }
			  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
			  printf("|");
			  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
				  printf("*");
			  }else{
				  printf(" ");
			  }
			  }
			  printf("\n\n");
			  /*
			  printf("\nquery_2\t\t\t\t\t");
			  for(i=401; i<500; i++){
				  printf("%c",aln->result_b[i]);
			  }
			  printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			  for(i=401; i<500; i++){
				  printf("%c",aln->result_a[i]);
			  }
			  printf("\n\t\t\t\t\t");
			  for(i=401; i<500; i++){
				  if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					  printf("|");
					  }else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						  printf("*");
						  }else{
							  printf(" ");
						  }
			  }
			  printf("\n\n");
			  printf("\nquery_2\t\t\t\t\t");
			for(i=501;i<600; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=501;i<600; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=501;i<600; i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						printf("*");
						}else{
							printf(" ");
						}
			}
			printf("\n\n");
			printf("\nquery_2\t\t\t\t\t");
			for(i=601;i<700; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=601;i<700; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=601;i<700; i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						printf("*");
						}else{
							printf(" ");
						}
			}
			printf("\n\n");
			printf("\nquery_2\t\t\t\t\t");
			for(i=701;i<800; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=701;i<800; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=701;i<800; i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='N' && aln->result_b[i]!='-'){
						printf("*");
						}else{
							printf(" ");
						}
			}
			printf("\n\n");
			printf("\nquery_2\t\t\t\t\t");
			for(i=801;i<900; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=801;i<900; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=801;i<900; i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='-' && aln->result_b[i]!='-' && aln->result_a[i]!='N'){
						printf("*");
						}else{
							printf(" ");
						}
			}
			printf("\n\n");
			printf("\nquery_2\t\t\t\t\t");
			for(i=901; aln->result_b[i]!='\0'; i++){
				printf("%c",aln->result_b[i]);
			}
			printf("\nmatch(%d) root(%d) %s\t\t",match,leaf_coordinates[match][0],treeArr[leaf_coordinates[match][0]][leaf_coordinates[match][1]].name);
			for(i=901;aln->result_b[i]!='\0'; i++){
				printf("%c",aln->result_a[i]);
			}
			printf("\n\t\t\t\t\t");
			for(i=901;aln->result_b[i]!='\0'; i++){
				if(aln->result_a[i]==aln->result_b[i] && aln->result_a[i]!='-' && aln->result_b[i]!='-'){
					printf("|");
					}else if (aln->result_a[i]!='-' && aln->result_b[i]!='-' && aln->result_a[i]!='N'){
						printf("*");
						}else{
							printf(" ");
						}
			}
			printf("\n\n");*/
		//}
		//clock_gettime(CLOCK_MONOTONIC, &tstart);
		//printf("assigningscores to 2nd pair...\n");
		}
		FILE* site_scores_file;
		if ( print_all_nodes == 1 ){
			if (( site_scores_file = fopen("site_scores.txt","a")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
			fprintf(site_scores_file,"%s\t",reverse_name);
		}
		assignScores_Arr_paired(leaf_coordinates[match][0],rootArr[leaf_coordinates[match][0]],locQuery,positions,nodeScores,alength,match,print_all_nodes,site_scores_file);
		if ( print_all_nodes == 1){
			fclose(site_scores_file);
		}
		//clock_gettime(CLOCK_MONOTONIC, &tend);
		//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
		//free(leaf_sequence);
	}
		}
	}
	type_of_PP maximum=-9999999999999999;
	int minRoot=0;
	int minNode=0;
	int match_number=0;
	//clock_gettime(CLOCK_MONOTONIC, &tstart);
	//printf("finding minimum score...\n");
	FILE* node_scores_file;
	if ( print_all_nodes == 1 ){
		if (( node_scores_file = fopen("scores_all_nodes.txt","w")) == (FILE *) NULL ) fprintf(stderr, "File could not be opened.\n");
		fprintf(node_scores_file,"Tree_Number\tNode_Number\tScore\n");
	}
	for (i=0; i<number_of_matches;i++){
		for (j=leaf_coordinates[i][0]; j<leaf_coordinates[i][0]+1; j++){
			for(k=0; k<2*numspecArr[j]-1; k++){
				if ( maximum < nodeScores[i][j][k]){
					maximum=nodeScores[i][j][k];
					match_number=i;
					minRoot=j;
					minNode=k;
				}
				if ( print_all_nodes == 1){
					fprintf(node_scores_file,"%d\t%d\t%lf\n",j,k,nodeScores[i][j][k]);
				}
			}
		}
	}
	if (print_all_nodes == 1){
		fclose(node_scores_file);
	}
	//clock_gettime(CLOCK_MONOTONIC, &tend);
	//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	//printf("match number: %d\n",match_number);
	//printf("minimum score: %lf\n",maximum);
	//printf("minimum root: %d\n",minRoot);
	//printf("minimum node: %d\n",minNode);
	//printf("C interval: %lf\n",Cinterval);
	//clock_gettime(CLOCK_MONOTONIC, &tstart);
	//printf("clearing voteroot...\n");
	/*for(i=0; i<numberOfTotalRoots; i++){
		for(j=0;j<2*numspecArr[i]-1;j++){
			voteRoot[i][j]=0;
		}
	}*/
	//clock_gettime(CLOCK_MONOTONIC, &tend);
	//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	int index = 0;
	//clock_gettime(CLOCK_MONOTONIC, &tstart);
	//printf("filling out voteroot...\n");
	for(i=0; i<number_of_matches; i++){
		//for(j=leaf_coordinates[i][0]; j<leaf_coordinates[i][0]+1; j++){
			for(k=0; k<2*numspecArr[leaf_coordinates[i][0]]-1; k++){
				if ( nodeScores[i][leaf_coordinates[i][0]][k] >= (maximum-Cinterval) && nodeScores[i][leaf_coordinates[i][0]][k] <= (maximum+Cinterval) ){
					//printf("Match : %d Min Root: %d Min node: %d, score: %lf\n",i,j,k,nodeScores[i][leaf_coordinates[i][0]][k]);
					voteRoot[leaf_coordinates[i][0]][k]=1;
					index++;
				}
			}
		//}
	}
	//clock_gettime(CLOCK_MONOTONIC, &tend);
	//printf("finished... %.5f\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
	minimum_score[0] = maximum;
	minimum_score[1] = forward_mismatch;
	minimum_score[2] = reverse_mismatch;
	//printf("%lf\t",minimum);
	//free(leaf_sequence);
	//free(positionsInRoot);
}
