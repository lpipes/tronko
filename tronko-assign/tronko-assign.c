#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <pthread.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "global.h"
#include "placement.h"
#include "needleman_wunsch.h"
#include "readreference.h"
#include "options.h"
#include "printAlignments.h"
#include "bwa_source_files_include.h"
#include "hashmap.h"
#include "allocateMemoryForResults.h"
#include "WFA2/wavefront_align.h"
#include "hashmap_base.h"
//int numspec, numbase, root/***seq, numundspec[MAXNUMBEROFINDINSPECIES+1]*/;
char ****taxonomyArr;
struct node **treeArr;
struct queryMatPaired *pairedQueryMat;
struct queryMatSingle *singleQueryMat;
type_of_PP Cinterval;
//struct hashmap_base map;
//struct map;
int *numspecArr, *numbaseArr, *rootArr;
struct timespec tstart={0,0}, tend={0.0};

char* strupr(char* s){
	unsigned char* t = (unsigned char*) s;
	while( *t){
		*t = toupper(*t);
		t++;
	}
	return s;
}
void store_PPs_Arr(int numberOfRoots, double c){
	int i, j, k, l;
	for(i=0; i<numberOfRoots; i++){
		for (j=0; j<2*numspecArr[i]-1; j++){
			for (k=0; k<numbaseArr[i]; k++){
				for (l=0; l<4; l++){
					if ( treeArr[i][j].posteriornc[k][l] == -1 ){   //Missing data
						treeArr[i][j].posteriornc[k][l]=1;
					}else{
							//treeArr[i][j].posteriornc[k][l] = log((c/3) + ((1-(((4*c)/3) * treeArr[i][j].posteriornc[k][l]))));
							double d = 1-c;
							double e = c/3;
							double f = d * treeArr[i][j].posteriornc[k][l];
							double g = e * (1-treeArr[i][j].posteriornc[k][l]);
							treeArr[i][j].posteriornc[k][l] = log( (f + g) );
					}
				}
			}
		}
	}
}
void assignDepthArr(int node0, int node1, int depth, int whichPartitions){
	if( node0 != -1 && node1 != -1){
		treeArr[whichPartitions][node0].depth = depth;
		treeArr[whichPartitions][node1].depth = depth;
		assignDepthArr(treeArr[whichPartitions][node0].up[0], treeArr[whichPartitions][node0].up[1],depth+1,whichPartitions);
		assignDepthArr(treeArr[whichPartitions][node1].up[0], treeArr[whichPartitions][node1].up[1],depth+1,whichPartitions);
	}
}
void printTreeInfo(int whichPartition, int node, FILE* file){
	int i,j;
	if (treeArr[whichPartition][node].up[0]==-1 && treeArr[whichPartition][node].up[1]==-1){
		fprintf(file,"%s\t%d\t%d\n",treeArr[whichPartition][node].name,whichPartition,node);
		return;
	}else{
		printTreeInfo(whichPartition,treeArr[whichPartition][node].up[0],file);
		printTreeInfo(whichPartition,treeArr[whichPartition][node].up[1],file);
	}
}
void run_bwa(int start, int end, bwaMatches* bwa_results, int concordant, int numberOfTrees, char *databasefile, int paired, int max_query_length, int max_readname_length, int max_acc_name){
	int i,j;
	int number_of_threads=1;
	if (paired != 0){
		main_mem(databasefile,end-start,number_of_threads, bwa_results, concordant, numberOfTrees, start, paired, start, end, max_query_length, max_readname_length, max_acc_name);
	}else{
		main_mem(databasefile,end-start,number_of_threads, bwa_results, concordant, numberOfTrees, start, paired, start, end, max_query_length, max_readname_length, max_acc_name);
	}
}
int getLCA_Arr(int node1, int node2, int whichRoot){
	if (node1 == node2){ return node1; }
	if (treeArr[whichRoot][node1].depth > treeArr[whichRoot][node2].depth){
		int tmp = node1;
		node1 = node2;
		node2 = tmp;
	}
	node2 = treeArr[whichRoot][node2].down;
	return getLCA_Arr(node1,node2,whichRoot);
}
int getKeysCount(int whichRoot, int node, int* minNodes, int matching_nodes, int* ancestors, int numMinNodes){
	int child0 = treeArr[whichRoot][node].up[0];
	int child1 = treeArr[whichRoot][node].up[1];
	if ( child0 != -1 && child1 != -1 ){
		matching_nodes += getKeysCount(whichRoot,child0,minNodes,matching_nodes,ancestors,numMinNodes) + getKeysCount(whichRoot,child1,minNodes,matching_nodes,ancestors,numMinNodes);
	}
	int i;
	for(i=0; i<numMinNodes; i++){
		if ( minNodes[i] == node ){
			matching_nodes++;
		}
	}
	if ( matching_nodes == numMinNodes ){
		for(i=0; i<2*numspecArr[whichRoot]-1; i++){
			if ( ancestors[i] == -1 ){
				break;
			}
		}
		ancestors[i] = node;
	}
	return matching_nodes;
}
int LCA_of_nodes(int whichRoot, int root_node, int* minNodes, int numMinNodes){
	int* ancestors = (int*)malloc((2*numspecArr[whichRoot]-1)*sizeof(int));
	int i;
	for(i=0; i<2*numspecArr[whichRoot]-1; i++){
		ancestors[i] = -1;
	}
	int matching_nodes = 0;
	getKeysCount(whichRoot, root_node, minNodes, matching_nodes, ancestors, numMinNodes);
	int LCA = ancestors[0];
	free(ancestors);
	return LCA;
}
int getLCAofArray_Arr(int *minNodes,int whichRoot,int maxNumSpec, int number_of_total_nodes){
	int LCA = minNodes[0];
	int maxDepth = 1000000000;
	int i;
	for( i=1; i<number_of_total_nodes; i++){
		if ( minNodes[i] == -1 ){ return LCA; }
		if (treeArr[whichRoot][minNodes[i]].depth < maxDepth){
			LCA = getLCA_Arr(LCA, minNodes[i], whichRoot);
		}
	}
	return LCA;
}
int getLCAofArray_Arr_Multiple(int *voteroot,int whichRoot, int maxNumSpec, int number_of_total_nodes){
	/*int LCA = minNodes[0];
	int maxDepth = 1000000000;
	int i;
	for( i=1; i<number_of_total_nodes; i++){
		if ( minNodes[i] == 0 ){ return LCA; }
		if (treeArr[whichRoot][minNodes[i]].depth < maxDepth){
			LCA = getLCA_Arr(LCA, minNodes[i], whichRoot);
		}
	}
	return LCA;*/
	int i;
	int* minNodes = (int*)malloc((2*numspecArr[whichRoot]-1)*sizeof(int));
	for(i=0; i<2*numspecArr[whichRoot]-1; i++){
		minNodes[i]=-1;
	}
	int count=0;
	for(i=0; i<2*numspecArr[whichRoot]-1; i++){
		if ( voteroot[i] == 1 ){
			minNodes[count]=i;
			count++;
		}
	}
	int LCA = LCA_of_nodes(whichRoot,rootArr[whichRoot],minNodes,count);
	free(minNodes);
	return LCA;
}
void *runAssignmentOnChunk_WithBWA(void *ptr){
	struct mystruct *mstr = (mystruct *) ptr;
	resultsStruct *results=mstr->str;
	char **rootSeqs=mstr->rootSeqs;
	int numberOfTrees = mstr->ntree;
	char query_1[mstr->max_query_length];
	char query_2[mstr->max_query_length];
	int maxNumSpec = mstr->maxNumSpec;
	int iter = 0;
	int i,j,k,lineNumber;
	int end=mstr->end;
	int *minNodes=results->minNodes;
	char **LCAnames = results->LCAnames;
	int paired = mstr->paired;
	int print_alignments = results->print_alignments;
	int use_nw = mstr->use_nw;
	int print_alignments_to_file = mstr->print_alignments_to_file;
	int print_unassigned = mstr->print_unassigned;
	int use_leaf_portion = mstr->use_leaf_portion;
	int padding = mstr->padding;
	int max_query_length = mstr->max_query_length;
	int max_readname_length = mstr->max_readname_length;
	int max_acc_name = mstr->max_acc_name;
	int max_numbase = mstr->max_numbase;
	int number_of_total_nodes = mstr->number_of_total_nodes;
	int print_all_nodes = mstr->print_all_nodes;
	/*affine_penalties_t affine_penalties = {
		.match = 0,
		.mismatch = 4,
		.gap_opening = 6,
		.gap_extension = 2,
	};*/
	/*mm_allocator_t* mm_allocator;
	affine_wavefronts_t* affine_wavefronts;
	char* pattern_alg;
	char* text_alg;*/
	/*if (use_nw==0){
		if (use_leaf_portion==1){
			mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
			pattern_alg=mm_allocator_calloc(mm_allocator,max_query_length+max_query_length+padding+50+1,char,true);
			text_alg=mm_allocator_calloc(mm_allocator,max_query_length+max_query_length+padding+50+1,char,true);
			affine_wavefronts = affine_wavefronts_new_complete(max_query_length+max_query_length+padding+50+1,max_query_length+max_query_length+padding+50+1,&affine_penalties,NULL,mm_allocator);
		}
	}*/
	char* resultsPath = (char*)malloc((max_readname_length+mstr->max_lineTaxonomy+120)*sizeof(char));
	bwaMatches* bwa_results = (bwaMatches *)malloc((end-(mstr->start))*sizeof(bwaMatches));
	for (i=0; i<end-mstr->start; i++){
		//bwa_results[i].readname = (char*)malloc((max_readname_length+1)*sizeof(char));
		bwa_results[i].concordant_matches_roots = (int *)malloc(MAX_NUM_BWA_MATCHES*sizeof(int));
		bwa_results[i].concordant_matches_nodes = (int *)malloc(MAX_NUM_BWA_MATCHES*sizeof(int));
		bwa_results[i].discordant_matches_roots = (int *)malloc(MAX_NUM_BWA_MATCHES*sizeof(int));
		bwa_results[i].discordant_matches_nodes = (int *)malloc(MAX_NUM_BWA_MATCHES*sizeof(int));
		bwa_results[i].use_portion = use_leaf_portion;
		if ( use_leaf_portion == 1 ){
			bwa_results[i].cigars_forward = (char **)malloc(MAX_NUM_BWA_MATCHES*sizeof(char *));
			bwa_results[i].starts_forward = (int *)malloc(MAX_NUM_BWA_MATCHES*sizeof(int));
			if (paired == 1){
				bwa_results[i].cigars_reverse = (char **)malloc(MAX_NUM_BWA_MATCHES*sizeof(char *));
				bwa_results[i].starts_reverse = (int *)malloc(MAX_NUM_BWA_MATCHES*sizeof(int));
			}
		}
		for(j=0; j<MAX_NUM_BWA_MATCHES; j++){
			bwa_results[i].discordant_matches_roots[j] = -1;
			bwa_results[i].discordant_matches_nodes[j] = -1;
			bwa_results[i].concordant_matches_roots[j] = -1;
			bwa_results[i].concordant_matches_nodes[j] = -1;
			if (use_leaf_portion == 1){
				bwa_results[i].starts_forward[j] = -1;
				if(paired==1){
					bwa_results[i].starts_reverse[j] = -1;
				}
			}
		}
		for(j=0; j<MAX_NUM_BWA_MATCHES; j++){
			if (use_leaf_portion == 1){
				bwa_results[i].cigars_forward[j] = (char *)malloc(MAX_CIGAR*sizeof(char));
				memset(bwa_results[i].cigars_forward[j],'\0',MAX_CIGAR);
				if (paired==1){
					bwa_results[i].cigars_reverse[j] = (char *)malloc(MAX_CIGAR*sizeof(char));
					memset(bwa_results[i].cigars_reverse[j],'\0',MAX_CIGAR);
				}
			}
		}
		bwa_results[i].n_matches = 0;
	}
	int trees_search[mstr->ntree];
	int *leaf_coord_arr;
	char *leaf_sequence;
	int *positionsInRoot;
	if (use_leaf_portion == 1){
		leaf_sequence = (char *)malloc((max_query_length+max_query_length+2*padding+1)*sizeof(char));
		positionsInRoot = (int *)malloc((max_query_length+max_query_length+2*padding+1)*sizeof(int));
	}else{
		leaf_sequence = (char *)malloc((max_query_length+mstr->max_numbase+1)*sizeof(char));
		positionsInRoot = (int *)malloc((max_query_length+mstr->max_numbase+1)*sizeof(int));
	}
	struct leafMap *leaf_map;	
	run_bwa(mstr->start, end, bwa_results, mstr->concordant, mstr->ntree, mstr->databasefile, paired, max_query_length, max_readname_length, max_acc_name);
	for ( lineNumber=mstr->start; lineNumber<end; lineNumber++){
		for(i=0; i<mstr->ntree; i++){
			trees_search[i]=-1;
		}
		j=0;
		int hashValue;
		int no_add=0;
		int leaf_iter=0;
		if (bwa_results[iter].concordant_matches_roots[0]==-1 && mstr->concordant==1){
			for (i=0; i<mstr->ntree; i++){
				if (bwa_results[iter].discordant_matches_roots[0] < 0 ){
					//for(j=0; j<mstr->ntree;j++){
					//	results->leaf_coordinates[j][0]=j;
					//	results->leaf_coordinates[j][1]=rootArr[j];
					//}
					//leaf_iter=mstr->ntree;
					//i=mstr->ntree;
					break;
				}else if ( bwa_results[iter].discordant_matches_roots[i]==-1){
					break;
				}else{
					results->leaf_coordinates[leaf_iter][0]=bwa_results[iter].discordant_matches_roots[i];
					results->leaf_coordinates[leaf_iter][1]=bwa_results[iter].discordant_matches_nodes[i];
						if (use_leaf_portion==1){
							results->starts_forward[leaf_iter] = bwa_results[iter].starts_forward[i];
							strcpy(results->cigars_forward[leaf_iter],bwa_results[iter].cigars_forward[i]);
							if ( paired==1){
								results->starts_reverse[leaf_iter] = bwa_results[iter].starts_reverse[i];
								strcpy(results->cigars_reverse[leaf_iter],bwa_results[iter].cigars_reverse[i]);
							}
						}
				}
				int index1=mstr->ntree-1;
				for(k=mstr->ntree-1; k>=0; k--){
					if (trees_search[k]==-1){
						index1=k;
					}
				}
				int found=0;
				for(k=0; k<index1; k++){
					if (trees_search[k] == results->leaf_coordinates[leaf_iter][0]){
						found=1;
					}
				}
				if (found==0){
					trees_search[index1]=results->leaf_coordinates[leaf_iter][0];
					leaf_iter++;
				}
			}
		}else if (mstr->concordant==1){
			for(i=0; i<mstr->ntree; i++){
				if (bwa_results[iter].concordant_matches_roots[i]==-1){
				//if (strlen(bwa_results[iter].concordant_leaf_matches[i])<=3){
					break;
				}else{
					results->leaf_coordinates[leaf_iter][0]=bwa_results[iter].concordant_matches_roots[i];
					results->leaf_coordinates[leaf_iter][1]=bwa_results[iter].concordant_matches_nodes[i];
					if (use_leaf_portion==1){
						results->starts_forward[leaf_iter] = bwa_results[iter].starts_forward[i];
						strcpy(results->cigars_forward[leaf_iter],bwa_results[iter].cigars_forward[i]);
						if ( paired==1){
							results->starts_reverse[leaf_iter] = bwa_results[iter].starts_reverse[i];
							strcpy(results->cigars_reverse[leaf_iter],bwa_results[iter].cigars_reverse[i]);
						}
					}
					int index1=mstr->ntree-1;
					for(k=mstr->ntree-1; k>=0; k--){
						if (trees_search[k]==-1){
							index1=k;
						}
					}
					int found=0;
					for(k=0; k<index1; k++){
						if (trees_search[k] == results->leaf_coordinates[leaf_iter][0]){
							found=1;
						}
					}
					if (found==0){
						trees_search[index1]=results->leaf_coordinates[leaf_iter][0];
						leaf_iter++;
					}
				}
			}
		}else{
			j=0;
			for(i=0; i<mstr->ntree; i++){
				//if(strlen(bwa_results[iter].discordant_leaf_matches[i])<=3){
				if(bwa_results[iter].discordant_matches_roots[i]==-1){
					break;
				}else{
					if (no_add==0){
						results->leaf_coordinates[leaf_iter][0]=bwa_results[iter].discordant_matches_roots[i];
						results->leaf_coordinates[leaf_iter][1]=bwa_results[iter].discordant_matches_nodes[i];
						if (use_leaf_portion==1){
							results->starts_forward[leaf_iter] = bwa_results[iter].starts_forward[i];
							strcpy(results->cigars_forward[leaf_iter],bwa_results[iter].cigars_forward[i]);
							if (paired==1){
								results->starts_reverse[leaf_iter] = bwa_results[iter].starts_reverse[i];
								strcpy(results->cigars_reverse[leaf_iter],bwa_results[iter].cigars_reverse[i]);
							}
						}
						int index1=mstr->ntree-1;
						for(k=mstr->ntree-1; k>=0; k--){
							if (trees_search[k]==-1){
								index1=k;
							}
						}
						int found=0;
						for(k=0; k<index1; k++){
							if (trees_search[k] == results->leaf_coordinates[leaf_iter][0]){
								found=1;
							}
						}
						if (found==0){
							trees_search[index1]=results->leaf_coordinates[leaf_iter][0];
							leaf_iter++;
						}
						j++;
					}
					no_add=0;
				}
			}
			for(i=0; i<mstr->ntree; i++){
				if (bwa_results[iter].concordant_matches_roots[i]==-1){
				//if (strlen(bwa_results[iter].concordant_leaf_matches[i])<=3){
					break;
				}else{
					if (no_add==0){
						results->leaf_coordinates[leaf_iter][0]=bwa_results[iter].discordant_matches_roots[i];
						results->leaf_coordinates[leaf_iter][1]=bwa_results[iter].discordant_matches_nodes[i];
						if (use_leaf_portion == 1){
							results->starts_forward[leaf_iter] = bwa_results[iter].starts_forward[i];
							strcpy(results->cigars_forward[leaf_iter],bwa_results[iter].cigars_forward[i]);
							if (paired==1){
								results->starts_reverse[leaf_iter] = bwa_results[iter].starts_reverse[i];
								strcpy(results->cigars_reverse[leaf_iter],bwa_results[iter].cigars_reverse[i]);
							}
						}
					int index1=mstr->ntree-1;
					for(k=mstr->ntree-1; k>=0; k--){
						if (trees_search[k]=-1){
							index1=k;
						}
					}
					int found=0;
					for(k=0; k<index1; k++){
						if (trees_search[k] == results->leaf_coordinates[leaf_iter][0]){
							found=1;
						}
					}
					if (found==0){
						trees_search[index1]=results->leaf_coordinates[leaf_iter][0];
						leaf_iter++;
					}
						j++;
					}
					no_add=0;
				}
			}
		}
		results->minimum[0]=0;
		if (leaf_iter > 0 ){
			//strcmp(results->cigars_forward[0],"*")!=0 && strcmp(results->cigars_forward[0]," ")!=0){
		if (paired != 0){
			if ( use_nw==0 ){
				place_paired(pairedQueryMat->query1Mat[lineNumber],pairedQueryMat->query2Mat[lineNumber],rootSeqs,mstr->ntree,results->positions,results->locQuery,results->nodeScores,results->voteRoot, leaf_iter, results->leaf_coordinates,paired,results->minimum,mstr->alignmentsdir,pairedQueryMat->forward_name[lineNumber],pairedQueryMat->reverse_name[lineNumber],print_alignments,leaf_sequence,positionsInRoot,maxNumSpec,results->starts_forward,results->cigars_forward,results->starts_reverse,results->cigars_reverse,print_alignments_to_file,use_leaf_portion,padding,max_query_length,max_numbase,print_all_nodes);
			}else{
				place_paired_with_nw(pairedQueryMat->query1Mat[lineNumber],pairedQueryMat->query2Mat[lineNumber],rootSeqs,mstr->ntree,results->positions,results->locQuery,results->nw,results->aln,results->scoring,results->nodeScores,results->voteRoot, leaf_iter, results->leaf_coordinates,paired,results->minimum,mstr->alignmentsdir,pairedQueryMat->forward_name[lineNumber],pairedQueryMat->reverse_name[lineNumber],print_alignments,leaf_sequence,positionsInRoot,maxNumSpec,results->starts_forward,results->cigars_forward,results->starts_reverse,results->cigars_reverse,print_alignments_to_file,use_leaf_portion,padding,max_query_length,max_numbase,print_all_nodes);
			}
		}else{
			if (use_nw==0){
				place_paired(singleQueryMat->queryMat[lineNumber],NULL,rootSeqs,mstr->ntree,results->positions,results->locQuery,results->nodeScores,results->voteRoot, leaf_iter, results->leaf_coordinates,paired,results->minimum,mstr->alignmentsdir,singleQueryMat->name[lineNumber],NULL,print_alignments,leaf_sequence,positionsInRoot,maxNumSpec,results->starts_forward,results->cigars_forward,results->starts_reverse,results->cigars_reverse,print_alignments_to_file,use_leaf_portion,padding,max_query_length,max_numbase,print_all_nodes);
			}else{
				place_paired_with_nw(singleQueryMat->queryMat[lineNumber],NULL,rootSeqs,mstr->ntree,results->positions,results->locQuery,results->nw,results->aln,results->scoring,results->nodeScores,results->voteRoot,leaf_iter, results->leaf_coordinates,paired,results->minimum,mstr->alignmentsdir,singleQueryMat->name[lineNumber],NULL,print_alignments,leaf_sequence,positionsInRoot,maxNumSpec,results->starts_forward,results->cigars_forward,results->starts_reverse,results->cigars_reverse,print_alignments_to_file,use_leaf_portion,padding,max_query_length,max_numbase,print_all_nodes);
			}
		}
		numberOfTrees = leaf_iter;
		for(i=0;i<number_of_total_nodes;i++){
			minNodes[i]=-1;
		}
		for(i=0; i<leaf_iter; i++){
			for(j=0; j<2*numspecArr[trees_search[i]]-1; j++){
				results->nodeScores[i][trees_search[i]][j]=0;
			}
			results->leaf_coordinates[i][0]=-1;
			results->leaf_coordinates[i][1]=-1;
		}
			if (use_leaf_portion==1){
				results->starts_forward[leaf_iter] = -1;
				for(j=0; j<MAX_CIGAR; j++){
					results->cigars_forward[leaf_iter][j]='\0';
				}
				if(paired==1){
					results->starts_reverse[leaf_iter]=-1;
					for(j=0; j<MAX_CIGAR; j++){
						results->cigars_reverse[leaf_iter][j]='\0';
					}
				}
			}
		}
		int countVotes[mstr->ntree];
		int count=0;
		for(i=0; i<mstr->ntree; i++){
			countVotes[i]=0;
			for(j=0;j<2*numspecArr[i]-1;j++){
				if (results->voteRoot[i][j]==1){
					countVotes[i]++;
					minNodes[count]=j;
					count++;
				}
			}
		}
		int numMinNodes=count;
		int max=0;
		int maxRoot=-1;
		count=0;
		for(i=0;i<mstr->ntree;i++){
			if (countVotes[i]>max){
				max=countVotes[i];
				maxRoot=i;
			}
			if (countVotes[i]>0){
				count++;
			}
		}
		int LCA, LCAs[count], maxRoots[count];
		int count2=0;
		for(i=0;i<mstr->ntree;i++){
			if ( countVotes[i]>0){
				maxRoots[count2]=i;
				count2++;
			}
		}
		int unassigned=0;
		int minLevel=0;
		int taxRoot,taxIndex0,taxIndex1,taxNode;
		if ( count == 1 ){
			clock_gettime(CLOCK_MONOTONIC, &tstart);
			//LCA=getLCAofArray_Arr(minNodes,maxRoot,maxNumSpec,number_of_total_nodes);
			LCA = LCA_of_nodes(maxRoot,rootArr[maxRoot],minNodes,numMinNodes);
		}else if (count != 0){
			for(i=0;i<count;i++){
				for(j=0; j<mstr->max_lineTaxonomy; j++){
					LCAnames[i][j]='\0';
				}
			}
			for(i=0;i<count;i++){
				LCAs[i]=getLCAofArray_Arr_Multiple(results->voteRoot[maxRoots[i]],maxRoots[i],maxNumSpec,number_of_total_nodes);
				if ( treeArr[maxRoots[i]][LCAs[i]].taxIndex[1]!=-1){ 
					strcpy(LCAnames[i],taxonomyArr[maxRoots[i]][treeArr[maxRoots[i]][LCAs[i]].taxIndex[0]][treeArr[maxRoots[i]][LCAs[i]].taxIndex[1]]);
					if ( treeArr[maxRoots[i]][LCAs[i]].taxIndex[1] > minLevel ){
						minLevel=treeArr[maxRoots[i]][LCAs[i]].taxIndex[1];
					}
				}else{
					unassigned=1;
				}
			}
			LCA = LCAs[0];
			int correctTax=0;
			int stop=0;
			while(stop==0 && minLevel<=6){
			for(i=0; i<count;i++){
				if (treeArr[maxRoots[i]][LCAs[i]].taxIndex[0] != -1){
				for(j=i+1; j<count; j++){
					if ( treeArr[maxRoots[j]][LCAs[j]].taxIndex[0] != -1){
					if ( strcmp(taxonomyArr[maxRoots[i]][treeArr[maxRoots[i]][LCAs[i]].taxIndex[0]][minLevel],taxonomyArr[maxRoots[j]][treeArr[maxRoots[j]][LCAs[j]].taxIndex[0]][minLevel])==0 && taxonomyArr[maxRoots[i]][treeArr[maxRoots[i]][LCAs[i]].taxIndex[0]][minLevel] != "NA" ){
						correctTax++;
					}
					}
				}
				}
			}
			if (correctTax>=count-1){
				stop=1;
				taxRoot=maxRoots[0];
				taxNode=LCAs[0];
				taxIndex0=treeArr[maxRoots[0]][LCAs[0]].taxIndex[0];
				taxIndex1=minLevel;
			}else{
				minLevel++;
			}
			}
			if (minLevel>6){ unassigned=1;}
		}
		for(i=0; i<max_readname_length+mstr->max_lineTaxonomy+120;i++){
			resultsPath[i] = '\0';
		}
		int print_un=1;
		if (count==1){
			if ( paired != 0 ){
				strcpy(resultsPath,pairedQueryMat->forward_name[lineNumber]);
			}else{
				strcpy(resultsPath,singleQueryMat->name[lineNumber]);
			}
			strcat(resultsPath,"\t");
			int taxIndex1 = treeArr[maxRoot][LCA].taxIndex[1];
			if ( taxIndex1 == -1 ){
				strcat(resultsPath,"unassigned Euk or Bac");
				//if (print_unassigned==0){
				//	print_un = 0;
				//}
			}else{
			for(i=6;i>=taxIndex1;i--){
				if (i==taxIndex1){
					strcat(resultsPath,taxonomyArr[maxRoot][treeArr[maxRoot][LCA].taxIndex[0]][i]);
				}else{
					strcat(resultsPath,taxonomyArr[maxRoot][treeArr[maxRoot][LCA].taxIndex[0]][i]);
					strcat(resultsPath,";");
				}
			}
			}
			strcat(resultsPath,"\t");
			char *num = NULL;
			asprintf(&num,"%lf",results->minimum[0]);
			strcat(resultsPath,num);
			strcat(resultsPath,"\t");
			//free(num);
			char *num2 = NULL;
			asprintf(&num2,"%lf",results->minimum[1]);
			strcat(resultsPath,num2);
			strcat(resultsPath,"\t");
			//free(num2);
			char *num3 = NULL;
			asprintf(&num3,"%lf",results->minimum[2]);
			strcat(resultsPath,num3);
			strcat(resultsPath,"\t");
			//free(num3);
			char *num4 = NULL;
			asprintf(&num4,"%d",maxRoot);
			strcat(resultsPath,num4);
			strcat(resultsPath,"\t");
			//free(num4);
			char *num5 = NULL;
			asprintf(&num5,"%d",LCA);
			strcat(resultsPath,num5);
			//free(num5);
			//char* appendScores = (char*)malloc(30*sizeof(char));
			//sprintf(appendScores,"%lf\t%lf\t%lf\t%d\t%d",results->minimum[0],results->minimum[1],results->minimum[2],maxRoot,LCA);
			//strcat(resultsPath,appendScores);
			//free(appendScores);
			//printf("%s\n",resultsPath);
			//if ( print_un == 1){
				strcpy(results->taxonPath[iter], resultsPath);
			//}
			free(num);
			free(num2);
			free(num3);
			free(num4);
			free(num5);
		}else if (count==0 /*&& print_unassigned==1*/){
			if (paired != 0){
				strcpy(resultsPath,pairedQueryMat->forward_name[lineNumber]);
			}else{
				strcpy(resultsPath,singleQueryMat->name[lineNumber]);
			}
			strcat(resultsPath,"\tunassigned");
			strcpy(results->taxonPath[iter],resultsPath);
		}else{
			if (paired != 0){
				strcpy(resultsPath,pairedQueryMat->forward_name[lineNumber]);
			}else{
				strcpy(resultsPath,singleQueryMat->name[lineNumber]);
			}
			strcat(resultsPath,"\t");
			if (unassigned==1 ){
				strcat(resultsPath, "unassigned Eukaryote or Bacteria");
				//if (print_unassigned==0){
				//	print_un=0;
				//}
			}else{
			int taxIndex1=minLevel;
			for(i=6;i>=taxIndex1;i--){
				if (i==taxIndex1){
					strcat(resultsPath,taxonomyArr[taxRoot][treeArr[taxRoot][taxNode].taxIndex[0]][i]);
				}else{
					strcat(resultsPath,taxonomyArr[taxRoot][treeArr[taxRoot][taxNode].taxIndex[0]][i]);
					strcat(resultsPath,";");
				}
			}
			}
			strcat(resultsPath,"\t");
			//char* appendScores = (char*)malloc(18*sizeof(char));
			//char *appendScores = NULL;
			//asprintf(&appendScores,"%lf\t%d\t%d\t%d\t%d",results->minimum[0],results->minimum[1],results->minimum[2],maxRoot,LCA);
			char *num = NULL;
			asprintf(&num,"%lf",results->minimum[0]);
			strcat(resultsPath,num);
			strcat(resultsPath,"\t");
			//free(num);
			char *num2 = NULL;
			asprintf(&num2,"%lf",results->minimum[1]);
			strcat(resultsPath,num2);
			strcat(resultsPath,"\t");
			//free(num2);
			char *num3 = NULL;
			asprintf(&num3,"%lf",results->minimum[2]);
			strcat(resultsPath,num3);
			strcat(resultsPath,"\t");
			//free(num3);
			char *num4 = NULL;
			asprintf(&num4,"%d",maxRoot);
			strcat(resultsPath,num4);
			strcat(resultsPath,"\t");
			//free(num4);
			char *num5 = NULL;
			asprintf(&num5,"%d",LCA);
			strcat(resultsPath,num5);
			//strcat(resultsPath,appendScores);
			//if (print_un==1){
				strcpy(results->taxonPath[iter],resultsPath);
			//}
			//free(appendScores);
		}
		//free(bwa_results[iter].concordant_matches);
		//free(bwa_results[iter].discordant_matches);
		if (use_leaf_portion == 1){
			free(bwa_results[iter].starts_forward);
			if(paired==1){
				free(bwa_results[iter].starts_reverse);
			}
		}
		for(i=0; i<MAX_NUM_BWA_MATCHES;i++){
			//free(bwa_results[iter].discordant_leaf_matches);
			//free(bwa_results[iter].concordant_leaf_matches);
			if (use_leaf_portion == 1){
				free(bwa_results[iter].cigars_forward[i]);
				if(paired==1){
						free(bwa_results[iter].cigars_reverse[i]);
				}
			}
		}
		//free(bwa_results[iter].discordant_leaf_matches);
		//free(bwa_results[iter].concordant_leaf_matches);
		if (use_leaf_portion == 1){
			free(bwa_results[iter].cigars_forward);
			if (paired==1){
				free(bwa_results[iter].cigars_reverse);
			}
		}
		free(bwa_results[iter].concordant_matches_roots);
		free(bwa_results[iter].concordant_matches_nodes);
		free(bwa_results[iter].discordant_matches_roots);
		free(bwa_results[iter].discordant_matches_nodes);
		//free(bwa_results[iter].readname);
		iter++;
		results->minimum[0] = -1;
		results->minimum[1] = -1;
		results->minimum[2] = -1;
		LCA = -1;
		if (leaf_iter > 0){
			for(i=0; i<mstr->ntree; i++){
				for (j=0; j<2*numspecArr[i]-1; j++){
					results->voteRoot[i][j]=0;
				}
			}
		}
	}
	/*if (use_nw == 0){
		affine_wavefronts_delete(affine_wavefronts);
		mm_allocator_free(mm_allocator,pattern_alg);
		mm_allocator_free(mm_allocator,text_alg);
		mm_allocator_delete(mm_allocator);
	}*/
	/*if (use_nw == 0){
		mm_allocator_delete(mm_allocator);
	}*/
	free(leaf_sequence);
	free(positionsInRoot);
	free(bwa_results);
	free(resultsPath);
	pthread_exit(NULL);
}
int main(int argc, char **argv){
	int i, j, k, numberOfTrees;
	Options opt;
	opt.use_nw=0;
	opt.reference_mode=1; // default use reference
	opt.reverse_single_read=0; // default no reverse single read
	opt.reverse_second_of_paired_read=0;
	opt.print_alignments = 0;
	opt.print_node_info[0] = '\0';
	opt.results_file[0] = '\0';
	opt.reference_file[0] = '\0';
	opt.print_trees_dir[0] = '\0';
	opt.fastq=0; //default is FASTA
	opt.unassigned=0; //don't print unassigned sequences
	opt.print_alignments_to_file=0; //don't print alignments to file
	opt.use_leaf_portion=0;
	opt.padding=0;
	opt.skip_build=0;
	opt.number_of_cores=1;
	opt.number_of_lines_to_read=50000;
	opt.score_constant = 0.01;
	opt.print_all_nodes=0;
	parse_options(argc, argv, &opt);
	struct stat st = {0};
	if ( stat(opt.reference_file, &st) == -1 ){
		printf("Cannot find reference_tree.txt file. Exiting...\n");
		exit(-1);
	}
	if ( opt.fastq == 1 && opt.number_of_lines_to_read%4 != 0 ){
		printf("You chose FASTQ for your queries but the number of lines to read are not divisible by 4. Change -L to be divisible by 4. Exiting...\n");
	       exit(-1);	
	}
	if ( opt.fastq == 0 && opt.number_of_lines_to_read%2 != 0 ){
		printf("You chose FASTA for your queries but the number of lines to read are not divisible by 2. Change -L to be divisible by 2. Exiting...\n");
		exit(-1);
	}
	gzFile referenceTree = Z_NULL;
	referenceTree = gzopen(opt.reference_file,"r");
	assert(Z_NULL!=referenceTree);
	int* name_specs = (int*)malloc(3*sizeof(int));
	name_specs[0]=0;
	name_specs[1]=0;
	name_specs[2]=0;
	numberOfTrees = readReferenceTree(referenceTree,name_specs);
	gzclose(referenceTree);
	int max_nodename = name_specs[0];
	int max_taxname = name_specs[1];
	int max_lineTaxonomy = name_specs[2];
	free(name_specs);
	if ( opt.print_node_info[0] != '\0' ){
		printf("Printing Accession IDs, Tree numbers, and leaf numbers...\n");
		FILE* tree_info = fopen(opt.print_node_info,"w");
		if (tree_info == NULL ){ printf("Error opening node info file!\n"); exit(1); }
		for(i=0; i<numberOfTrees; i++){
			printTreeInfo(i,rootArr[i],tree_info);
		}
		fclose(tree_info);
		exit(1);
	}
	/*if (opt.print_trees_dir[0] = '\0' ){
		printf("Printing Newick trees used for assignment...\n");
		char newickbuffer[3000];
		for(i=0; i<numberOfTrees; i++){
			FILE *newick_out;
			snprintf(buffer,3000,"%s/%s.nwk",opt.print_trees_dir,i);
			if (NULL==(newick_out=fopen(buffer,"w"))){ puts("Cannot open newick file!\n"); exit(-1); }
			printNewick(newick_out,i,rootArr[i]);
			fclose(newick_out);
		}
		exit(1);
	}*/
	opt.print_leave_seqs=0;
	if (opt.print_leave_seqs == 1){
		printf("Printingleavesfile\n");
		FILE* leaves_file = fopen("leaves.fasta","w");
		if(leaves_file == NULL ){ printf("Error opening leaves file!\n"); exit(1); }
		for(i=0; i<numberOfTrees; i++){
			printLeaveSeqsToFile(leaves_file,rootArr[i],i,numbaseArr[i]);	
		}
		exit(1);
	}
	store_PPs_Arr(numberOfTrees,opt.score_constant);
	int maxNumSpec=0;
	int maxNumBase=0;
	int *specs = (int*)malloc(2*sizeof(int));
	specs[0]=0;
	specs[1]=0;
	int numspec_total=setNumbase_setNumspec(numberOfTrees,specs);
	maxNumSpec = specs[0];
	maxNumBase = specs[1];
	free(specs);
	Cinterval = opt.cinterval;
	//HASHMAP(char, leafMap) map;
	//hashmap_init(&map, hashmap_hash_string, strcmp);
	//for(i=0; i<numberOfTrees; i++){
	//	for(j=numspecArr[i]-1; j<2*numspecArr[i]-1; j++){
	//		struct leafMap *l;
	//		l = malloc(sizeof(*l));
	//		l->name = treeArr[i][j].name;
	//		l->root = i;
	//		l->node = j;
	//		hashmap_put(&map,l->name,l);
	//	}
	//}
	//leaf_map = (leafMap *)malloc(numspec_total*sizeof(leafMap));
	//k=0;
	//for(i=0; i<numberOfTrees; i++){
	//	for(j=numspecArr[i]-1; j<2*numspecArr[i]-1; j++){
	//		leaf_map[k].name = (char*)malloc(max_nodename*sizeof(char));
	//		strcpy(leaf_map[k].name,treeArr[i][j].name);
	//		leaf_map[k].root = i;
	//		leaf_map[k].node = j;
	//		k++;
	//	}
	//}
	int *read_specs = (int*)malloc(2*sizeof(int));
	read_specs[0] = 0;
	read_specs[1] = 0;
	int number_of_total_nodes = 0;
	for(i=0; i<numberOfTrees; i++){
		number_of_total_nodes += 2*numspecArr[i]-1;
	}
	int max_name_length = 0;
	int max_query_length = 0;
	int numberOfLinesToRead=opt.number_of_lines_to_read;
	mystruct mstr[opt.number_of_cores];//array of stuct that contains input and output for each thread
	if ( strcmp("single",opt.paired_or_single)==0){
		if (opt.skip_build==0){
			bwa_index(2,opt.fasta_file);
		}
		gzFile *reads_file =gzopen(opt.read1_file,"r");
		if ( reads_file == (gzFile) Z_NULL ){
			printf("**reads file could not be opened.\n");
		}
		find_specs_for_reads(read_specs,reads_file,opt.fastq);
		gzclose(reads_file);
		max_name_length = read_specs[0];
		max_query_length = read_specs[1];
		free(read_specs);
		singleQueryMat = malloc(sizeof(struct queryMatSingle));
		if ( opt.fastq == 0 ){
			singleQueryMat->queryMat = (char **)malloc(sizeof(char *)*(numberOfLinesToRead/2));
			singleQueryMat->name = (char **)malloc(sizeof(char *)*(numberOfLinesToRead/2));
			for (i=0; i<numberOfLinesToRead/2; i++){
				singleQueryMat->queryMat[i] = (char *)malloc(sizeof(char)*(max_query_length+1));
				singleQueryMat->name[i] = (char *)malloc(sizeof(char)*(max_name_length+1));
			}
		}else{
			singleQueryMat->queryMat = (char **)malloc(sizeof(char *)*(numberOfLinesToRead/4));
			singleQueryMat->name = (char **)malloc(sizeof(char *)*(numberOfLinesToRead/4));
			for (i=0; i<numberOfLinesToRead/4; i++){
				singleQueryMat->queryMat[i] = (char *)malloc(sizeof(char)*(max_query_length+1));
				singleQueryMat->name[i] = (char *)malloc(sizeof(char)*(max_name_length+1));
			}
		}
		FILE *results = fopen(opt.results_file,"w");
		if ( results == NULL ){ printf("Error opening output file!\n"); exit(1); }
		fprintf(results,"Readname\tTaxonomic_Path\tScore\tForward_Mismatch\tReverse_Mismatch\tTree_Number\tNode_Number\n");	
		int keepTrackOfReadLine=0;
		pthread_t threads[opt.number_of_cores];//array of our threads
		int divideFile, start, end;
		int returnLineNumber=0;
		gzFile *seqinfile = gzopen(opt.read1_file,"r");
		if (seqinfile == (gzFile) Z_NULL){
			printf("*** fasta file could not be opened.\n");
		}
		int first_iter=1;
		for(i=0; i<opt.number_of_cores; i++){
			mstr[i].str = malloc(sizeof(struct resultsStruct));
			if ( opt.fastq == 0){
				allocateMemForResults(mstr[i].str, numberOfLinesToRead/2, opt.number_of_cores, numberOfTrees, opt.print_alignments,maxNumSpec,0,opt.use_nw,max_lineTaxonomy,max_name_length,max_query_length,maxNumBase,opt.use_leaf_portion,opt.padding,number_of_total_nodes);
			}else{
				allocateMemForResults(mstr[i].str, numberOfLinesToRead/4, opt.number_of_cores, numberOfTrees, opt.print_alignments,maxNumSpec,0,opt.use_nw,max_lineTaxonomy,max_name_length,max_query_length,maxNumBase,opt.use_leaf_portion,opt.padding,number_of_total_nodes);
			}
			mstr[i].concordant=0;
			mstr[i].maxNumSpec=maxNumSpec;
			mstr[i].databasefile=opt.fasta_file;
			mstr[i].numspec_total=numspec_total;
			mstr[i].use_nw = opt.use_nw;
			mstr[i].print_alignments_to_file = opt.print_alignments_to_file;
			mstr[i].print_unassigned = opt.unassigned;
			mstr[i].use_leaf_portion = opt.use_leaf_portion;
			mstr[i].padding = opt.padding;
			mstr[i].max_query_length = max_query_length;
			mstr[i].max_readname_length = max_name_length;
			mstr[i].max_acc_name = max_nodename;
			mstr[i].max_numbase = maxNumBase;
			mstr[i].max_lineTaxonomy = max_lineTaxonomy;
			mstr[i].number_of_total_nodes = number_of_total_nodes;
			mstr[i].print_all_nodes = opt.print_all_nodes;
		}
		while (1){
			if (opt.fastq==0){
				returnLineNumber=readInXNumberOfLines(numberOfLinesToRead/2,seqinfile,0,opt,max_query_length,max_name_length);
			}else{
				returnLineNumber=readInXNumberOfLines_fastq(numberOfLinesToRead/4,seqinfile,0,opt,max_query_length,max_name_length,first_iter);
			}
			if (returnLineNumber==0){
				break;
			}
			divideFile = returnLineNumber/opt.number_of_cores;
			first_iter=0;
			j=0;
			for (i=0; i<opt.number_of_cores; i++){
				start=j;
				end=j+divideFile;
				if ( i==opt.number_of_cores-1){
					end=returnLineNumber;
				}
				mstr[i].start=start;
				mstr[i].end=end;
				mstr[i].paired=0;
				mstr[i].ntree = numberOfTrees;
				mstr[i].alignmentsdir = opt.print_alignments_dir;
				j=j+divideFile;
				mstr[i].str->taxonPath =(char**) malloc((end-start)*(sizeof(char *)));
				for(k=0; k<end-start; k++){
					mstr[i].str->taxonPath[k] = malloc((max_name_length+max_lineTaxonomy+120)*(sizeof(char)));
				}
			}
			for(i=0; i<opt.number_of_cores;i++){
				pthread_create(&threads[i], NULL, runAssignmentOnChunk_WithBWA, &mstr[i]);
			}
			for ( i=0; i<opt.number_of_cores;i++){
				pthread_join(threads[i], NULL);
			}
			for ( i=0; i<opt.number_of_cores; i++){
				for ( j=0; j<(mstr[i].end-mstr[i].start); j++){
					fprintf(results,"%s\n",mstr[i].str->taxonPath[j]);
				}
			}
			for ( i=0; i<opt.number_of_cores; i++){
				for(j=0; j<mstr[i].end-mstr[i].start; j++){
					free(mstr[i].str->taxonPath[j]);
				}
				free(mstr[i].str->taxonPath);
			}
		}
		fclose(results);
		gzclose(seqinfile);
		if ( opt.fastq == 0 ){
			for(i=0; i<numberOfLinesToRead/2; i++){
				free(singleQueryMat->queryMat[i]);
				free(singleQueryMat->name[i]);
			}
		}else{
			for(i=0; i<numberOfLinesToRead/4; i++){
				free(singleQueryMat->queryMat[i]);
				free(singleQueryMat->name[i]);
			}
		}
		free(singleQueryMat->queryMat);
		free(singleQueryMat->name);
		free(singleQueryMat);
	}else{
		gzFile *seqinfile_1 = gzopen(opt.read1_file,"r");
		gzFile *seqinfile_2 = gzopen(opt.read2_file,"r");
		if (seqinfile_1 == (gzFile) Z_NULL){
			printf("*** fasta/fastq file could not be opened.\n");
		}
		if (seqinfile_2 == (gzFile) Z_NULL){
			printf("*** fasta/fastq file could not be opened.\n");
		}
		find_specs_for_reads(read_specs,seqinfile_1,opt.fastq);
		gzclose(seqinfile_1);
		find_specs_for_reads(read_specs,seqinfile_2,opt.fastq);
		gzclose(seqinfile_2);
		max_name_length = read_specs[0];
		max_query_length = read_specs[1];
		free(read_specs);
		if (opt.skip_build==0){
			bwa_index(2,opt.fasta_file);
		}
		pairedQueryMat = malloc(sizeof(struct queryMatPaired));
		if (opt.fastq==0){
			pairedQueryMat->query1Mat = (char **)malloc(sizeof(char *)*numberOfLinesToRead/2);
			pairedQueryMat->query2Mat = (char **)malloc(sizeof(char *)*numberOfLinesToRead/2);
			pairedQueryMat->forward_name = (char **) malloc(sizeof(char *)*numberOfLinesToRead/2);
			pairedQueryMat->reverse_name = (char **) malloc(sizeof(char *)*numberOfLinesToRead/2);
			for ( i=0; i<numberOfLinesToRead/2; i++){
				pairedQueryMat->query1Mat[i] = (char *)malloc(sizeof(char)*max_query_length+1);
				pairedQueryMat->query2Mat[i] = (char *)malloc(sizeof(char)*max_query_length+1);
				pairedQueryMat->forward_name[i] = (char *)malloc(sizeof(char)*max_name_length+1);
				pairedQueryMat->reverse_name[i] = (char *)malloc(sizeof(char)*max_name_length+1);
				for(j=0; j<max_name_length+1; j++){
					pairedQueryMat->forward_name[i][j]='\0';
					pairedQueryMat->reverse_name[i][j]='\0';
				}
				for(j=0; j<max_query_length+1; j++){
					pairedQueryMat->query1Mat[i][j] = '\0';
					pairedQueryMat->query2Mat[i][j] = '\0';
				}
			}
		}else{
			pairedQueryMat->query1Mat = (char **)malloc(sizeof(char *)*numberOfLinesToRead/4);
			pairedQueryMat->query2Mat = (char **)malloc(sizeof(char *)*numberOfLinesToRead/4);
			pairedQueryMat->forward_name = (char **) malloc(sizeof(char *)*numberOfLinesToRead/4);
			pairedQueryMat->reverse_name = (char **) malloc(sizeof(char *)*numberOfLinesToRead/4);
			for ( i=0; i<numberOfLinesToRead/4; i++){
				pairedQueryMat->query1Mat[i] = (char *)malloc(sizeof(char)*max_query_length+1);
				pairedQueryMat->query2Mat[i] = (char *)malloc(sizeof(char)*max_query_length+1);
				pairedQueryMat->forward_name[i] = (char *)malloc(sizeof(char)*max_name_length+1);
				 pairedQueryMat->reverse_name[i] = (char *)malloc(sizeof(char)*max_name_length+1);
				for(j=0; j<max_name_length+1; j++){
					pairedQueryMat->forward_name[i][j]='\0';
					pairedQueryMat->reverse_name[i][j]='\0';
				}
				for(j=0; j<max_query_length+1; j++){
					pairedQueryMat->query1Mat[i][j]='\0';
					pairedQueryMat->query2Mat[i][j]='\0';
				}
			}
		}
		FILE *results = fopen(opt.results_file,"w");
		if ( results == NULL ){ printf("Error opening output file!\n"); exit(1); }
		fprintf(results,"Readname\tTaxonomic_Path\tScore\tForward_Mismatch\tReverse_Mismatch\tTree_Number\tNode_Number\n");	
		int keepTrackOfReadLine=0;
		pthread_t threads[opt.number_of_cores];//array of our thrads
		int divideFile, start, end;
		int returnLineNumber=0;//<- this is the number of records that we have read.
		int returnLineNumber2 =0;
		int concordant=1;
		seqinfile_1 = gzopen(opt.read1_file,"r");
		seqinfile_2 = gzopen(opt.read2_file,"r");
		if (seqinfile_1 == (gzFile) Z_NULL){
			printf("*** fasta/fastq file could not be opened.\n");
		}
		if (seqinfile_2 == (gzFile) Z_NULL){
			printf("*** fasta/fastq file could not be opened.\n");
		}
		int first_iter=1;
		for(i=0;i<opt.number_of_cores;i++){
			mstr[i].str = malloc(sizeof(struct resultsStruct));
			if ( opt.fastq == 0){
				allocateMemForResults(mstr[i].str, numberOfLinesToRead/2, opt.number_of_cores, numberOfTrees, opt.print_alignments,maxNumSpec,1,opt.use_nw,max_lineTaxonomy,max_name_length,max_query_length,maxNumBase,opt.use_leaf_portion,opt.padding,number_of_total_nodes);
			}else{
				allocateMemForResults(mstr[i].str, numberOfLinesToRead/4, opt.number_of_cores, numberOfTrees, opt.print_alignments,maxNumSpec,1,opt.use_nw,max_lineTaxonomy,max_name_length,max_query_length,maxNumBase,opt.use_leaf_portion,opt.padding,number_of_total_nodes);
			}
			mstr[i].concordant=concordant;
			mstr[i].maxNumSpec=maxNumSpec;
			mstr[i].databasefile=opt.fasta_file;
			mstr[i].numspec_total=numspec_total;
			mstr[i].use_nw = opt.use_nw;
			mstr[i].print_alignments_to_file = opt.print_alignments_to_file;
			mstr[i].print_unassigned = opt.unassigned;
			mstr[i].use_leaf_portion = opt.use_leaf_portion;
			mstr[i].padding = opt.padding;
			mstr[i].max_query_length = max_query_length;
			mstr[i].max_readname_length = max_name_length;
			mstr[i].max_acc_name = max_nodename;
			mstr[i].max_numbase = maxNumBase;
			mstr[i].max_lineTaxonomy = max_lineTaxonomy;
			mstr[i].number_of_total_nodes = number_of_total_nodes;
			mstr[i].print_all_nodes = opt.print_all_nodes;
		}
		while (1){
			if (opt.fastq==0){
				returnLineNumber=readInXNumberOfLines(numberOfLinesToRead/2,seqinfile_1,1,opt,max_query_length,max_name_length);
			}else{
				returnLineNumber=readInXNumberOfLines_fastq(numberOfLinesToRead/4,seqinfile_1,1,opt,max_query_length,max_name_length,first_iter);
			}
			if (returnLineNumber==0)
				break;
			if (opt.fastq==0){
				returnLineNumber2 = readInXNumberOfLines ( numberOfLinesToRead/2, seqinfile_2, 2, opt,max_query_length,max_name_length);
			}else{
				returnLineNumber2 = readInXNumberOfLines_fastq(numberOfLinesToRead/4,seqinfile_2, 2, opt,max_query_length,max_name_length,first_iter);
			}
			returnLineNumber = returnLineNumber2;
			first_iter=0;
			divideFile= returnLineNumber/opt.number_of_cores;
			j=0;
			for ( i=0; i<opt.number_of_cores; i++){
				start=j;
				end=j+divideFile;
				if(i==opt.number_of_cores-1)
					end=returnLineNumber;
				mstr[i].start=start;
				mstr[i].end=end;
				mstr[i].paired = 1;
				mstr[i].ntree = numberOfTrees;
				mstr[i].alignmentsdir = opt.print_alignments_dir;
				j=j+divideFile;
				mstr[i].str->taxonPath =(char**) malloc((end-start)*(sizeof(char *)));
				for(k=0; k<end-start; k++){
					mstr[i].str->taxonPath[k] = malloc((max_name_length+max_lineTaxonomy+120)*(sizeof(char)));
				}
			}
			for (i=0; i<opt.number_of_cores;i++){
				pthread_create(&threads[i], NULL, runAssignmentOnChunk_WithBWA, &mstr[i]);
			}
			for ( i=0; i<opt.number_of_cores;i++){
				pthread_join(threads[i], NULL);
			}
			for ( i=0; i<opt.number_of_cores; i++){
				for ( j=0; j<(mstr[i].end-mstr[i].start); j++){
					fprintf(results,"%s\n",mstr[i].str->taxonPath[j]);
				}
			}
			for ( i=0; i<opt.number_of_cores; i++){
				for(j=0; j<mstr[i].end-mstr[i].start; j++){
					free(mstr[i].str->taxonPath[j]);
				}
				free(mstr[i].str->taxonPath);
			}
		}
		fclose(results);
		gzclose(seqinfile_1);
		gzclose(seqinfile_2);
		if (opt.fastq==0){
			for(i=0; i<numberOfLinesToRead/2; i++){
				free(pairedQueryMat->query1Mat[i]);
				free(pairedQueryMat->query2Mat[i]);
				free(pairedQueryMat->forward_name[i]);
				free(pairedQueryMat->reverse_name[i]);
			}
		}else{
			for(i=0; i<numberOfLinesToRead/4; i++){
				free(pairedQueryMat->query1Mat[i]);
				free(pairedQueryMat->query2Mat[i]);
				free(pairedQueryMat->forward_name[i]);
				free(pairedQueryMat->reverse_name[i]);
			}
		}
		free(pairedQueryMat->query1Mat);
		free(pairedQueryMat->query2Mat);
		free(pairedQueryMat->forward_name);
		free(pairedQueryMat->reverse_name);
		free(pairedQueryMat);
	}
	for(i=0; i<opt.number_of_cores; i++){
		if (opt.fastq == 0){
			if (strcmp("single",opt.paired_or_single)==0){
				freeMemForResults(mstr[i].str, numberOfLinesToRead/2, opt.number_of_cores, numberOfTrees, 0, opt.use_nw, opt.use_leaf_portion,maxNumSpec,number_of_total_nodes);
			}else{
				freeMemForResults(mstr[i].str, numberOfLinesToRead/2, opt.number_of_cores, numberOfTrees, 1, opt.use_nw, opt.use_leaf_portion,maxNumSpec,number_of_total_nodes);
			}
		}else{
			if (strcmp("single",opt.paired_or_single)==0){
				freeMemForResults(mstr[i].str, numberOfLinesToRead/4, opt.number_of_cores, numberOfTrees, 0, opt.use_nw, opt.use_leaf_portion,maxNumSpec,number_of_total_nodes);
			}else{
				freeMemForResults(mstr[i].str, numberOfLinesToRead/4, opt.number_of_cores, numberOfTrees, 1, opt.use_nw, opt.use_leaf_portion,maxNumSpec,number_of_total_nodes);
			}
		}
	}
	for(i=0; i<numberOfTrees; i++){
		for(j=0; j<2*numspecArr[i]-1; j++){
			for(k=0; k<numbaseArr[i]; k++){
				free(treeArr[i][j].posteriornc[k]);
			}
			free(treeArr[i][j].posteriornc);
		}
		for(j=numspecArr[i]-1; j<(2*numspecArr[i]-1); j++){
			free(treeArr[i][j].name);
		}
		free(treeArr[i]);
	}
	free(treeArr);
	free(numbaseArr);
	free(rootArr);
	free(numspecArr);
}
