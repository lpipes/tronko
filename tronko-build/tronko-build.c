#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <sys/mman.h>
#include <fcntl.h>
#include "global.h"
#include "opt.h"
#include "options.h"
#include "hashmap.h"
HASHMAP(int, struct masterArr) mastermap;
struct node **treeArr;
int numspec, numbase, **seq, numundspec[MAXNUMBEROFINDINSPECIES+1];
char **nodeIDs;
char ***nodeIDsArr;
int root,tip,comma=0; /*globals used to read in the tree*/
double Logfactorial[MAXNUMBEROFINDINSPECIES];
double LRVEC[STATESPACE][STATESPACE], RRVEC[STATESPACE][STATESPACE], RRVAL[STATESPACE], PMAT1[STATESPACE][STATESPACE], PMAT2[STATESPACE][STATESPACE];
double LRVECnc[4][4], RRVECnc[4][4], RRVALnc[4], PMATnc[2][4][5];
double *statevector, UFC, *UFCnc, **templike_nc;
int COUNT2;
int COUNT;
double *localpi;
int localnode;
double parameters[10] = {0.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
type_of_PP ***PP;
type_of_PP ***PPcopy;
type_of_PP ****PP_Arr;
char ****taxonomyArr;
int *SPscoreArr;
int *numspecArr, *numbaseArr, ***seqArr, *rootArr;
double minVariance = 99999999999999;
int minVarNode = -1;
int returnNode=-1;
int *partitionSizes;
int *nodesToCut;
int *nodesToCutMinVar;
node **treeArr;


/*This function calculates the number of  desdencents of each node in the tree stored in tree[node].nd*/
int get_number_descendantsArr(int node, struct masterArr *m){
	if (m->tree[node].up[0]==-1) return (m->tree[node].nd=1);
	else return (m->tree[node].nd=(get_number_descendantsArr(m->tree[node].up[0],m)+get_number_descendantsArr(m->tree[node].up[1],m)));
}
void changePP_Arr (int node, int whichRoot){
	int child0 = treeArr[whichRoot][node].up[0];
	int child1 = treeArr[whichRoot][node].up[1];
	int i,j;
	if ( treeArr[whichRoot][node].up[0]==-1 && treeArr[whichRoot][node].up[1]==-1){
		for (i=0; i<numbaseArr[whichRoot]; i++){
			if ( seqArr[whichRoot][node-numspecArr[whichRoot]+1][i] == 4 ){
				for(j=0;j<4;j++){
					treeArr[whichRoot][node].posteriornc[i][j]=-1;
				}
			}
		}
		return;
	}else{
		changePP_Arr(child0,whichRoot);
		changePP_Arr(child1,whichRoot);
		return;
	}
}
void changePP_parents_Arr(int node, int whichRoot){
	int parent = treeArr[whichRoot][node].down;
	int child0 = treeArr[whichRoot][node].up[0];	
	int child1 = treeArr[whichRoot][node].up[1];
	if ( parent==-1 ){
		return;
	}
	if ( child0== -1 && child1== -1 ){
		changePP_parents_Arr(parent,whichRoot);
		return;
	}
	int i,j;
	for( i=0; i<numbaseArr[whichRoot]; i++){
		for (j=0; j<4; j++){
			if ( treeArr[whichRoot][child0].posteriornc[i][j] == -1 && treeArr[whichRoot][child1].posteriornc[i][j] == -1 ){
				treeArr[whichRoot][node].posteriornc[i][j]=-1;
			}
		}
	}
	changePP_parents_Arr(parent,whichRoot);
}
void assignDepthArr(int node0, int node1, int depth, struct masterArr *m){
	if( node0 != -1 && node1 != -1){
		m->tree[node0].depth = depth;
		m->tree[node1].depth = depth;
		assignDepthArr(m->tree[node0].up[0], m->tree[node0].up[1],depth+1,m);
		assignDepthArr(m->tree[node1].up[0], m->tree[node1].up[1],depth+1,m);
	}
}
void findMaxTaxName(FILE* file, int* specs){
	int max_tax_name=0;
	char buffer[BUFFER_SIZE];
	char *s;
	char *lineTaxonomy;
	int max_line = 0;
	while( fgets(buffer,BUFFER_SIZE,file) != NULL ){
		s = strtok(buffer,"\t");
		lineTaxonomy = strtok(NULL,"\n");
		int taxsize = strlen(lineTaxonomy);
		if ( max_line < taxsize ){
			max_line = taxsize;
		}
		int j=6;
		while(j>-1){
			s = strtok_r(lineTaxonomy,";",&lineTaxonomy);
			int size = strlen(s);
			if ( max_tax_name < size ){
				max_tax_name = size;
			}
			j--;
		}
	}
	if ( specs[0] < max_tax_name ){
		specs[0] = max_tax_name;
	}
	if (specs[1] < max_line ){
		specs[1] = max_line;
	}
}
void assignTaxonomyToLeavesArr(int node,char *tax,struct masterArr *m, int max_nodename, int max_tax_name){
	int child0 = m->tree[node].up[0];
	int child1 = m->tree[node].up[1];
	char buffer[BUFFER_SIZE];
	char *s;
	char *lineAccession, *lineTaxonomy, *taxLevelName;
	char fullTaxonomy[BUFFER_SIZE];
	int startIndex;
	char accessionID[max_nodename];
	int i;
	FILE *taxonomy_file;
	if ( child0 == -1 && child1 == -1 ){
		if (( taxonomy_file = fopen(tax,"r")) == (FILE *) NULL) printf("*** taxonomy file could not be opened.\n");
		while( fgets(buffer,BUFFER_SIZE,taxonomy_file) != NULL){
			lineAccession = strtok(buffer,"\t");
			lineTaxonomy = strtok(NULL,"\n");
			assert(strlen(lineAccession) <= max_nodename);
			strncpy(&accessionID[0], lineAccession, max_nodename);
			if ( strcmp(accessionID,m->tree[node].name)==0 ){
				m->tree[node].taxIndex[0]=node-m->numspec+1;
				m->tree[node].taxIndex[1]=0;
				int j=6;
				while(j>-1){
					taxLevelName = strtok_r(lineTaxonomy,";",&lineTaxonomy);
					strncpy(m->taxonomy[node-m->numspec+1][j], taxLevelName, max_tax_name); 
					j--;
				}
			}
		}
		fclose(taxonomy_file);
		return;
	}else{
		m->tree[node].taxIndex[0] = -1;
		m->tree[node].taxIndex[1] = -1;
	}
	if(child0 != -1){
		assignTaxonomyToLeavesArr(child0,tax,m,max_nodename,max_tax_name);
	}
	if(child1 != -1){
		assignTaxonomyToLeavesArr(child1,tax,m,max_nodename,max_tax_name);
	}
}
// by default calloc returns a pointer to allocated memory that has been
// initialized to zero (vs malloc, which just returns a pointer to allocated
// memory and doesn't initialize it), but it's good practice to check that the
// returned pointer is valid, since calloc (and malloc) will return NULL if
// the system is out of memory.
void *calloc_check(size_t nmemb, size_t size){
	void *ptr = calloc(nmemb, size);
	if(!ptr){   // if calloc returns NULL, we've run out of memory
		fprintf(stderr, "Out of memory!\n");
		exit(1);
	}
	return ptr;
}
void allocateMemoryForTaxArr(int whichPartitions, int max_tax_name){
	int phylogenyLevels = 7;
	int max_tax_name_len = max_tax_name;
	int i,j,k;
	// allocate the array to contain the different 'base' taxonomies / lists
	taxonomyArr = (char ****)calloc_check(whichPartitions,sizeof(char ***));
	for(k=0;k<whichPartitions;k++){
		taxonomyArr[k] = (char ***)calloc_check(numspecArr[k], sizeof(char **));
		for(i=0; i<numspecArr[k]; i++){
			// first, allocate the space for each taxonomy  (= list of taxonomic
			// classifications, from most specific to least specific)
			taxonomyArr[k][i] = (char **)calloc_check(phylogenyLevels, sizeof(char *));
			for(j=0; j<phylogenyLevels; j++){
				taxonomyArr[k][i][j] = (char *)calloc_check(max_tax_name_len, sizeof(char));
			}
		}
	}
}
int* getTaxonomyArr(int node, struct masterArr *m){
	int* taxIndexA = NULL;
	int* taxIndexB = NULL;
	int* taxonomyOfNode = malloc(sizeof(int)*2);
	taxonomyOfNode[0] = -1;
	taxonomyOfNode[1] = -1;
	if ( m->tree[node].up[0] != -1 ){
		if ( m->tree[node].taxIndex[0] == -1 ){
			taxIndexA = getTaxonomyArr(m->tree[node].up[0],m);
			taxIndexB = getTaxonomyArr(m->tree[node].up[1],m);
			if ( taxIndexA[0]==-1 || taxIndexB[0]==-1 ){
				m->tree[node].taxIndex[0]=-1;
				m->tree[node].taxIndex[1]=-1;
			}else{
				int i=0;
				int phylogenyLevel=0;
				int maxABLevel = (taxIndexA[1] > taxIndexB[1]) ? taxIndexA[1] : taxIndexB[1];
				//int maxABLevel = 0;
				for(i=maxABLevel;i<7;i++){
					if ( strcmp(m->taxonomy[taxIndexA[0]][i],m->taxonomy[taxIndexB[0]][i])==0){
						phylogenyLevel = i;
						taxonomyOfNode[0] = taxIndexA[0];
						taxonomyOfNode[1] = phylogenyLevel;
						m->tree[node].taxIndex[0] = taxIndexA[0];
						m->tree[node].taxIndex[1] = phylogenyLevel;
						break;
					}
				}
			}
		}
	}else{
		taxonomyOfNode[0] = m->tree[node].taxIndex[0];
		taxonomyOfNode[1] = m->tree[node].taxIndex[1];
	}
	if(taxIndexA != NULL) free(taxIndexA);
	if(taxIndexB != NULL) free(taxIndexB);
	return taxonomyOfNode;
}
int* getTaxonomyArr_UsePartitions(int node, int whichPartitions, char**** taxonomyArr_heap){
	int* taxIndexA = NULL;
	int* taxIndexB = NULL;
	int* taxonomyOfNode = malloc(sizeof(int)*2);
	taxonomyOfNode[0] = -1;
	taxonomyOfNode[1] = -1;
	if ( treeArr[whichPartitions][node].up[0] != -1 ){
		if ( treeArr[whichPartitions][node].taxIndex[0] == -1 ){
			taxIndexA = getTaxonomyArr_UsePartitions(treeArr[whichPartitions][node].up[0],whichPartitions,taxonomyArr_heap);
			taxIndexB = getTaxonomyArr_UsePartitions(treeArr[whichPartitions][node].up[1],whichPartitions,taxonomyArr_heap);
			if ( taxIndexA[0]==-1 || taxIndexB[0]==-1 ){
				treeArr[whichPartitions][node].taxIndex[0]=-1;
				treeArr[whichPartitions][node].taxIndex[1]=-1;
			}else{
				int i=0;
				int phylogenyLevel=0;
				int maxABLevel = (taxIndexA[1] > taxIndexB[1]) ? taxIndexA[1] : taxIndexB[1];
				//int maxABLevel = 0;
				for(i=maxABLevel;i<7;i++){
					if ( strcmp(taxonomyArr_heap[whichPartitions][taxIndexA[0]][i],taxonomyArr_heap[whichPartitions][taxIndexB[0]][i])==0){
						phylogenyLevel = i;
						taxonomyOfNode[0] = taxIndexA[0];
						taxonomyOfNode[1] = phylogenyLevel;
						treeArr[whichPartitions][node].taxIndex[0] = taxIndexA[0];
						treeArr[whichPartitions][node].taxIndex[1] = phylogenyLevel;
						break;
					}
				}
			}
		}
	}else{
		taxonomyOfNode[0] = treeArr[whichPartitions][node].taxIndex[0];
		taxonomyOfNode[1] = treeArr[whichPartitions][node].taxIndex[1];
	}
	if(taxIndexA != NULL) free(taxIndexA);
	if(taxIndexB != NULL) free(taxIndexB);
	return taxonomyOfNode;
}
void readSeqArr(FILE *partitionsFile, int maxname, struct masterArr *master){
	char buffer[FASTA_MAXLINE];
	int i, j, m, k=0, row=0;
	char c;
	int size;
	char nodename[maxname+1];
	char *s;
	int firstIter=1;
	int index=0;
	while( fgets(buffer,FASTA_MAXLINE,partitionsFile) != NULL){
		s = strtok(buffer,"\n");
		size = strlen(s);
		if ( buffer[0] == '>'){
			if ( size > maxname ){ size = maxname; }
			for(i=1; buffer[i]!='\0'; i++){
				nodename[i-1]=buffer[i];
			}
			nodename[size-1]='\0';
			strcpy(master->names[index], nodename);
			index++;
		}else{
			if ( firstIter==1 ){ 
				int l=0;
				int length=0;
				for(l=0; buffer[l] != '\0'; l++){
					length++;
				}
				master->numbase=length;
				master->msa = (int**)malloc(master->numspec*sizeof(int*));
				for (i=0; i<master->numspec; i++){
					master->msa[i]=(int *)malloc(master->numbase*(sizeof(int)));
				}
				firstIter=0;
				i=0;
			}
			for( m=0; m<master->numbase; m++){
				c=tolower(buffer[m]);
				if ((c!=' ')&&(c!='\t')) {
					if (c=='a' || c=='A') master->msa[row][m] = 0;
					else if (c=='c' || c=='C') master->msa[row][m] = 1;
					else if (c=='g' || c=='G') master->msa[row][m] = 2;
					else if (c=='t' || c=='T') master->msa[row][m] = 3;
					else if (c=='n'||c=='-'||c=='~' || c=='N') master->msa[row][m] = 4;
					else{
						printf("\nBAD BASE (%c) in species %i base %i",c ,row+1,m+1);
						scanf("%i",&row);
						exit(-1);
					}
				}
			}
			row++;
		}
	}
	if (master->numspec < 3) {printf("This is for more than two seq.s only!\n"); exit(-1);}
}
int *findLeavesOfMinVarArr(int node, int *leafNodeList, int size, struct masterArr *m){
	int child0 = m->tree[node].up[0];
	int child1 = m->tree[node].up[1];
	int i=0;
	int currentPos=-1;
	if ( node == minVarNode ){ return leafNodeList; }
	if (child0 ==-1 && child1 == -1){
		for(i=size-1;i>=0;i--){
			if (leafNodeList[i] == -1){
				currentPos = i;
			}
		}
		leafNodeList[currentPos]=node;
		return leafNodeList;
	}
	findLeavesOfMinVarArr(child0,leafNodeList,size,m);
	findLeavesOfMinVarArr(child1,leafNodeList,size,m);
	return leafNodeList;
}
char *removeDashArr(int *sequenceWithDash,int length, char *sequenceWithoutDash){
	int i,j;
	int countDash;
	countDash=0;
	for(i=0;i<length;i++){
		if( sequenceWithDash[i] < 4){
			if ( sequenceWithDash[i]==0 ){ sequenceWithoutDash[countDash]='A'; }
				else if (sequenceWithDash[i]==1 ){ sequenceWithoutDash[countDash]='C'; }
				else if (sequenceWithDash[i]==2 ){ sequenceWithoutDash[countDash]='G'; }
				else { sequenceWithoutDash[countDash]='T'; }
			countDash++;
		}
	}
	sequenceWithoutDash[countDash]='\0';
	return sequenceWithoutDash;
}
void printPartitionsToFileArr(int *partition1,int partition1size, int *partition2,int partition2size, int *partition3, int partition3size,int whichPartitions, int partitionCount, Options opt, struct masterArr *m){
	int i=0;
	int j=0;
	char buf_fasta[BUFFER_SIZE];
	char buf_tax[BUFFER_SIZE];
	int which = partitionCount;
	char *seqWithoutDash;
	snprintf(buf_fasta,BUFFER_SIZE,"%s/partition%d.fasta",opt.partitions_directory,which);
	FILE *p1=fopen(buf_fasta,"w");
	snprintf(buf_tax,BUFFER_SIZE,"%s/partition%d_taxonomy.txt",opt.partitions_directory,which);				
	FILE *p1_tax=fopen(buf_tax,"w");
	seqWithoutDash = (char *)malloc((m->numbase+1)*sizeof(char));
	for(i=0; i<m->numbase+1; i++){
		seqWithoutDash[i]='\0';
	}
	for(i=0;i<partition1size;i++){
		fprintf(p1,">%s\n",m->names[partition1[i]-m->numspec+1]);
		seqWithoutDash = removeDashArr(m->msa[partition1[i]-m->numspec+1],m->numbase,seqWithoutDash);
		fprintf(p1,"%s\n",seqWithoutDash);
		int taxnode=0;
		fprintf(p1_tax,"%s\t%s;%s;%s;%s;%s;%s;%s\n",m->names[partition1[i]-m->numspec+1],m->taxonomy[m->tree[partition1[i]].taxIndex[0]][6],m->taxonomy[m->tree[partition1[i]].taxIndex[0]][5],m->taxonomy[m->tree[partition1[i]].taxIndex[0]][4],m->taxonomy[m->tree[partition1[i]].taxIndex[0]][3],m->taxonomy[m->tree[partition1[i]].taxIndex[0]][2],m->taxonomy[m->tree[partition1[i]].taxIndex[0]][1],m->taxonomy[m->tree[partition1[i]].taxIndex[0]][0]);
		for (j=0; j<m->numbase; j++){
			seqWithoutDash[j]='\0';
		}
	}
	fclose(p1);
	fclose(p1_tax);
	which = partitionCount+1;
	snprintf(buf_fasta,BUFFER_SIZE,"%s/partition%d.fasta",opt.partitions_directory,which);				
	FILE *p2=fopen(buf_fasta,"w");
	snprintf(buf_tax,BUFFER_SIZE,"%s/partition%d_taxonomy.txt",opt.partitions_directory,which);				
	FILE *p2_tax=fopen(buf_tax,"w");
	for(i=0;i<partition2size;i++){
		fprintf(p2,">%s\n",m->names[partition2[i]-m->numspec+1]);
		seqWithoutDash = removeDashArr(m->msa[partition2[i]-m->numspec+1],m->numbase,seqWithoutDash);
		fprintf(p2,"%s\n",seqWithoutDash);
		fprintf(p2_tax,"%s\t%s;%s;%s;%s;%s;%s;%s\n",m->names[partition2[i]-m->numspec+1],m->taxonomy[m->tree[partition2[i]].taxIndex[0]][6],m->taxonomy[m->tree[partition2[i]].taxIndex[0]][5],m->taxonomy[m->tree[partition2[i]].taxIndex[0]][4],m->taxonomy[m->tree[partition2[i]].taxIndex[0]][3],m->taxonomy[m->tree[partition2[i]].taxIndex[0]][2],m->taxonomy[m->tree[partition2[i]].taxIndex[0]][1],m->taxonomy[m->tree[partition2[i]].taxIndex[0]][0]);
		for(j=0; j<m->numbase+1; j++){
			seqWithoutDash[j]='\0';
		}
	}
	fclose(p2);
	fclose(p2_tax);
	which=partitionCount+2; 
	snprintf(buf_fasta,BUFFER_SIZE,"%s/partition%d.fasta",opt.partitions_directory,which);				
	FILE *p3=fopen(buf_fasta,"w");
	snprintf(buf_tax,BUFFER_SIZE,"%s/partition%d_taxonomy.txt",opt.partitions_directory,which);				
	FILE *p3_tax=fopen(buf_tax,"w");
	for(i=0;i<partition3size;i++){
		fprintf(p3,">%s\n",m->names[partition3[i]-m->numspec+1]);
		seqWithoutDash = removeDashArr(m->msa[partition3[i]-m->numspec+1],m->numbase,seqWithoutDash);
		fprintf(p3,"%s\n",seqWithoutDash);
		fprintf(p3_tax,"%s\t%s;%s;%s;%s;%s;%s;%s\n",m->names[partition3[i]-m->numspec+1],m->taxonomy[m->tree[partition3[i]].taxIndex[0]][6],m->taxonomy[m->tree[partition3[i]].taxIndex[0]][5],m->taxonomy[m->tree[partition3[i]].taxIndex[0]][4],m->taxonomy[m->tree[partition3[i]].taxIndex[0]][3],m->taxonomy[m->tree[partition3[i]].taxIndex[0]][2],m->taxonomy[m->tree[partition3[i]].taxIndex[0]][1],m->taxonomy[m->tree[partition3[i]].taxIndex[0]][0]);
		for(j=0; j<m->numbase; j++){
			seqWithoutDash[j]='\0';
		}
	}
	free(seqWithoutDash);
	fclose(p3);
	fclose(p3_tax);
}
double calculateSPArr(struct masterArr *m){
	int i,j,k;
	int numpairs=0;
	double SPscore=0;
	j=0;
	int *partition = (int*)malloc(sizeof(int)*m->numspec);
	for(i=m->numspec-1; i<2*m->numspec-1; i++){
		partition[j]=i;
		j++;
	}	
	for (i = 0; i<m->numbase; i++){
		for( j=0; j<m->numspec; j++){
			for( k=j+1; k<m->numspec; k++){
				if (m->msa[partition[j]-m->numspec+1][i]==m->msa[partition[k]-m->numspec+1][i]){
					SPscore=SPscore+3;
				}else if (m->msa[partition[j]-m->numspec+1][i] != 4 && m->msa[partition[k]-m->numspec+1][i] != 4 && m->msa[partition[j]-m->numspec+1][i]!=m->msa[partition[k]-m->numspec+1][i] ){
					SPscore=SPscore-2;
				}else{
					SPscore--;
				}
				numpairs++;
			}
		}
	}
	free(partition);
	//printf("raw SPscore: %lf\n",SPscore);
	//printf("numpairs: %d\n",numpairs);
	double exponent;
	exponent = (-5.0/2.0);
	//SPscore=SPscore/(m->numspec*exp(exponent));
	//SPscore=SPscore/(m->numspec*exp(-5/2));
	SPscore=SPscore/m->numspec;
	SPscore=SPscore/numpairs;
	printf("SPscore: %lf\n",SPscore);
	return SPscore;
}
void findMinVarianceArr(int node, int size, struct masterArr *m){
	int child0 = m->tree[node].up[0];
	int child1 = m->tree[node].up[1];
	int parent = m->tree[node].down;
	int i=0;
	if (node==-1){ return; }
	if ( parent == -1 ){
		findMinVarianceArr(child0,size,m);
		findMinVarianceArr(child1,size,m);
		return;
	}
	if (child0 == -1 ){ return; }
	if (child1 == -1 ){ return; }	
	int num_children0 = m->tree[child0].nd;
	int num_children1 = m->tree[child1].nd;
	int num_ancestors = size-num_children1-num_children0;
	double mean = (double)(num_children0 + num_children1 + num_ancestors )/3;
	double variance = (double)((num_ancestors-mean)*(num_ancestors-mean) + (num_children0-mean)*(num_children0-mean) + (num_children1-mean)*(num_children1-mean))/3;
	if (minVariance > variance){
		minVariance = variance;
		minVarNode = node;
	}
	findMinVarianceArr(child0,size,m);
	findMinVarianceArr(child1,size,m);
	return;
}
void createNewRoots(int rootCount, Options opt, int max_nodename, int max_lineTaxonomy, struct masterArr *m){
	int i,j,k,count;
	char buffer[BUFFER_SIZE];
	int *leaves;
	FILE *partition, *partitionTree;
	minVariance = 99999999999999;
	findMinVarianceArr(m->root,m->numspec,m);
	int *partition1 = (int *)malloc(sizeof(int)*m->tree[m->tree[minVarNode].up[0]].nd);
	for(i=0; i<m->tree[m->tree[minVarNode].up[0]].nd; i++){
		partition1[i]=-1;
	}
	int *partition2 = (int *)malloc(sizeof(int)*m->tree[m->tree[minVarNode].up[1]].nd);
	for(i=0; i<m->tree[m->tree[minVarNode].up[1]].nd; i++){
		partition2[i]=-1;
	}
	int *partition3 = (int *)malloc(sizeof(int)*(m->tree[m->root].nd-m->tree[minVarNode].nd));
	for(i=0; i<(m->tree[m->root].nd-m->tree[minVarNode].nd); i++){
		partition3[i]=-1;
	}
	partition1 = findLeavesOfMinVarArr(m->tree[minVarNode].up[0],partition1,m->tree[m->tree[minVarNode].up[0]].nd,m);
	partition2 = findLeavesOfMinVarArr(m->tree[minVarNode].up[1],partition2,m->tree[m->tree[minVarNode].up[1]].nd,m);
	partition3 = findLeavesOfMinVarArr(m->root,partition3,(m->tree[m->root].nd-m->tree[minVarNode].nd),m);
	int *partitionSizes = (int *)malloc(sizeof(int)*3);
	partitionSizes[0]=m->tree[m->tree[minVarNode].up[0]].nd;
	partitionSizes[1]=m->tree[m->tree[minVarNode].up[1]].nd;
	partitionSizes[2]=(m->tree[m->root].nd-m->tree[minVarNode].nd);
	if ( partitionSizes[0] < 4 || partitionSizes[1] < 4 || partitionSizes[2] < 4){
		SPscoreArr[rootCount]=1;
		return;
	}
	int partitionCount=-1;
	for(i=MAX_NUMBEROFROOTS-1;i>=0;i--){
		if(SPscoreArr[i]==-1){
			partitionCount=i;
		}
	}
	SPscoreArr[partitionCount]=1;
	SPscoreArr[partitionCount+1]=1;
	SPscoreArr[partitionCount+2]=1;
	partitionCount++;
	if (partitionCount > opt.restart){
		printPartitionsToFileArr(partition1,partitionSizes[0],partition2,partitionSizes[1],partition3,partitionSizes[2],rootCount,partitionCount,opt,m);
	}
	for(i=0; i<m->numspec; i++){
		for(j=0; j<7; j++){
			free(m->taxonomy[i][j]);
		}
		free(m->taxonomy[i]);
		free(m->msa[i]);
		free(m->names[i]);
	}
	for(i=0; i<2*m->numspec-1; i++){
		free(m->tree[i].like);
		free(m->tree[i].posterior);
		free(m->tree[i].name);
	}
	free(m->taxonomy);
	free(m->msa);
	free(m->names);
	free(m->tree);
	hashmap_remove(&mastermap,m->index);
	free(m);
	rootCount=partitionCount;
	printf("rootCount=%d\n",rootCount);
	int which,status;
	char buf[BUFFER_SIZE];
	for(which=rootCount;which<rootCount+3;which++){
		if (rootCount > opt.restart){
			pid_t pid;
			int ret = 1;
			char buf3[BUFFER_SIZE];	
			char buf2[BUFFER_SIZE];
			snprintf(buf2,BUFFER_SIZE,"%s/partition%d.fasta",opt.partitions_directory,which);
			printf("buf2 is %s\n",buf2);
			snprintf(buf3,BUFFER_SIZE,"%s/partition%d_MSA.fasta",opt.partitions_directory,which);
			printf("buf3 is %s\n",buf3);
			pid=fork();
			if (pid==-1){
			//pid==-1 means error occured
				printf("can't fork, error occured\n");
			}else if (pid==0){
				printf("child process, pid = %u\n",getpid());
				char *arguments[] = {"famsa",buf2,buf3,NULL};
				printf("ARGUMENTS: %s %s %s\n",arguments[0],arguments[1],arguments[2]);
				execvp("famsa",arguments);
				exit(0);
			}else{
				printf("parent process, pid = %u\n",getppid());
				if (waitpid(pid, &status, 0) > 0){
					if (WIFEXITED(status) && !WEXITSTATUS(status))
					printf("program execution successful\n");
				else if (WIFEXITED(status) && WEXITSTATUS(status)) {
					if (WEXITSTATUS(status) == 127) {
						printf("execv failed\n");
					}else printf("program terminated normally but returned a non-zero status\n");
					}else printf("program didn't terminate normally\n");
					}else{
						printf("waitpid() failed\n");
					}
			}
			snprintf(buf,BUFFER_SIZE,"sed -i ':a; $!N; /^>/!s/\\n\\([^>]\\)/\\1/; ta; P; D' %s/partition%d_MSA.fasta",opt.partitions_directory,which);
			status = system(buf);
			snprintf(buf,BUFFER_SIZE,"fasta2phyml.pl %s/partition%d_MSA.fasta",opt.partitions_directory,which);
			status = system(buf);
			snprintf(buf,BUFFER_SIZE,"raxmlHPC-PTHREADS --silent -m GTRGAMMA -w %s/ -n partition%d -p 1234 -T 8 -s %s/partition%d_MSA.phymlAln",opt.partitions_directory,which,opt.partitions_directory,which);
			status = system(buf);
			char buf4[BUFFER_SIZE];
			char buf5[BUFFER_SIZE];
			snprintf(buf4,BUFFER_SIZE,"%s/RAxML_bestTree.partition%d",opt.partitions_directory,which);
			snprintf(buf5,BUFFER_SIZE,"%s/RAxML_bestTree.partition%d.reroot",opt.partitions_directory,which);
			pid=fork();
			if (pid==-1){
				printf("can't fork, error occured\n");
			}else if (pid==0){	
				printf("child process reroot, pid = %u\n",getpid());
				int reroot_file;
				reroot_file = open(buf5, O_WRONLY | O_TRUNC | O_CREAT, S_IRUSR | S_IRGRP | S_IWGRP | S_IWUSR);
				dup2(reroot_file,1);
				close(reroot_file);
				char *arguments[3];
				arguments[0]="/space/s1/lenore/software/newick-utils-1.6/src/nw_reroot";
				arguments[1]=buf4;
				arguments[2]=NULL;
				ret=execv("/space/s1/lenore/software/newick-utils-1.6/src/nw_reroot",arguments);
				exit(0);
			}else{
				printf("parent process, pid = %u\n",getppid());
				if (waitpid(pid, &status, 0) > 0){
					if (WIFEXITED(status) && !WEXITSTATUS(status))
						printf("program execution successful\n");
					else if (WIFEXITED(status) && WEXITSTATUS(status)) {
						if (WEXITSTATUS(status) == 127) {
							printf("execv failed\n");
						}else printf("program terminated normally but returned a non-zero status\n");
					}else printf("program didn't terminate normally\n");
				}else{
					printf("waitpid() failed\n");
				}
			}
		}
		snprintf(buf,BUFFER_SIZE,"%s/partition%d_MSA.fasta",opt.partitions_directory,which);
		printf("buffer: %s\n",buf);
		struct masterArr *t=malloc(sizeof(masterArr));;
		//itoa(which-1,t->index,10);
		sprintf(t->index,"%d",which-1);
		if (NULL==(partition=fopen(buf,"r"))){ puts("Cannot open partition file 1!"); exit(-1);}
		t->numspec=setNumspecArr(partition);
		fclose(partition);
		t->tree=(struct node*)malloc((2*t->numspec-1)*sizeof(struct node));
		t->names=(char**)malloc(sizeof(char*)*t->numspec);
		for(i=0; i<t->numspec; i++){
			t->names[i]=(char*)malloc(sizeof(char)*(max_nodename+1));
		}
		t->msa=(int**)malloc(t->numspec*sizeof(int*));
		t->taxonomy=(char***)malloc(t->numspec*sizeof(char**));
		if (NULL==(partition=fopen(buf,"r"))){ puts("Cannot open partition file 2!"); exit(-1);}
		readSeqArr(partition,max_nodename,t);
		fclose(partition);
		allocateTreeArrMemory(t,max_nodename);
		snprintf(buf,BUFFER_SIZE,"%s/RAxML_bestTree.partition%d.reroot",opt.partitions_directory,which);
		if (NULL==(partitionTree=fopen(buf,"r"))){ puts("Cannot open partition tree file!"); exit(-1);}
		comma=0;
		tip=0;
		t->root=getcladeArr(partitionTree,t,max_nodename)-1;
		fclose(partitionTree);
		t->tree[t->root].down = -1;
		get_number_descendantsArr(t->root,t);
		int child0 = t->tree[t->root].up[0];
		int child1 = t->tree[t->root].up[1];
		t->tree[t->root].depth = 0;
		assignDepthArr(child0,child1,1,t);
		printf("NumbaseArr[%d]=%d, numspecArr[%d]=%d, rootArr[%d]=%d\n",t->index,t->numbase,t->index,t->numspec,t->index,t->root);
		snprintf(buf,BUFFER_SIZE,"%s/partition%d_taxonomy.txt",opt.partitions_directory,which);
		printf("buffer: %s\n",buf);
		for(j=0; j<t->numspec; j++){
			t->taxonomy[j] = (char **)calloc_check(7, sizeof(char *));
			for(k=0; k<7; k++){
				t->taxonomy[j][k] = (char *)calloc_check(max_lineTaxonomy, sizeof(char));
			}
		}
		assignTaxonomyToLeavesArr(t->root,buf,t,max_nodename,max_lineTaxonomy);
		getTaxonomyArr(t->root,t);
		hashmap_put(&mastermap, t->index, t);
		double initialSPscore=-1;
		if ( opt.use_spscore==1 && opt.use_min_leaves==0){
				initialSPscore = calculateSPArr(t);
				if (initialSPscore < /*0.05*/opt.sp_score){
					SPscoreArr[which-1]=0;
					createNewRoots(which-1,opt,max_nodename,max_lineTaxonomy,t);
				}else{
					SPscoreArr[which-1]=1;
				}
		}
		if ( opt.use_spscore==0 && opt.use_min_leaves==1  ){
			if ( t->numspec > opt.min_leaves ){
				SPscoreArr[which-1]=0;
				createNewRoots(which-1,opt,max_nodename,max_lineTaxonomy,t);
			}else{
				SPscoreArr[which-1]=1;
			}
		}
		if ( opt.use_spscore==1 && opt.use_min_leaves==1 ){
			initialSPscore = calculateSPArr(t);
			if (initialSPscore < /*0.05*/opt.sp_score && t->numspec > opt.min_leaves){
				createNewRoots(which-1,opt,max_nodename,max_lineTaxonomy,t);
			}else{
				SPscoreArr[which-1]=1;
			}
		}
	}
	free(partition1);
	free(partition2);
	free(partition3);
	free(partitionSizes);
}
void populate(int* nodes, char** taxnames, int node){
	int child0 = treeArr[0][node].up[0];
	int child1 = treeArr[0][node].up[1];
	int parent = treeArr[0][node].down;
	int i,j;
	int no_add=0;
	//if ( treeArr[0][node].taxIndex[1] == 2 ){ 
		for(i=0; i<100; i++){
			//if (strcmp(taxonomyArr[0][treeArr[0][node].taxIndex[0]][treeArr[0][node].taxIndex[1]],taxnames[i])==0){
			if (strcmp(taxonomyArr[0][treeArr[0][node].taxIndex[0]][2],taxnames[i])==0){
				if ( treeArr[0][node].depth > treeArr[0][nodes[i]].depth ){
					nodes[i]=node;
				}
				no_add=1;
			}
		}
		if ( no_add == 0 ){
			for(i=0; i<100; i++){
				if ( nodes[i]==-1 ){
					break;
				}
			}
			nodes[i] = node;
			//strcpy(taxnames[i],taxonomyArr[0][treeArr[0][node].taxIndex[0]][treeArr[0][node].taxIndex[1]]);
			strcpy(taxnames[i],taxonomyArr[0][treeArr[0][node].taxIndex[0]][2]);
		}
	//}else{
		if ( child0 != -1 && child1 != -1 ){
			populate(nodes,taxnames,child0);
			populate(nodes,taxnames,child1);
		}
	//}

}
/*void printFamilySeqs(){
	int* nodes = (int*)malloc(sizeof(int)*100);
	int i, j;
	for(i=0; i<100; i++){
		nodes[i]=-1;
	}
	char** taxnames = (char**)malloc(sizeof(char*)*100);
	for(i=0; i<100; i++){
		taxnames[i]=malloc(sizeof(char)*100);
		for(j=0; j<100; j++){
			taxnames[i][j]='\0';
		}
	}
	populate(nodes,taxnames,rootArr[0]);
	for(i=0; i<100; i++){
		if ( nodes[i] == -1 ){ break; }
		char* seq = (char*)malloc(sizeof(char)*1000);
		for(j=0; j<1000; j++){
			seq[j] = '\0';
		}
		printf(">%s\n",taxnames[i]);
		getSequenceInNode(0,nodes[i],seq);
		printf("%s\n",seq);
	}
}*/
int main(int argc, char **argv){
	Options opt;
	opt.min_leaves=0;
	opt.use_partitions=0;
	opt.use_spscore=0;
	opt.use_min_leaves=0;
	opt.sp_score = 0.05;
	opt.missing_data=1;
	opt.restart = 0;
	opt.number_of_partitions = 0;
	parse_options(argc, argv, &opt);
	int i, j, k, numberOfTrees;
	int max_nodename = 0;
	int max_tax_name = 0;
	int max_lineTaxonomy = 0;
	/*gzFile ref = Z_NULL;
	if (( ref = gzopen("../cluster_quality/PITS_new/clusters2_wo18_tronkodb/reference_tree.txt","r")) == (FILE *) NULL) fprintf(stderr,"MSA file could not be opened.\n");
	i=readReferenceTree(ref);
	gzclose(ref);
	printf("pass\n");
	exit(1);*/
	hashmap_init(&mastermap,hashmap_hash_string,strcmp);
	if (opt.number_of_trees==1 && opt.use_partitions==0){
		printf("Using a single tree... \n");
		numberOfTrees=1;
		struct masterArr *m = malloc(sizeof(masterArr));
		m->tree = malloc(sizeof(node *));
		//numspecArr = (int *)malloc(sizeof(int));
		//numbaseArr = (int *)malloc(sizeof(int));
		//rootArr = (int *)malloc(sizeof(int));
		int *specifications = (int*)malloc(3*sizeof(int));
		specifications[0]=0;
		specifications[1]=0;
		specifications[2]=0;
		FILE* infile;
		if (( infile = fopen(opt.msa_file,"r")) == (FILE *) NULL) fprintf(stderr,"MSA file could not be opened.\n");
		setNumspec(infile,specifications);
		fclose(infile);
		m->numspec = specifications[0];
		max_nodename = specifications[1];
		m->numbase = specifications[2];
		free(specifications);
		initlogfactorial();
		//nodeIDsArr = (char ***)malloc(sizeof(char**));
		//itoa(0,m->index,10);
		sprintf(m->index,"%d",0);
		m->tree = (struct node*)malloc((2*m->numspec-1)*sizeof(struct node));
		m->names = (char **)malloc(sizeof(char *)*m->numspec);
		for(i=0; i<m->numspec;i++) {
			m->names[i] = malloc((max_nodename+1)*sizeof(char));
		}
		//seqArr = (int ***)malloc(sizeof(int **));
		m->msa = (int **)malloc(m->numspec*sizeof(int *));
		for(i=0; i<m->numspec; i++){
			m->msa[i]=(int *)malloc(m->numbase*sizeof(int));
		}
		if (( infile = fopen(opt.msa_file,"r")) == (FILE *) NULL) fprintf(stderr,"MSA file could not be opened.\n");
		readseq(infile,max_nodename,m);
		fclose(infile);
		allocateTreeArrMemory(m,max_nodename);
		comma=0;
		tip=0;
		FILE* treefile;
		if (( treefile = fopen(opt.tree_file,"r")) == (FILE *) NULL) fprintf(stderr,"*** tree file could not be opened.\n");
		m->root=getcladeArr(treefile,m,max_nodename)-1;
		close(treefile);
		m->tree[m->root].down = -1;
		get_number_descendantsArr(m->root,m);
		int child0 = m->tree[m->root].up[0];
		int child1 = m->tree[m->root].up[1];
		m->tree[m->root].depth=0;
		assignDepthArr(child0,child1,1,m);
		int* tax_specs = (int*)malloc(2*sizeof(int));
		tax_specs[0]=0;
		tax_specs[1]=0;
		FILE* taxfile;
		if (( taxfile = fopen(opt.taxonomy_file,"r")) == (FILE *) NULL) fprintf(stderr,"Taxonomy file could not be opened.\n");
		findMaxTaxName(taxfile,tax_specs);
		fclose(taxfile);
		max_tax_name = tax_specs[0];
		max_lineTaxonomy = tax_specs[1];
		free(tax_specs);
		m->taxonomy = (char ***)calloc_check(m->numspec, sizeof(char **));
		for(j=0; j<m->numspec; j++){
			m->taxonomy[j] = (char **)calloc_check(7, sizeof(char *));
			for(k=0; k<7; k++){
				m->taxonomy[j][k] = (char *)calloc_check(max_lineTaxonomy, sizeof(char));
			}
		}
		assignTaxonomyToLeavesArr(m->root,opt.taxonomy_file,m,max_nodename,max_tax_name);
		getTaxonomyArr(m->root,m);
		hashmap_put(&mastermap,m->index,m);
		/*treeArr = malloc(sizeof(node*));
		treeArr[0] = m->tree;
		taxonomyArr = (char****)malloc(sizeof(char***));
		taxonomyArr[0] = m->taxonomy;
		seqArr = (int***)malloc(sizeof(int**));
		seqArr[0] = m->msa;
		numbaseArr = (int*)malloc(sizeof(int));
		numbaseArr[0] = m->numbase;
		numspecArr = (int*)malloc(sizeof(int));
		numspecArr[0] = m->numspec;
		rootArr = (int*)malloc(sizeof(int));
		rootArr[0] = m->root;*/
	}else{
		printf("Using %d trees... \n",opt.number_of_partitions);
		//treeArr = malloc(sizeof(node *)*opt.number_of_partitions);
		SPscoreArr = malloc(sizeof(int)*MAX_NUMBEROFROOTS);
		for(i=0;i<MAX_NUMBEROFROOTS;i++){
			SPscoreArr[i]=-1;
		}
		//numspecArr = (int *)malloc(sizeof(int)*opt.number_of_partitions);
		//numbaseArr = (int *)malloc(sizeof(int)*opt.number_of_partitions);
		//rootArr = (int *)malloc(sizeof(int)*opt.number_of_partitions);
		partition_files* pf = malloc(sizeof(partition_files));
		pf->tree_files = (char **)malloc(sizeof(char *)*opt.number_of_partitions);
		pf->msa_files = (char **)malloc(sizeof(char *)*opt.number_of_partitions);
		pf->tax_files = (char **)malloc(sizeof(char *)*opt.number_of_partitions);
		for (i=0; i<opt.number_of_partitions; i++){
			pf->tree_files[i] = (char *)malloc(sizeof(char)*MAXFILENAME);
			pf->msa_files[i] = (char *)malloc(sizeof(char)*MAXFILENAME);
			pf->tax_files[i] = (char *)malloc(sizeof(char)*MAXFILENAME);
		}
		readFilesInDir(opt.readdir,opt.number_of_partitions,pf);
		printf("Input files...\n");
		for (i=0; i<opt.number_of_partitions; i++){
			printf("%d : %s\n",i,pf->msa_files[i]);
			printf("%d : %s\n",i,pf->tree_files[i]);
			printf("%d : %s\n",i,pf->tax_files[i]);
		}
		int* tax_specs = (int*)malloc(2*sizeof(int));
		tax_specs[0]=0;
		tax_specs[1]=0;
		char buffer[BUFFER_SIZE];
		int* specifications = (int*)malloc(3*sizeof(int));
		specifications[0]=0;
		specifications[1]=0;
		specifications[2]=0;
		for(i=0; i<opt.number_of_partitions; i++){
			FILE *taxfiles;
			snprintf(buffer,BUFFER_SIZE,"%s/%s",opt.readdir,pf->tax_files[i]);
			if (NULL==(taxfiles=fopen(buffer,"r"))){ puts ("Cannot open tax file!"); exit(-1); }
			findMaxTaxName(taxfiles,tax_specs);
			fclose(taxfiles);
			FILE *msafiles;
			snprintf(buffer,BUFFER_SIZE,"%s/%s",opt.readdir,pf->msa_files[i]);
			if (NULL==(msafiles=fopen(buffer,"r"))){ puts ("Cannot open tax file!"); exit(-1); }
			setNumspec(msafiles,specifications);
			fclose(msafiles);
		}
		max_tax_name = tax_specs[0];
		max_lineTaxonomy = tax_specs[1];
		free(tax_specs);
		int max_numspec = specifications[0];
		max_nodename = specifications[1];
		int max_numbase = specifications[2];
		free(specifications);
		FILE *partition, *partitionTree;
		int status;
		int partition_count = opt.number_of_partitions;
		for(i=0; i<partition_count; i++){
			SPscoreArr[i]=0;
		}
		//taxonomyArr = (char ****)malloc(opt.number_of_partitions*sizeof(char ***));
		//nodeIDsArr = (char ***)malloc(opt.number_of_partitions*sizeof(char**));
		//seqArr = (int ***)malloc(opt.number_of_partitions*sizeof(int **));
		/*for(i=0; i<MAX_NUMBEROFROOTS; i++){
			nodeIDsArr[i] = (char**)malloc(sizeof(char*)*max_numspec);
			seqArr[i] = (int **)malloc(max_numspec*sizeof(int*));
			for(j=0; j<max_numspec; j++){
				nodeIDsArr[i][j]=malloc(max_nodename*sizeof(char));
				seqArr[i][j]=malloc(max_numbase*sizeof(int));
				for(k=0; k<max_numbase; k++){
					nodeIDsArr[i][j][k]='\0';
					seqArr[i][j][k]=0;
				}
			}
		}*/
		for(i=0; i<partition_count; i++){
			struct masterArr *m = malloc(sizeof(masterArr));
			//itoa(i,m->index,10);
			sprintf(m->index,"%d",i);
			snprintf(buffer,BUFFER_SIZE,"%s/%s",opt.readdir,pf->msa_files[i]);
			if (NULL==(partition=fopen(buffer,"r"))){ puts("Cannot open partition file!"); exit(-1);}
			m->numspec = setNumspecArr(partition);
			printf("m->numspec: %d\n",m->numspec);
			fclose(partition);
			m->tree=(struct node*)malloc((2*m->numspec-1)*sizeof(struct node));
			m->msa=(int**)malloc(m->numspec*sizeof(int*));
			m->taxonomy=(char***)malloc(m->numspec*sizeof(char**));
			m->names=(char**)malloc(m->numspec*sizeof(char*));
			for(j=0;j<m->numspec;j++){
				m->names[j]=(char*)malloc(sizeof(char)*(max_nodename+1));
			}
			if (NULL==(partition=fopen(buffer,"r"))){ puts("Cannot open partition file!"); exit(-1);}
			readSeqArr(partition,max_nodename,m);
			fclose(partition);
			allocateTreeArrMemory(m,max_nodename);
			snprintf(buffer,BUFFER_SIZE,"%s/%s",opt.readdir,pf->tree_files[i]);
			if (( partitionTree = fopen(buffer,"r")) == (FILE *) NULL) printf("*** tree file could not be opened.\n");
			comma=0;
			tip=0;
			m->root=getcladeArr(partitionTree,m,max_nodename)-1;
			fclose(partitionTree);
			m->tree[m->root].down = -1;
			get_number_descendantsArr(m->root,m);
			int child0 = m->tree[m->root].up[0];
			int child1 = m->tree[m->root].up[1];
			m->tree[m->root].depth=0;
			assignDepthArr(child0,child1,1,m);
			snprintf(buffer,BUFFER_SIZE,"%s/%s",opt.readdir,pf->tax_files[i]);
			m->taxonomy = (char ***)calloc_check(m->numspec, sizeof(char **));
			for(j=0; j<m->numspec; j++){
				m->taxonomy[j] = (char **)calloc_check(7, sizeof(char *));
				for(k=0; k<7; k++){
					m->taxonomy[j][k] = (char *)calloc_check(max_lineTaxonomy, sizeof(char));
				}
			}
			assignTaxonomyToLeavesArr(m->root,buffer,m,max_nodename,max_lineTaxonomy);
			getTaxonomyArr(m->root,m);
			hashmap_put(&mastermap, m->index, m);
			if ( opt.use_partitions==1 && m->numspec > opt.min_leaves ){
				SPscoreArr[i]=0;
				createNewRoots(i,opt,max_nodename,max_lineTaxonomy,m);
			}
		}
	}
	free(SPscoreArr);
	numberOfTrees=0;
	int key;
	struct masterArr* final;
	hashmap_foreach(key,final,&mastermap){
		numberOfTrees++;
	}
	printf("Number of trees: %d\n",numberOfTrees);
	treeArr = malloc(numberOfTrees*sizeof(node*));
	taxonomyArr = (char****)malloc(numberOfTrees*sizeof(char***));
	seqArr = (int***)malloc(numberOfTrees*sizeof(int**));
	numbaseArr = (int*)malloc(numberOfTrees*sizeof(int));
	numspecArr = (int*)malloc(numberOfTrees*sizeof(int));
	rootArr = (int*)malloc(numberOfTrees*sizeof(int));
	int index=0;
	hashmap_foreach(key,final,&mastermap){
		treeArr[index] = final->tree;
		taxonomyArr[index] =  final->taxonomy;
		seqArr[index] = final->msa;
		numbaseArr[index] = final->numbase;
		numspecArr[index] = final->numspec;
		rootArr[index] = final->root;
		//printtreeArr(index);
		index++;
	}
	allocatetreememory_for_nucleotide_Arr(numberOfTrees);
	for(i=0; i<numberOfTrees; i++){
		estimatenucparameters_Arr(parameters,i);
		getposterior_nc_Arr(parameters,i);
		//set_posteriors(i);
	}
	//FILE *for_monica = fopen("/space/s1/lenore/trout_copy/s2_copy/for_monica/OU061397_1_PP.txt","w");
	//if ( for_monica == NULL ){ printf("Error opening file!\n"); exit(1); }
	//for(i=0; i<numbaseArr[0]; i++){
	//	fprintf(for_monica,"%d",i);
	//	for(j=0; j<4; j++){
	//		fprintf(for_monica,"\t%.17g",treeArr[0][200].posteriornc[i][j]);
	//	}
	//	fprintf(for_monica,"\n");
	//}
	//fclose(for_monica);
	//exit(1);
	if (opt.missing_data==1){
		for(i=0; i<numberOfTrees; i++){
			changePP_Arr(rootArr[i],i);
			for(j=numspecArr[i]-1;j<2*numspecArr[i]-1;j++){
				changePP_parents_Arr(j,i);
			}
		}
	}
	printTreeFile(numberOfTrees,max_nodename,max_tax_name,max_lineTaxonomy,opt);
	/*hashmap_foreach(key,final,&mastermap){
		for(i=0; i<final->numspec; i++){
			free(final->names[i]);
			for(j=0; j<7; j++){
				free(final->taxonomy[i][j]);
			}
			free(final->taxonomy[i]);
			free(final->msa[i]);
		}
		free(final->taxonomy);
		free(final->names);
		free(final->msa);
	}*/
	for(i=0; i<numberOfTrees; i++){
		freeTreeMemory(i);
	}
	free(treeArr);
	//hashmap_cleanup(&mastermap);
	free(taxonomyArr);
	free(seqArr);
	free(numbaseArr);
	free(numspecArr);
	free(rootArr);
}
