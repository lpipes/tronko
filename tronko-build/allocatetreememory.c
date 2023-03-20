#include "allocatetreememory.h"

/*void allocatetreememmory_for_nucleotide()
{
  int i, j;

  for (i=0; i<(numspec*2-1); i++){
    tree[i].likenc = malloc(numbase*(sizeof(double *)));
    tree[i].posteriornc = malloc(numbase*(sizeof(double *)));
    for (j=0; j<numbase; j++){
      tree[i].likenc[j] = malloc(4*(sizeof(double)));
      tree[i].posteriornc[j] = malloc(4*(sizeof(double)));
    }
  }
  for (i=0;i<4;i++){
    PMATnc[MAX_NUMBEROFROOTS][0][i][4]=1.0;//missing data
    PMATnc[MAX_NUMBEROFROOTS][1][i][4]=1.0;//missing data
  }
}*/
void allocatetreememory_for_nucleotide_Arr(int numberOfTrees){
  int i, j, k;
  for(k=0; k<numberOfTrees; k++){
  	for (i=0; i<(numspecArr[k]*2-1); i++){
    	treeArr[k][i].likenc = malloc(numbaseArr[k]*(sizeof(double *)));
    	treeArr[k][i].posteriornc = malloc(numbaseArr[k]*(sizeof(double *)));
    	for (j=0; j<numbaseArr[k]; j++){
      		treeArr[k][i].likenc[j] = malloc(4*(sizeof(double)));
      		treeArr[k][i].posteriornc[j] = malloc(4*(sizeof(double)));
    	}
  	}
  }
  //for(k=0;k<numberOfTrees; k++){
  	for (i=0;i<4;i++){
    	PMATnc[0][i][4]=1.0;//missing data
    	PMATnc[1][i][4]=1.0;//missing data
  	}
  //}
}
void allocateTreeArrMemory(struct masterArr *m, int max_nodename){
	int i;
	for (i=0; i<m->numspec*2-1; i++){
		m->tree[i].like = malloc(STATESPACE*(sizeof(double)));
		m->tree[i].posterior = malloc(STATESPACE*(sizeof(double)));
		m->tree[i].name = malloc((max_nodename+1)*sizeof(char));
		m->tree[i].up[0] = -2;  // -2 is uninitialized / NULL value
		m->tree[i].up[1] = -2;
		m->tree[i].down = -2;
		m->tree[i].nd = -2;
		m->tree[i].depth = -2;
		m->tree[i].bl = -2.0;
		m->tree[i].likenc = NULL;
		m->tree[i].posteriornc = NULL;
		m->tree[i].s = -2;
		m->tree[i].numsites = -2;
		m->tree[i].spec = -2;
		m->tree[i].mrca = -2;
		memset(m->tree[i].name, '\0', max_nodename+1);  // name is 20 bytes
	}
}
void freeTreeMemory(int whichPartition){
	int i,j,k;
	for(i=0; i<2*numspecArr[whichPartition]-1; i++){
		for(j=0; j<numbaseArr[whichPartition]; j++){
			free(treeArr[whichPartition][i].likenc[j]);
			free(treeArr[whichPartition][i].posteriornc[j]);
		}
		free(treeArr[whichPartition][i].like);
		free(treeArr[whichPartition][i].posterior);
		free(treeArr[whichPartition][i].likenc);
		free(treeArr[whichPartition][i].posteriornc);
	}
	free(treeArr[whichPartition]);
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
