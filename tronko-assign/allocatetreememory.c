#include "allocatetreememory.h"

void allocatetreememory_for_nucleotide_Arr(int numberOfTrees){
	int i, j, k, l;
	for(i=0; i<numberOfTrees; i++){
		for (j=0; j<(numspecArr[i]*2-1); j++){
			treeArr[i][j].posteriornc = (double**)malloc(numbaseArr[i]*(sizeof(double *)));
			for (k=0; k<numbaseArr[i]; k++){
				treeArr[i][j].posteriornc[k] = (double*)malloc(4*(sizeof(double)));
				for(l=0; l<4; l++){
					treeArr[i][j].posteriornc[k][l]=0;
				}		
			}
		}
	}
}
void allocateTreeArrMemory(int whichPartition, int max_nodename){
	int i,j,k,l;
	treeArr[whichPartition]=malloc((numspecArr[whichPartition]*2-1)*(sizeof(struct node)));
	for (i=0; i<(numspecArr[whichPartition]*2-1); i++){
		treeArr[whichPartition][i].up[0] = -2;  // -2 is uninitialized / NULL value
		treeArr[whichPartition][i].up[1] = -2;
		treeArr[whichPartition][i].down = -2;
		treeArr[whichPartition][i].depth = -2;
		//treeArr[whichPartition][i].posteriornc = NULL;
		treeArr[whichPartition][i].posteriornc = (double**)malloc(numbaseArr[whichPartition]*sizeof(double*));
		for ( k=0; k<numbaseArr[whichPartition]; k++){
			treeArr[whichPartition][i].posteriornc[k] = (double*)malloc(4*(sizeof(double)));
			for(l=0; l<4; l++){
				treeArr[whichPartition][i].posteriornc[k][l]=0.0;
			}
		}
	}
	for(i=numspecArr[whichPartition]-1; i<(numspecArr[whichPartition]*2-1); i++){
		treeArr[whichPartition][i].name = (char *)malloc(max_nodename*sizeof(char));
		for(j=0; j<max_nodename; j++){
			treeArr[whichPartition][i].name[j]='\0';
		}
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
void allocateMemoryForTaxArr(int whichPartitions, int max_tax_name_len){
	int phylogenyLevels = 7;
	max_tax_name_len++;
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
void getReverseComplement(char *read, char *reverseComplement, int max_query_length){
	static const unsigned char basemap[256] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O', 'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z', 91, 92, 93, 94, 95, 64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o', 'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,  160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255 };
	int length=0;
	int i;
	for (i=0; i<max_query_length; i++){
		if ( read[i] != '\0' ){
			length++;
		}else{
			break;
		}
	}
	int j=0;
	for (i=length-1; i>0; i--){
		reverseComplement[j] = basemap[(int)read[i]];
		j++;
	}
	reverseComplement[length]='\0';
}
