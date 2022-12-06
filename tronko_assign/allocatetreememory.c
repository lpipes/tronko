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
