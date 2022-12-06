#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "global.h"

/*void getSequenceInRoot(){
	int i;
	int j;
	type_of_PP minimum;
	int index;
	//seqInRoot = malloc(sizeof(char) * (numbase+1));
	if(!seqInRoot){
		fprintf(stderr, "Out of memory");
		exit(1);
	}
	for(i=0;i<numbase;i++){
		minimum=PP[root][i][0];
		index=0;
		for(j=0; j<4; j++){
			if (minimum > PP[root][i][j]){
				minimum=PP[root][i][j];
				index=j;
			}
			//printf("minimum is %lf i is %d j is %d PP[%d][i][j] is %lf\t",minimum,i,j,root,PP[root][i][j]);
		}
		//printf("root is %d\n",root);
		//printf("index is %d\n",index);
		char base;
		if( index == 0 ){ base='a'; }
			else if( index == 1 ){ base='c'; }
			else if( index == 2 ){ base='g'; }
			else{ base='t'; }
			seqInRoot[i]=base;
	}
	seqInRoot[numbase]='\0';
	printf("seqInRoot is %s\n",seqInRoot);
	//return seqInRoot;
}*/
int getStartPosition(int start, int rootNum, int node, int padding){
	int i, j, k;
	j=0;
	k=0;
	if ( start - padding <= 0 ){
		return 0;
	}
	for(i=0; i<numbaseArr[rootNum]; i++){
		if (j==start){ break; }
		if ( treeArr[rootNum][node].posteriornc[i][0] == 1 ){
			k++;
		}else{
			j++;
		}
	}
	return start+k-1-padding;
}
int getEndPosition(char* s, int rootNum, int node, int start, int padding){
	int i, j, k;
	char *cigar_string = (char*)malloc(MAX_CIGAR*sizeof(char));
	strcpy(cigar_string,s);
	char *copy = (char*)malloc(MAX_CIGAR*sizeof(char));
	strcpy(copy,s);
	char *context = NULL;
	char *res = strtok_r(cigar_string, "MIDHS",&context);
	int cigar [MAX_CIGAR];
	char cigar_chars[MAX_CIGAR];
	int index = 0;
	while(res){
		int from = res - cigar_string + strlen(res);
		int cigar_count = 0;
		sscanf(res, "%d", &cigar_count);
		res = strtok_r( NULL, "MIDHS", &context);
		int to = res != NULL ? res-cigar_string : strlen(copy);
		char cigar_char = '\0';
		sscanf(copy+from, "%c", &cigar_char);
		cigar[index] = cigar_count;
		cigar_chars[index] = cigar_char;
		index++;
	}
	free(copy);
	free(cigar_string);
	int cigar_char_count = index;
	int alignment_length = 0;
	for (i=0; i<cigar_char_count; i++){
		if ( cigar_chars[i] != 'I' ){
			alignment_length += cigar[i];
		}
	}
	j=0;
	k=0;
	alignment_length += padding;
	if ( start + alignment_length >= numbaseArr[rootNum] ){
		return numbaseArr[rootNum];
	}
	for(i=start; i<numbaseArr[rootNum]; i++){
		if (j==alignment_length){ break; }
		if (treeArr[rootNum][node].posteriornc[i][0] == 1){
			k++;
		}else{
			j++;
		}
	}
	return start+k+j;
}
void getSequenceInNode(int rootNum, int node, char *sequenceInNode){
	type_of_PP maximum;
	//char *sequenceInNode;
	int index,i,j;
	//if ((sequenceInNode = (char *)malloc(sizeof(char)*(numbaseArr[rootNum]))) == NULL){
	//	fprintf(stderr, "malloc failed\n");
	//	exit(1);
	//}
	for(i=0;i<numbaseArr[rootNum];i++){	
		//maximum=PP_Arr[rootNum][node][i][0];
		maximum=treeArr[rootNum][node].posteriornc[i][0];
		//minimum=0;
		index=0;
		char base;
		//if (PP_Arr[rootNum][node][i][0] != 1 ){
		if (treeArr[rootNum][node].posteriornc[i][0] != 1 ){
		for(j=0;j<4;j++){
			//if (maximum < PP_Arr[rootNum][node][i][j]){
			if (maximum < treeArr[rootNum][node].posteriornc[i][j]){
				//maximum=PP_Arr[rootNum][node][i][j];
				maximum=treeArr[rootNum][node].posteriornc[i][j];
				index=j;
			}
		}
		if (index==0){ base='A'; }
		else if (index ==1 ){ base='C'; }
		else if (index==2 ){ base='G'; }
		else if (index==3 ){base='T'; }
		else { base='N'; }
		}else{
			base='N';
		}
		sequenceInNode[i]=base;
	}
	sequenceInNode[numbaseArr[rootNum]]='\0';
	//return sequenceInNode;
}
void getSequenceInNodeWithoutNs(int rootNum, int node, char *sequenceInNode, int *positionsInRoot, int start_position, int end_position){
	type_of_PP maximum;
	int index,i,j;
	int count=0;
	//for(i=0; i<numbaseArr[rootNum];i++){
	for(i=start_position; i<end_position; i++){
		//maximum=PP_Arr[rootNum][node][i][0];
		maximum=treeArr[rootNum][node].posteriornc[i][0];
		index=0;
		char base;
		//if (PP_Arr[rootNum][node][i][0] != 1 ){
		if (treeArr[rootNum][node].posteriornc[i][0] != 1 ){
			for(j=0; j<4; j++){
				//if (maximum < PP_Arr[rootNum][node][i][j]){
				if (maximum < treeArr[rootNum][node].posteriornc[i][j]){
					//maximum=PP_Arr[rootNum][node][i][j];
					maximum=treeArr[rootNum][node].posteriornc[i][j];
					index=j;
				}
			}
		}else{
			index=4;
		}
		if (index==0){ base='A'; sequenceInNode[count]=base; positionsInRoot[count]=i; count++;}
		else if (index==1){ base='C'; sequenceInNode[count]=base; positionsInRoot[count]=i; count++;}
		else if (index==2){ base='G'; sequenceInNode[count]=base; positionsInRoot[count]=i; count++;}
		else if (index==3){ base='T'; sequenceInNode[count]=base; positionsInRoot[count]=i; count++;}
	}
	sequenceInNode[count]='\0';
}
