#include "printAlignments.h"

void createNewFile( char *referenceName, char *alignments_directory, alignment_t *aln, char *readname){
	FILE *alignmentFile;
	char fileName[1000];
	snprintf(fileName,1000,"%s/%s.fasta",alignments_directory, referenceName);	
	alignmentFile = fopen(fileName, "w");
	if (alignmentFile == NULL){ printf("Error!"); exit(1); }
	fprintf(alignmentFile,">%s\n",referenceName);
	int i;
	for(i=0; aln->result_a[i]!='\0';i++){
		if (aln->result_a[i]!='-'){
			if (aln->result_a[i]=='N'){
				fprintf(alignmentFile,"-");
			}else{
				fprintf(alignmentFile,"%c",aln->result_a[i]);
			}
		}
	}
	fprintf(alignmentFile,"\n");
	fprintf(alignmentFile,">%s\n",readname);
	for(i=0; aln->result_b[i]!='\0';i++){
		if ( aln->result_a[i]!='-'){
			if (aln->result_b[i]=='N'){
				fprintf(alignmentFile,"-");
			}else{
				fprintf(alignmentFile,"%c",aln->result_b[i]);
			}
		}
	}
	fprintf(alignmentFile,"\n");
	fclose(alignmentFile);
}
void printToFile (char *referenceName, char *alignments_directory, alignment_t *aln, char *readname){
	FILE *alignmentFile;
	char fileName[1000];
	snprintf(fileName,1000,"%s/%s.fasta",alignments_directory, referenceName);
	alignmentFile = fopen(fileName, "a");
	 if (alignmentFile == NULL){ printf("Error!"); exit(1); }
	fprintf(alignmentFile,">%s\n",readname);
	int i;
	for(i=0; aln->result_b[i]!='\0';i++){
		if (aln->result_a[i]!='-'){
			if (aln->result_b[i]=='N'){
				fprintf(alignmentFile,"-");
			}else{
				fprintf(alignmentFile,"%c",aln->result_b[i]);
			}
		}
	}
	fprintf(alignmentFile,"\n");
	fclose(alignmentFile);
}
void createNewFile2( char *referenceName, char *alignments_directory, alignment_t *aln, char *readname, int length_of_read){
	FILE *alignmentFile;
	char fileName[1000];
	snprintf(fileName,1000,"%s/%s.fasta",alignments_directory, referenceName);	
	alignmentFile = fopen(fileName, "w");
	if (alignmentFile == NULL){ printf("Error!"); exit(1); }
	/*fprintf(alignmentFile,">%s\n",referenceName);
	int i;
	for(i=0; aln->result_a[i]!='\0';i++){
		if (aln->result_a[i]!='-'){
			if (aln->result_a[i]=='N'){
				fprintf(alignmentFile,"-");
			}else{
				fprintf(alignmentFile,"%c",aln->result_a[i]);
			}
		}
	}
	fprintf(alignmentFile,"\n");*/
	fprintf(alignmentFile,"%s\t%s\t%d\n",readname,referenceName,length_of_read);
	int counter=0;
	int iter=0;
	int offset=0;
	//b is query
	//a is reference
	while(aln->result_b[counter]!='\0'){
		if ( aln->result_b[counter]!='-'){
			if (aln->result_a[counter]!='-'){
				if ( iter != 0 ){
					fprintf(alignmentFile,";%d",counter-offset);
				}else{
					fprintf(alignmentFile,"%d",counter-offset);
					iter=1;
				}
			}else{
				if (iter != 0){
					fprintf(alignmentFile,";-1");
				}else{
					fprintf(alignmentFile,"-1");
					iter=1;
				}
				offset++;
			}
		}
		counter++;
	}
	fprintf(alignmentFile,"\n");
	fclose(alignmentFile);
}
void printToFile2 (char *referenceName, char *alignments_directory, alignment_t *aln, char *readname, int length_of_read, int* positionsInRoot){
	FILE *alignmentFile;
	char fileName[1000];
	//snprintf(fileName,1000,"%s/%s.fasta",alignments_directory, referenceName);
	snprintf(fileName,1000,"%s",alignments_directory);
	alignmentFile = fopen(fileName, "a");
	 if (alignmentFile == NULL){ printf("Error!"); exit(1); }
	fprintf(alignmentFile,">%s\t%s\t%d\n",readname,referenceName,length_of_read);
	int i;
	int counter=0;
	int iter=0;
	int offset=0;
	while(aln->result_b[counter]!='\0'){
		if (aln->result_b[counter]!='-'){
			if (aln->result_a[counter]!='-'){
				if ( iter!= 0){
					//fprintf(alignmentFile,";%d",counter-offset);
					fprintf(alignmentFile,";%d",positionsInRoot[counter-offset]);
				}else{
					//fprintf(alignmentFile,"%d",counter-offset);
					fprintf(alignmentFile,"%d",positionsInRoot[counter-offset]);
					iter=1;
				}
			}else{
				if (iter != 0){
					fprintf(alignmentFile,";-1");
				}else{
					fprintf(alignmentFile,"-1");
					iter=1;
				}
				offset++;
			}
		}
		counter++;
	}
	fprintf(alignmentFile,";\n");
	fclose(alignmentFile);
}
void printToFile_WFA(char *referenceName, char *alignments_directory, char* const pattern_alg, char* const text_alg, char* readname, int length_of_read, int* positionsInRoot){
	FILE *alignmentFile;
	char fileName[1000];
	snprintf(fileName,1000,"%s",alignments_directory);
	alignmentFile = fopen(fileName, "a");
	if (alignmentFile == NULL){ printf("Error No Alignment File!"); exit(1); }
	fprintf(alignmentFile,">%s\t%s\t%d\n",readname,referenceName,length_of_read);
	int i;
	int counter=0;
	int iter=0;
	int offset=0;
	while(text_alg[counter]!='\0'){
		if (text_alg[counter]!='-'){
			if (pattern_alg[counter]!='-'){
				if (iter != 0){
					fprintf(alignmentFile,";%d",positionsInRoot[counter-offset]);
				}else{
					fprintf(alignmentFile,"%d",positionsInRoot[counter-offset]);
					iter=1;
				}
			}else{
				if (iter != 0){
					fprintf(alignmentFile,";-1");
				}else{
					fprintf(alignmentFile,"-1");
					iter=1;
				}
				offset++;
			}
		}
		counter++;
	}
	fprintf(alignmentFile,";\n");
	fclose(alignmentFile);
}	
void printLeaveSeqsToFile(FILE *output, int node, int whichTree, int numbase){
	int child0 = treeArr[whichTree][node].up[0];
	int child1 = treeArr[whichTree][node].up[1];
	if (child0 == -1 && child1 == -1){
		char *leaf_sequence = (char *)malloc((numbaseArr[whichTree]+1)*sizeof(char));
		int i;
		for(i=0;i<numbaseArr[whichTree]+1;i++){
			leaf_sequence[i]='\0';
		}
		getSequenceInNode(whichTree,node,leaf_sequence);
		fprintf(output,">%s\n",treeArr[whichTree][node].name);
		fprintf(output,"%s\n",leaf_sequence);
		free(leaf_sequence);
		return;
	}else{
		printLeaveSeqsToFile(output,child0,whichTree,numbase);
		printLeaveSeqsToFile(output,child1,whichTree,numbase);
	}
}
void printTreesNumbersToFile(FILE *output, int node, int whichTree){
	int child0 = treeArr[whichTree][node].up[0];
	int child1 = treeArr[whichTree][node].up[1];
	if (child0 == -1 && child1 == -1){
		fprintf(output,"%s\t%d\t%d\n",treeArr[whichTree][node].name,whichTree,node);		
		return;
	}else{
		printTreesNumbersToFile(output,child0,whichTree);
		printTreesNumbersToFile(output,child1,whichTree);
	}
}
