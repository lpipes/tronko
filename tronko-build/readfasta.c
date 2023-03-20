#include "readfasta.h"

/*some old code for reading a fasta file and storing DNA as ints*/
void setNumspec(FILE* infile, int* specs){
	char buffer[FASTA_MAXLINE];
	int countLeafNodes = 0;
	int max_numbase=0;
	int max_nodename=0;
	int iter=0;
	int i=0;
	while( fgets(buffer,FASTA_MAXLINE,infile) != NULL ){
		if ( buffer[0] == '>' ){
			countLeafNodes++;
			if ( max_nodename < strlen(buffer) ){
				max_nodename = strlen(buffer);
			}
		}else if ( iter == 0 ){
			int length=0;
			for(i=0; buffer[i]!='\n'; i++){
				length++;
			}
			max_numbase=length;
			iter=1;
		}
	}
	if (specs[0] < countLeafNodes){
		specs[0] = countLeafNodes;
	}
	if (specs[1] < max_nodename){
		specs[1] = max_nodename+1;
	}
	if (specs[2] < max_numbase){
		specs[2] = max_numbase;
	}
}
int setNumspecArr(FILE *partitionsFile){
	char buffer[FASTA_MAXLINE];
	int countLeafNodes = 0;
	while( fgets(buffer,FASTA_MAXLINE,partitionsFile) != NULL ){
		if ( buffer[0] == '>' ){
			countLeafNodes++;
		}
	}
	return countLeafNodes;
}
/*void readSeqArr_UsePartitions(FILE *partitionsFile, int whichPartitions, int*** seqArr_heap,char*** nodeIDsArr_heap){
	char buffer[FASTA_MAXLINE];
	int i, j, m, k=0, row=0;
	char c;
	int size;
	char nodename[MAX_NODENAME];
	char *s;
	int firstIter=1;
	int index=0;
	//nodeIDsArr[whichPartitions]=(char **)malloc(sizeof(char *)*numspecArr[whichPartitions]);
	unsigned int object_size = 3*(numspecArr[whichPartitions]*sizeof(char **));
	nodeIDsArr_heap[whichPartitions]=mmap(NULL, object_size, PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
	if (madvise(nodeIDsArr_heap[whichPartitions], object_size, MADV_DONTFORK)==-1){
		perror("madvise error");
	}
	for(i=0;i<numspecArr[whichPartitions];i++){
		unsigned int object_size2 = 3*(MAX_NODENAME*sizeof(char *));
		nodeIDsArr_heap[whichPartitions][i]=mmap(NULL, object_size2, PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
		if (madvise(nodeIDsArr_heap[whichPartitions][i], object_size2, MADV_DONTFORK)==-1){
			perror("madvise error");
		}
		//nodeIDsArr[whichPartitions][i] = (char *)malloc(sizeof(char)*MAX_NODENAME);
	}
	while( fgets(buffer,FASTA_MAXLINE,partitionsFile) != NULL){
		s = strtok(buffer,"\n");
		size = strlen(s);
		if ( buffer[0] == '>'){
			if ( size > MAX_NODENAME ){ size = MAX_NODENAME; }
			for(i=1; i<size; i++){
				nodename[i-1]=buffer[i];
			}
			nodename[size-1]='\0';
			strcpy(nodeIDsArr_heap[whichPartitions][index], nodename);
			index++;
		}else{
			if ( firstIter==1 ){ 
				numbaseArr[whichPartitions]=size;
				unsigned int object_size = 3*(numspecArr[whichPartitions]*(sizeof(int *)));
				//seqArr_heap[whichPartitions] = mmap(NULL, object_size, PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
				//if (madvise(seqArr_heap[whichPartitions], object_size, MADV_DONTFORK)==-1){
				//	perror("madvise error");
				//	exit(0);
				//}
				seqArr_heap[whichPartitions] = (int **)malloc(numspecArr[whichPartitions]*(sizeof(int *)));
				//seqArr[whichPartitions] =(int **)malloc(numspecArr[whichPartitions]*(sizeof(int *)));
				for (i=0; i<numspecArr[whichPartitions]; i++){
					unsigned int object_size2 = 4*(numbaseArr[whichPartitions]*(sizeof(int)));
					//seqArr_heap[whichPartitions][i]=mmap(NULL, object_size2, PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
					//if(madvise(seqArr_heap[whichPartitions][i], object_size2, MADV_DONTFORK==-1)){
					//	perror("madvise error");
					//	exit(0);
					//}
					seqArr_heap[whichPartitions][i] = (int **)malloc(numbaseArr[whichPartitions]*(sizeof(int)));
					//seqArr[whichPartitions][i]=(int *)malloc(numbaseArr[whichPartitions]*(sizeof(int)));
				}
				firstIter=0;
				i=0;
			}
			for( m=0; m<numbaseArr[whichPartitions]; m++){
				c=tolower(buffer[m]);
				if ((c!=' ')&&(c!='\t')) {
					if (c=='a') seqArr_heap[whichPartitions][row][m] = 0;
					else if (c=='c') seqArr_heap[whichPartitions][row][m] = 1;
					else if (c=='g') seqArr_heap[whichPartitions][row][m] = 2;
					else if (c=='t') seqArr_heap[whichPartitions][row][m] = 3;
					else if (c=='n'||c=='-'||c=='~') seqArr_heap[whichPartitions][row][m] = 4;
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
	if (numspecArr[whichPartitions] < 3) {printf("This is for more than two seq.s only!\n"); exit(-1);}
}*/
void readseq(FILE* infile, int max_nodename, struct masterArr *master){
	char buffer[FASTA_MAXLINE];
	int i, j, test, m, k=0, row=0;
	for(i=0; i<FASTA_MAXLINE; i++){
		buffer[i] = '\0';
	}
	char c;
	int size;
	char nodename[max_nodename];
	char *s;
	int firstIter=1;
	int index=0;
	int tmp1=0;
	while( fgets(buffer,FASTA_MAXLINE,infile) != NULL){
		char* buffer2;
		s = strtok_r(buffer,"\n\0",&buffer2);
		size = strlen(buffer);
		if ( buffer[0] == '>'){
			if ( size > max_nodename ){ size = max_nodename; }
			for(i=1; i<size; i++){
				nodename[i-1]=buffer[i];
			}
			nodename[size-1]='\0';
			strcpy(master->names[index], nodename);
			index++;
		}else{
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
	int tmp2=0;
	if (master->numspec< 3) {printf("This is for more than two seqs only!\n"); exit(-1);}
}

void readfasta(FILE *infile){
	char buffer[FASTA_MAXLINE];
	int i;
	char *seqname;
	seqname = malloc(sizeof(char) * 20);
	char *s;
	int size;
	char *locQuery;
	locQuery = malloc(sizeof(char) *250);
	//int alength=150;
	//infile = openFasta("Todiramphus_chloris_COI_reads.fasta");
	//infile = openFasta("testreads.fasta");
	while( fgets(buffer,FASTA_MAXLINE,infile) != NULL){
		s = strtok(buffer,"\n");
		size = strlen(s);
		if ( buffer[0] == '>'){
			if ( size > 20 ){ size = 20; }
			for(i=1; i<size; i++){
				seqname[i-1]=buffer[i];
			}
		}else{
			//for(i=0; i<size; i++){
				//printf("char is %c\n",buffer[i]);
			//	locQuery[i]=buffer[i];
			strcpy(locQuery,buffer);
			locQuery[strlen(buffer)]='\0';
		}
	}
	printf("seqname is %s\n",seqname);
	fclose(infile);
}
