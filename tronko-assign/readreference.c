#include "readreference.h"
int readInXNumberOfLines_fastq(int numberOfLinesToRead, gzFile query_reads, int whichPair, Options opt, int max_query_length, int max_readname_length,int first_iter){
	char* buffer;
	char* query;
	char* reverse;
	int buffer_size = 0;
	if ( max_query_length > max_readname_length ){
		buffer_size = max_query_length + 2;
		buffer = (char *)malloc(sizeof(char)*(max_query_length+2));
	}else{
		buffer_size = max_readname_length + 2;
		buffer = (char *)malloc(sizeof(char)*(max_readname_length+2));
	}
	char* s;
	char seqname[max_readname_length+1];
	int size=0;
	int i=0;
	int iter=0;
	int next=0;
	query = (char *)malloc(sizeof(char)*max_query_length+2);
	reverse = (char *)malloc(sizeof(char)*max_query_length+2);
	for(i=0; i<max_query_length+2; i++){
		query[i] = '\0';
		reverse[i] = '\0';
	}
	while(gzgets(query_reads,buffer,buffer_size)!=NULL){
		s = strtok(buffer,"\n");
		size = strlen(s);
		if ( first_iter == 1 ){
			if ( buffer[0] != '@' ){
				printf("Query reads are not in FASTQ format\n");
				exit(-1);
			}
			first_iter=0;
		}
		if ( buffer[0] == '@' && whichPair==1){
			for(i=1; i<size; i++){
				if (buffer[i]==' '){ buffer[i] = '_'; }
				seqname[i-1]=buffer[i];
			}
			seqname[i-1] = '\0';
			memset(pairedQueryMat->forward_name[iter],'\0',max_readname_length);
			strcpy(pairedQueryMat->forward_name[iter],seqname);
			memset(seqname,'\0',max_readname_length);
			next=1;
		}else if ( buffer[0] == '@' && whichPair==2){
			for(i=1; i<size; i++){
				if(buffer[i]==' '){ buffer[i]='_'; }
				seqname[i-1]=buffer[i];
			}
			seqname[i-1] = '\0';
			char tempname[max_readname_length];
			memset(tempname,'\0',max_readname_length);
			for(i=0; i<size-1; i++){
				if ( pairedQueryMat->forward_name[iter][i] == '1' && pairedQueryMat->forward_name[iter][i-1] == '_'){
					tempname[i] = '2';
				}else{
					tempname[i] = seqname[i];
				}
			}
			tempname[size-1]='\0';
			int skipped=iter;
			/*while(strcmp(tempname,seqname)!=0){
				for(i=0; i<size-1; i++){
					if ( pairedQueryMat->forward_name[skipped][i] == '1' && pairedQueryMat->forward_name[skipped][i-1] == '_'){
						tempname[i] = '2';
					}else{
						tempname[i] = seqname[i];
					}
				}
				tempname[size-1]='\0';
				skipped++;
				if (skipped == numberOfLinesToRead){ break; }
			}*/
			if (skipped == iter){
				memset(pairedQueryMat->reverse_name[iter],'\0',max_readname_length);
				strcpy(pairedQueryMat->reverse_name[iter],seqname);
				memset(seqname,'\0',max_readname_length);
				next=1;
			}else{
				shiftUp(iter,skipped-iter,numberOfLinesToRead);
				next=0;
			}
			/*if ( pairedQueryMat->forward_name[iter+1][0] == '\0' ){
				iter++;
				break;
			}*/
		}else if ( buffer[0] == '@' && whichPair==0){
			for(i=1; i<size; i++){
				if ( buffer[i]==' '){ buffer[i] = '_'; }
				seqname[i-1]=buffer[i];
			}
			seqname[i-1]='\0';
			for(i=0; i<max_readname_length; i++){
				singleQueryMat->name[iter][i]='\0';
			}
			strcpy(singleQueryMat->name[iter],seqname);
			for (i=0; i<size-1; i++){
				seqname[i]='\0';
			}
			next=1;
		}else if (next==1){
			for(i=0; i<size; i++){
				query[i]=toupper(buffer[i]);
			}
			query[size]='\0';
			if (whichPair == 0){
				if (opt.reverse_single_read != 1){
					strcpy(singleQueryMat->queryMat[iter],query);
				}else{
					getReverseComplement(query,reverse,max_query_length);
					strcpy(singleQueryMat->queryMat[iter],reverse);
				}
			}
			if (whichPair == 1){
				strcpy(pairedQueryMat->query1Mat[iter],query);
			}
			if (whichPair == 2){
				if (opt.reverse_second_of_paired_read != 1){
					strcpy(pairedQueryMat->query2Mat[iter],query);
				}else{
					getReverseComplement(query,reverse,max_query_length);
					strcpy(pairedQueryMat->query2Mat[iter],reverse);
					for(i=0; i<size; i++){
						reverse[i] = '\0';
					}
				}
			}
			for(i=0; i<size; i++){
				query[i] = '\0';
			}
			iter++;
			next=0;
			if(iter==numberOfLinesToRead){ break; }
		}
	}
	free(buffer);
	free(query);
	free(reverse);
	return iter;
}
//this functions returns the number of 'reads'. This is half the number of actual lines.
//if eof it returns zero.
int readInXNumberOfLines(int numberOfLinesToRead, gzFile query_reads, int whichPair, Options opt, int max_query_length, int max_readname_length){
	char* buffer;
	char* query;
	char* reverse;
	int buffer_size = 0;
	if ( max_query_length > max_readname_length ){
		buffer_size = max_query_length +2;
		buffer = (char *)malloc(sizeof(char)*(max_query_length+2));
	}else{
		buffer_size = max_readname_length +2;
		buffer = (char *)malloc(sizeof(char)*(max_readname_length+2));
	}
	char seqname[max_readname_length];
	int size;
	char *s;
	int i;
	int iter=0;
	int next=0;
	query = (char *)malloc(sizeof(char)*max_query_length+2);
	reverse = (char *)malloc(sizeof(char)*max_query_length+2);
	for(i=0; i<max_query_length+2; i++){
		query[i]='\0';
		reverse[i]='\0';
	}
	int first_iter=1;
	while(gzgets(query_reads,buffer,buffer_size)!=NULL){
		s = strtok(buffer,"\n");
		size = strlen(s);
		if (first_iter==1){
			if ( buffer[0] != '>' ){
				printf("Query reads are not in FASTA format. Try specifying -q if using FASTQ reads.\n");
				exit(-1);
			}
			first_iter=0;
		}
		if ( buffer[0] == '>' && whichPair==1 ){
			for(i=1; i<size; i++){
				if ( buffer[i]==' '){ buffer[i] = '_'; }
				seqname[i-1]=buffer[i];
			}
			seqname[i-1] = '\0';
			memset(pairedQueryMat->forward_name[iter],'\0',max_readname_length);
			strcpy(pairedQueryMat->forward_name[iter],seqname);
			memset(seqname,'\0',max_readname_length);
			next=1;
		}else if (buffer[0] == '>' && whichPair==2) {
			for(i=1; i<size; i++){
				if ( buffer[i]==' '){ buffer[i] = '_'; }
				seqname[i-1]=buffer[i];
			}
			seqname[i-1] = '\0';
			char tempname[max_readname_length];
			memset(tempname,'\0',max_readname_length);
			for(i=0; i<size-1; i++){
				if (pairedQueryMat->forward_name[iter][i] == '1' && pairedQueryMat->forward_name[iter][i-1] == '_'){
					tempname[i] ='2';
				}else{
					tempname[i] = seqname[i];
				}
			}
			tempname[size-1]='\0';
			int skipped=iter;
			/*while(strcmp(tempname,seqname)!=0){
				for(i=0; i<size-1; i++){
					if ( pairedQueryMat->forward_name[skipped][i] == '1' && pairedQueryMat->forward_name[skipped][i-1] == '_'){
						tempname[i] = '2';
					}else{
						tempname[i] = seqname[i];
					}
					tempname[size-1]='\0';
					skipped++;
					if (skipped == numberOfLinesToRead){ break; }
				}
			}*/
			if (skipped == iter){
				memset(pairedQueryMat->reverse_name[iter],'\0',max_readname_length);
				strcpy(pairedQueryMat->reverse_name[iter],seqname);
				memset(seqname,'\0',max_readname_length);
				next=1;
			}else{
				shiftUp(iter,skipped-iter,numberOfLinesToRead);
				next=0;
			}
			/*if ( pairedQueryMat->forward_name[iter+1][0] == '\0' ){
				iter++;
				break;
			}*/
			/*for(i=0; i < MAXREADNAME; i++){
				pairedQueryMat->reverse_name[iter][i]='\0';
			}
			strcpy(pairedQueryMat->reverse_name[iter],seqname);
			for (i=0; i<size-1; i++){
				seqname[i]='\0';
			}*/
		}else if (buffer[0] == '>' && whichPair==0){
			for(i=1; i<size; i++){
				if ( buffer[i]==' '){ buffer[i] = '_'; }
				seqname[i-1]=buffer[i];
			}
			seqname[i-1] = '\0';
			for(i=0; i<max_readname_length; i++){
				singleQueryMat->name[iter][i]='\0';
			}
			strcpy(singleQueryMat->name[iter],seqname);
			for (i=0; i<size-1; i++){
				seqname[i]='\0';
			}
			next=1;
		}else{
			//char query[size];
			//char* query;
			//query = (char *)malloc(sizeof(char)*size);
			for(i=0; i<size; i++){
				query[i]=toupper(buffer[i]);
			}
			query[size]='\0';
			if ( whichPair == 0 ){
				if (opt.reverse_single_read != 1){
					strcpy(singleQueryMat->queryMat[iter],query);
				}else{
					getReverseComplement(query,reverse,max_query_length);
					strcpy(singleQueryMat->queryMat[iter],reverse);
				}
			}
			if ( whichPair == 1 ){
				strcpy(pairedQueryMat->query1Mat[iter],query);
			}
			if ( whichPair == 2 ){
				//char* reverse;
				//reverse = (char *)malloc(sizeof(char)*size);
				if (opt.reverse_second_of_paired_read != 1){
					strcpy(pairedQueryMat->query2Mat[iter],query);
				}else{
					getReverseComplement(query,reverse,max_query_length+2);
					strcpy(pairedQueryMat->query2Mat[iter],reverse);
				//free(reverse);
				//strcpy(pairedQueryMat->query2Mat[iter],query);
					for(i=0; i<size; i++){
						reverse[i] = '\0';
					}
				}
			}
			for(i=0; i<size; i++){
				query[i] = '\0';
			}
			//free(query);
			iter++;
			next=0;
			if(iter==numberOfLinesToRead)
				break;
			}
	}
	free(buffer);
	free(query);
	free(reverse);
	return iter;
}
void shiftUp(int iter, int jump, int numberOfLinesToRead){
	int i;
	for(i=iter; i<numberOfLinesToRead-jump; i++){
		strcpy(pairedQueryMat->forward_name[i],pairedQueryMat->forward_name[i+jump]);
		strcpy(pairedQueryMat->query1Mat[i],pairedQueryMat->query1Mat[i+jump]);
	}
	for(i=numberOfLinesToRead-jump; i<numberOfLinesToRead; i++){
		memset(pairedQueryMat->forward_name[i],'\0',MAXREADNAME);
	}
}
int readReferenceTree( gzFile referenceTree, int* name_specs){
	char buffer[BUFFER_SIZE];
	char acc_name[30];
	int up[2], taxIndex[2];
	char *s;
	int max_lineTaxonomy, max_tax_name, max_nodename, treeNumber, nodeNumber, down, depth, success, i, j, k, firstIter, numberOfTrees;
	firstIter = 1;
	char* refTreeFlag = buffer;
	while (refTreeFlag != NULL ){
		if ( firstIter == 1 ){
			refTreeFlag = gzgets(referenceTree,buffer,BUFFER_SIZE);
			if(refTreeFlag == NULL) {
				break;
			}
			s = strtok(buffer, "\n");
			if ( s == NULL ){
				success = 0;
			}else{
				success = sscanf(s, "%d", &numberOfTrees);
			}
			refTreeFlag = gzgets(referenceTree,buffer,BUFFER_SIZE);
			if(refTreeFlag == NULL) {
				break;
			}
			s = strtok(buffer, "\n");
			if ( s == NULL ){
				success = 0;
			}else{
				success = sscanf(s, "%d", &max_nodename);
				name_specs[0]=max_nodename;
			}
			refTreeFlag = gzgets(referenceTree,buffer,BUFFER_SIZE);
			s = strtok(buffer, "\n");
			if ( s == NULL ){
				success = 0;
			}else{
				success = sscanf(s, "%d", &max_tax_name);
				name_specs[1]=max_tax_name;
			}
			refTreeFlag = gzgets(referenceTree,buffer,BUFFER_SIZE);
			s = strtok(buffer, "\n");
			if ( s == NULL ){
				success = 0;
			}else{
				success = sscanf(s, "%d", &max_lineTaxonomy);
				name_specs[2]=max_lineTaxonomy;
			}
			numbaseArr = (int*)malloc(numberOfTrees*(sizeof(int)));
			rootArr = (int*)malloc(numberOfTrees*(sizeof(int)));
			numspecArr= (int*)malloc(numberOfTrees*(sizeof(int)));
			for(i=0; i<numberOfTrees; i++){	
				refTreeFlag = gzgets(referenceTree,buffer,BUFFER_SIZE);
				if ( refTreeFlag == NULL ){
					break;
				}
				s = strtok(buffer, "\n");
				if ( s == NULL ){
					success = 0;
				}else{
					success = sscanf(s, "%d\t%d\t%d",&(numbaseArr[i]),&(rootArr[i]),&(numspecArr[i]));
				}	
			}
			//for(i=0;i<numberOfTrees;i++){
			//	printf("Tree %d Numbase: %d, Root: %d, Numspec %d\n",i,numbaseArr[i],rootArr[i],numspecArr[i]);
			//}
			allocateMemoryForTaxArr(numberOfTrees,max_tax_name);
			for(i=0;i<numberOfTrees;i++){
				for(j=0; j<numspecArr[i]; j++){
					refTreeFlag = gzgets(referenceTree,buffer,BUFFER_SIZE);
					if(refTreeFlag == NULL) {
						break;
					}
					s = strtok(buffer,";\n");
					taxonomyArr[i][j][0]=strcpy(taxonomyArr[i][j][0],s);
					for(k=1;k<7;k++){
						s = strtok(NULL,";\n");
						taxonomyArr[i][j][k]=strcpy(taxonomyArr[i][j][k],s);
					}
				}
			}
			treeArr = malloc(numberOfTrees*sizeof(struct node *));
			for (i=0; i<numberOfTrees; i++){
				allocateTreeArrMemory(i,max_nodename);
			}
			//allocatetreememory_for_nucleotide_Arr(numberOfTrees);
			firstIter=0;
		}
		refTreeFlag = gzgets(referenceTree,buffer,BUFFER_SIZE);
		if (refTreeFlag == NULL){
			break;
		}
		s = strtok(buffer, "\n");
		sscanf(s, "%d %d %d %d %d %d %d %d %s",&treeNumber, &nodeNumber, &(up[0]), &(up[1]), &down, &depth, &(taxIndex[0]), &(taxIndex[1]), acc_name);
		treeArr[treeNumber][nodeNumber].up[0] = up[0];
		treeArr[treeNumber][nodeNumber].up[1] = up[1];
		treeArr[treeNumber][nodeNumber].down = down;
		treeArr[treeNumber][nodeNumber].depth = depth;
		treeArr[treeNumber][nodeNumber].taxIndex[0] = taxIndex[0];
		treeArr[treeNumber][nodeNumber].taxIndex[1] = taxIndex[1];
		if ( up[0] == -1 && up[1] == -1){
			strcpy(treeArr[treeNumber][nodeNumber].name,acc_name);
		}
		for(i=0; i<numbaseArr[treeNumber]; i++){
			refTreeFlag = gzgets(referenceTree,buffer,BUFFER_SIZE);
			if(refTreeFlag == NULL) {
				break;
			}
			for (j=0; j<4; j++){
				if (j==0){
					s = strtok(buffer,"\t");
				}else{
					s = strtok(NULL,"\t");
				}
				if ( s==NULL ){
					success = 0;
				}else{ 
					success = sscanf(s,"%lf",&(treeArr[treeNumber][nodeNumber].posteriornc[i][j]));
				}
			}
		}
	}
	return numberOfTrees;
}
int setNumbase_setNumspec(int numberOfPartitions, int* specs){
	int i, maxNumbase, maxNumSpec;
	maxNumbase=0;
	maxNumSpec=0;
	int numspec_total=0;
	for (i=0; i<numberOfPartitions; i++){
		if (maxNumbase < numbaseArr[i]){
			maxNumbase = numbaseArr[i];
		}
	}
	specs[1] = maxNumbase;
	for( i=0; i<numberOfPartitions; i++){
		numspec_total += numspecArr[i];
		if (numspecArr[i] > maxNumSpec){
			maxNumSpec = numspecArr[i];
		}
	}
	specs[0] = maxNumSpec;
	return numspec_total;
}
void find_specs_for_reads(int* specs, gzFile file, int format){
	int max_name_length = specs[0];
	int max_query_length = specs[1];
	char *buffer = (char *)malloc(sizeof(char)*FASTA_MAXLINE);
	char *s;
	int size = 0;
	while(gzgets(file,buffer,FASTA_MAXLINE) != NULL){
		s = strtok(buffer,"\n");
		if (buffer[0] == '>' || buffer[0] == '@' ){
			size = strlen(s);
			if (max_name_length < size ){
				max_name_length = size;
			}
		}else{
			size = strlen(s);
			if ( max_query_length < size ){
				max_query_length = size;
			}
		}
	}
	specs[0] = max_name_length;
	specs[1] = max_query_length;
	free(buffer);
}
