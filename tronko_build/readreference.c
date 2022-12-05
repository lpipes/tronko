#include "readreference.h"
void readFilesInDir(char *directory, int number_of_partitions, partition_files* pf){
	struct dirent *de, *de2;
	DIR *dr = opendir(directory);
	regex_t regex1;
	regex_t regex2;
	int reti;
	char msgbuf[100];
	int i, j;
	j=0;
	char name[300];
	char name2[300];
	char name3[300];
	reti = regcomp(&regex1, "\_MSA\.fasta$", 0);
	if (reti){
		printf(stderr, "Could not compile regex\n");
		exit(1);
	}
	if (dr==NULL){
		printf("Could not open directory for reads");
		return 0;
	}
	while((de=readdir(dr)) != NULL){
		reti = regexec(&regex1, de->d_name, 0, NULL, 0);
		if (!reti){
			strcpy(pf->msa_files[j],de->d_name);
			for(i=0; i<300; i++){
				name[i]='\0';
			}
			for(i=0; i<(int)strlen(de->d_name)-10; i++){
				name[i] = de->d_name[i];
			}
			name[i] = '\0';
			strcpy(name2,name);
			strncat(name2,"_taxonomy.txt",250);
			strcpy(pf->tax_files[j],name2);
			strcpy(name3,"RAxML_bestTree.");
			strncat(name3,name,250);
			strncat(name3,".reroot",250);
			strcpy(pf->tree_files[j],name3);
			j++;
		}else if (reti == REG_NOMATCH){
		}else{
			regerror(reti, &regex1, msgbuf, sizeof(msgbuf));
			printf(stderr, "Regex match failed: %s\n", msgbuf);
			exit(1);
		}
	}
	closedir(dr);
	return pf;
}
int readReferenceTree(gzFile referenceTree){
	char buffer[BUFFER_SIZE];
	char acc_name[30];
	int up[2], taxIndex[2];
	char *s;
	int max_tax_name, max_nodename, treeNumber, nodeNumber, down, depth, success, i, j, k, firstIter, numberOfTrees;
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
			if ( s ==NULL){
				success = 0;
			}else{
				success = sscanf(s, "%d", &max_nodename);
			}
			refTreeFlag = gzgets(referenceTree,buffer,BUFFER_SIZE);
			if(refTreeFlag == NULL) {
				break;
			}
			s = strtok(buffer, "\n");
			if ( s == NULL ){
				success = 0;
			}else{
				success = sscanf(s, "%d", &max_tax_name);
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
			for(i=0;i<numberOfTrees;i++){
				printf("Tree %d Numbase: %d, Root: %d, Numspec %d\n",i,numbaseArr[i],rootArr[i],numspecArr[i]);
			}
			allocateMemoryForTaxArr(numberOfTrees);
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
			treeArr = malloc(numberOfTrees*sizeof(node *));
			for (i=0 ; i<numberOfTrees; i++){
				allocateTreeArrMemory(i,max_nodename);
			}
			allocatetreememory_for_nucleotide_Arr(numberOfTrees);
			firstIter=0;
		}
		refTreeFlag = gzgets(referenceTree,buffer,BUFFER_SIZE);
		if (refTreeFlag == NULL){
			break;
		}
		s = strtok(buffer, "\n");
		sscanf(s, "%d %d %d %d %d %d %d %d %d %s",&treeNumber, &nodeNumber, &(up[0]), &(up[1]), &down, &depth, &(taxIndex[0]), &(taxIndex[1]), acc_name);
		treeArr[treeNumber][nodeNumber].up[0] = up[0];
		treeArr[treeNumber][nodeNumber].up[1] = up[1];
		treeArr[treeNumber][nodeNumber].down = down;
		treeArr[treeNumber][nodeNumber].depth = depth;
		treeArr[treeNumber][nodeNumber].taxIndex[0] = taxIndex[0];
		treeArr[treeNumber][nodeNumber].taxIndex[1] = taxIndex[1];
		if ( up[0] != -1 && up[1] != -1){
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
