#include "readreference.h"
int compare_strings(const void* a, const void* b) {
    return strcmp(*(const char**)a, *(const char**)b);
}
#include <ctype.h>
#include <string.h>

int compare_natural(const void *a, const void *b) {
    const char *ia = *(const char **)a;
    const char *ib = *(const char **)b;
    while (*ia && *ib) {
        if (isdigit(*ia) && isdigit(*ib)) {
            // Compare the numbers
            long na = strtol(ia, (char **)&ia, 10);
            long nb = strtol(ib, (char **)&ib, 10);
            if (na < nb) {
                return -1;
            } else if (na > nb) {
                return 1;
            }
        } else {
            // Compare the characters
            if (*ia < *ib) {
                return -1;
            } else if (*ia > *ib) {
                return 1;
            }
            ++ia;
            ++ib;
        }
    }
    if (*ia) {
        return 1;
    } else if (*ib) {
        return -1;
    } else {
        return 0;
    }
}

void readFilesInDir(char *directory, int number_of_partitions, partition_files* pf) {
    struct dirent *de;
    DIR *dr = opendir(directory);
    regex_t regex1;
    char msgbuf[100];
    int i, j;
    j = 0;
    char name[300];
    char name2[300];
    char name3[300];
    char **msa_files = (malloc)(number_of_partitions*sizeof(char*)); // intermediate array to store filenames
	for(i=0; i< number_of_partitions; i++){
		msa_files[i] = (malloc)(300*sizeof(char));
		for(j=0; j<300; j++){
			msa_files[i][j] = '\0';
		}
	}
    int reti = regcomp(&regex1, "_MSA\\.fasta$", 0);
    if (reti) {
        fprintf(stderr, "Could not compile regex\n");
        exit(-1);
    }

    if (dr == NULL) {
        printf("Could not open directory for reads");
        exit(-1);
    }

    int file_count = 0;
    while ((de = readdir(dr)) != NULL) {
        reti = regexec(&regex1, de->d_name, 0, NULL, 0);
        if (!reti) {
            msa_files[file_count] = strdup(de->d_name); // store filenames in the array
            file_count++;
        } else if (reti != REG_NOMATCH) {
            regerror(reti, &regex1, msgbuf, sizeof(msgbuf));
            fprintf(stderr, "Regex match failed: %s\n", msgbuf);
            exit(1);
        }
    }

    closedir(dr);

    // Sort the msa_files
    qsort(msa_files, file_count, sizeof(char*), compare_natural);
	j=0;
    // Now loop over the sorted list and process each file
    for (int k = 0; k < file_count; k++) {
        strcpy(pf->msa_files[j], msa_files[k]);

        for (i = 0; i < 300; i++) {
            name[i] = '\0';
        }

        for (i = 0; i < (int)strlen(msa_files[k]) - 10; i++) {
            name[i] = msa_files[k][i];
        }

        name[i] = '\0';
        strcpy(name2, name);
        strncat(name2, "_taxonomy.txt", 250);
        strcpy(pf->tax_files[j], name2);

        strcpy(name3, "RAxML_bestTree.");
        strncat(name3, name, 250);
        strncat(name3, ".reroot", 250);
        strcpy(pf->tree_files[j], name3);

        j++;
        free(msa_files[k]);
    }
	free(msa_files);
    regfree(&regex1);
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
