#include "printtree.h"

/*old code for printing a tree*/
void printtreeArr(int whichRoot){
  int i,j;
  printf("\n");
  j=whichRoot;
  for (i=0; i<2*numspecArr[j]-1; i++){
    if (treeArr[j][i].up[0] != -1){
      printf("Node %i: up: (%i, %i) down: %i, name: %s, taxIndex[0]: %d, taxIndex[1]: %d, depth: %d, nd: %d, ",i,treeArr[j][i].up[0],treeArr[j][i].up[1],treeArr[j][i].down,treeArr[j][i].name,treeArr[j][i].taxIndex[0],treeArr[j][i].taxIndex[1],treeArr[j][i].depth,treeArr[j][i].nd);
	if ( treeArr[j][i].taxIndex[0] != -1 && treeArr[j][i].taxIndex[1] != -1 ){
		printf("taxonomy, %s",taxonomyArr[j][treeArr[j][i].taxIndex[0]][treeArr[j][i].taxIndex[1]]);
	}else{
		printf("taxonomy, Domain");
	}
    }else{
	printf("(Node %i): up: (%i, %i) down: %i, name: %s, taxIndex[0]: %d, taxIndex[1]: %d, depth: %d, nd: %d, ",i,treeArr[j][i].up[0],treeArr[j][i].up[1],treeArr[j][i].down,treeArr[j][i].name,treeArr[j][i].taxIndex[0],treeArr[j][i].taxIndex[1],treeArr[j][i].depth,treeArr[j][i].nd);
	if ( treeArr[j][i].taxIndex[0] != -1 && treeArr[j][i].taxIndex[1] != -1 ){
	       	printf("taxonomy, %s",taxonomyArr[j][treeArr[j][i].taxIndex[0]][treeArr[j][i].taxIndex[1]]);
	}else{
		printf("taxonomy, Domain");
	}
	}
    if (i != rootArr[j] ){
      printf(" (bl: %f)\n",treeArr[j][i].bl);
	}else {
		printf("\n");
	}
  }
}
void printTreeFile(int numberOfTrees, int max_nodename, int max_tax_name, int max_lineTaxonomy, Options opt){
	int i, j, k, l;
	char buf[BUFFER_SIZE];
	struct stat st = {0};
	if ( stat(opt.partitions_directory, &st) == -1){
		mkdir(opt.partitions_directory, 0700);
	}
	snprintf(buf,BUFFER_SIZE,"%s/reference_tree.txt",opt.partitions_directory);
	FILE *outputTree = fopen(buf,"w");
	if  ( outputTree == NULL ){ printf("Error opening reference tree file!\n"); exit(1); }
	fprintf(outputTree,"%d\n",numberOfTrees);
	fprintf(outputTree,"%d\n",max_nodename);
	fprintf(outputTree,"%d\n",max_tax_name);
	fprintf(outputTree,"%d\n",max_lineTaxonomy);
	for (i=0; i<numberOfTrees;i++){
		fprintf(outputTree,"%d\t%d\t%d\n",numbaseArr[i],rootArr[i],numspecArr[i]);
	}
	for(i=0; i<numberOfTrees; i++){
		for(j=0; j<numspecArr[i]; j++){
			for(k=0; k<7; k++){
				if (k==6){
					fprintf(outputTree,"%s\n",taxonomyArr[i][j][k]);
				}else{
					fprintf(outputTree,"%s;",taxonomyArr[i][j][k]);
				}
			}
		}
	}
	for(i=0; i<numberOfTrees;i++){
		for(j=0;j<2*numspecArr[i]-1;j++){
			fprintf(outputTree,"%d\t",i);
			fprintf(outputTree,"%d\t",j);
			fprintf(outputTree,"%d\t",treeArr[i][j].up[0]);
			fprintf(outputTree,"%d\t",treeArr[i][j].up[1]);
			fprintf(outputTree,"%d\t",treeArr[i][j].down);
			fprintf(outputTree,"%d\t",treeArr[i][j].depth);
			fprintf(outputTree,"%d\t",treeArr[i][j].taxIndex[0]);
			fprintf(outputTree,"%d\t",treeArr[i][j].taxIndex[1]);
			if ( treeArr[i][j].up[0] != -1 && treeArr[i][j].up[1] != -1 ){
				fprintf(outputTree,"\n");
			}else{
				fprintf(outputTree,"%s\n",treeArr[i][j].name);
			}
			for (k=0; k<numbaseArr[i]-1; k++){
				for (l=0; l<3; l++){
					fprintf(outputTree,"%.17g\t",treeArr[i][j].posteriornc[k][l]);
				}
				fprintf(outputTree,"%.17g\n",treeArr[i][j].posteriornc[k][3]);
			}
			fprintf(outputTree,"%.17g\t%.17g\t%.17g\t%.17g\n",treeArr[i][j].posteriornc[numbaseArr[i]-1][0],treeArr[i][j].posteriornc[numbaseArr[i]-1][1],treeArr[i][j].posteriornc[numbaseArr[i]-1][2],treeArr[i][j].posteriornc[numbaseArr[i]-1][3]);
		}
	}
	fclose(outputTree);
}
void printTaxonomyArrToFile(int numberOfTrees){
	int i,j,k;
	FILE *outputTaxArr = fopen("/space/s2/lenore/partitions2/taxonomy_reference.txt","w");
	if ( outputTaxArr == NULL ){ printf("Error opening file!\n"); exit(1); }
	for(i=0; i<numberOfTrees;i++){
		fprintf(outputTaxArr,"TREE %d\n",i);
		for(j=0;j<numspecArr[i];j++){
			fprintf(outputTaxArr,"%d",j);
			for(k=0;k<7;k++){
				fprintf(outputTaxArr,"\t%s",taxonomyArr[i][j][k]);
			}
			fprintf(outputTaxArr,"\n");
		}
	}
	fclose(outputTaxArr);
}

void printPP(){
	int i,j;
	FILE *outputPP;
	outputPP=fopen("/space/s2/lenore/partitions2/initial_tree_PP.txt","w");
	if (outputPP==NULL){ printf("Error opening file!\n"); exit(1); }
	printf("numspec is %d\n",numspec);
	for(i=0; i<2*numspec-1; i++){
		printf("NODE\t%d\t%d\n",i,numbase);
		//fprintf(stderr,"NODE\t%d\n",i);
		for(j=0;j<numbase;j++){
			printf("%d\t%lf\t%lf\t%lf\t%lf\n",j,PP[i][j][0],PP[i][j][1],PP[i][j][2],PP[i][j][3]);
			fprintf(outputPP,"%d\t%lf\t%lf\t%lf\t%lf\n",j,PP[i][j][0],PP[i][j][1],PP[i][j][2],PP[i][j][3]);
		}
	}
	fclose(outputPP);
}
void printPP_Arr(int numberOfTrees){
	int i,j,k;
	FILE *outputPP;
	outputPP=fopen("/space/s2/lenore/partitions2/PP_Arr.txt","w");
	if (outputPP==NULL){ printf("Error opening file!\n"); exit(1); }
	for(i=0; i<numberOfTrees; i++){
		fprintf(outputPP,"TREE\t%d\n",i);
		for(j=0; j<2*numspecArr[i]-1; j++){
			fprintf(outputPP,"NODE\t%d\n",j);
			for(k=0; k<numbaseArr[i]; k++){
				fprintf(outputPP,"%d\t%lf\t%lf\t%lf\t%lf\n",k,PP_Arr[i][j][k][0],PP_Arr[i][j][k][1],PP_Arr[i][j][k][2],PP_Arr[i][j][k][3]);
			}
		}
	}
	fclose(outputPP);
}
