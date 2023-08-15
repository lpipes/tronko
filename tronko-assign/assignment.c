#include "assignment.h"

/*type_of_PP **assignScores_Arr(int rootNum, int node, char *locQuery, int *positions, type_of_PP **scores, int alength){

	int child0 = treeArr[rootNum][node].up[0];
	int child1 = treeArr[rootNum][node].up[1];
	if (child0 == -1 && child1 == -1){
		scores[rootNum][node] = getscore_Arr(alength,node,rootNum,locQuery,positions);
		//treeArr[rootNum][node].score = getscore_Arr(alength,node,rootNum,locQuery,positions);
		//printf("node %d, score %lf, root %d\n",node,tree[node].score[rootNum],nodesToCutMinVar[rootNum]);
		return scores;
	}else if (child0 != -1 && child1 != -1){
		scores[rootNum][node] = getscore_Arr(alength,node,rootNum,locQuery,positions);
		//treeArr[rootNum][node].score = getscore_Arr(alength,node,rootNum);
		assignScores_Arr(rootNum, child0, locQuery, positions, scores, alength);
		assignScores_Arr(rootNum, child1, locQuery, positions, scores, alength);
	}*/
	/*if (child1 != -1 ){
		scores[rootNum][node] = getscore_Arr(alength,node,rootNum,locQuery,positions);
		//treeArr[rootNum][node].score = getscore_Arr(alength,node,rootNum);
		assignScores_Arr(rootNum, child1, locQuery, positions, scores, alength);
	}*/
/*}*/
void assignScores_Arr_paired( int rootNum, int node, char *locQuery, int *positions, type_of_PP ***scores, int alength, int search_number, int print_all_nodes, FILE* site_scores_file){
	//if (positions[0]==-1){
	//	return;
	//}
	int child0 = treeArr[rootNum][node].up[0];
	int child1 = treeArr[rootNum][node].up[1];
	if (child0 == -1 && child1 == -1){
		scores[search_number][rootNum][node] += getscore_Arr(alength,node,rootNum,locQuery,positions,print_all_nodes,site_scores_file);
		//scores[search_number][rootNum].nodeNumber = node;
		//if ( pair == 1 ){
		//	scores[search_number][rootNum].score1 = getscore_Arr(alength,node,rootNum,locQuery,positions);
		//}else{
		//	scores[search_number][rootNum].score2 = getscore_Arr(alength,node,rootNum,locQuery,positions);
		//}
		//treeArr[rootNum][node].score += getscore_Arr(alength,node,rootNum,locQuery,positions);
		//printf("node %d, score %lf, tree %d\n",node,scores[search_number][rootNum][node],rootNum);
	}else if(child0 != -1 && child1 != -1){
		scores[search_number][rootNum][node] += getscore_Arr(alength,node,rootNum,locQuery,positions,print_all_nodes,site_scores_file);
		//treeArr[rootNum][node].score += getscore_Arr(alength,node,rootNum,locQuery,positions);
		//printf("node %d, score %lf, tree %d\n",node,scores[search_number][rootNum][node],rootNum);
		//scores[search_number][rootNum].nodeNumber = node;
		//if ( pair == 1 ){
		//	scores[search_number][rootNum].score1 = getscore_Arr(alength,node,rootNum,locQuery,positions);
		//}else{
		//	scores[search_number][rootNum].score2 = getscore_Arr(alength,node,rootNum,locQuery,positions);
		//}
		assignScores_Arr_paired(rootNum, child0,locQuery, positions, scores, alength, search_number,print_all_nodes,site_scores_file);
		assignScores_Arr_paired(rootNum, child1, locQuery, positions,scores,alength, search_number,print_all_nodes,site_scores_file);
	}
	/*if(child1 != -1 ){
		scores[search_number][rootNum][node] +=  getscore_Arr(alength,node,rootNum,locQuery,positions);
		//treeArr[rootNum][node].score += getscore_Arr(alength,node,rootNum,locQuery,positions);
		//printf("node %d, score %lf, tree %d\n",node,scores[search_number][rootNum][node],rootNum);
		//scores[search_number][rootNum].nodeNumber = node;
		//if ( pair == 1 ){
		//	scores[search_number][rootNum].score1 = getscore_Arr(alength,node,rootNum,locQuery,positions);
		//}else{
		//	scores[search_number][rootNum].score2 = getscore_Arr(alength,node,rootNum,locQuery,positions);
		//}
		assignScores_Arr_paired(rootNum, child1, locQuery, positions,scores,alength, search_number);
	}*/
}
/*type_of_PP getscore(int alength, int node, int rootNum){
	type_of_PP score;
	int pos, i;
	score=0;
	if (when==2){
	  char* level;
	  level=malloc(sizeof(char)*20);
	  if (tree[node].taxIndex[1]==0){
	  level = "species";
	  }
	  if (tree[node].taxIndex[1]==1){
	  level = "genus";
	  }
	  if (tree[node].taxIndex[1]==2){
	  level="family";
	  }
	  if (tree[node].taxIndex[1]==3){
	  level="order";
	  }
	  if (tree[node].taxIndex[1]==4){
	  level="class";
	  }
	  if (tree[node].taxIndex[1]==5){
	  level="phylum";
	  }
	  printf("LCA Node %i, taxonomy (%s): %s, alength: %i %c\n",node,level,taxonomy[tree[node].taxIndex[0]][tree[node].taxIndex[1]],alength,locQuery[alength]);
	  }
	  if (alength < 40 ){
		  return 999999999999999;
	  }
	  int j=0;
	for (i=0; i<numbase; i++){//We assume that we have already ensured that the sequence only contains a, c, t, and g.  Missing data is not aligned
		pos=positions[j];//test if this is faster rather than referencing multiple times.  Compiler optimization may make the latter faster.
		if (i==pos){
			if (locQuery[j]=='-' && tree[node].taxIndex[1]==0 && node==rootNum){
				score=score+1;
			}else{
				//might be better simply to delete the parts of the query that does not align. This probably slows doww a lot.
				if (locQuery[j]=='a') score += PP[node][pos][0]; //test if a swtich statement is faster in this case. Might be compiler dependent.
				else if (locQuery[j]=='c') score += PP[node][pos][1];
				else if (locQuery[j]=='g') score += PP[node][pos][2];
				else if (locQuery[j]=='t') score += PP[node][pos][3];
				else score = score+1;
				//oldscore=score;
			}
		//}
			j++;
		}else{
			score = score+1;
		}
	}
	return score;
}*/
int checkPolyA(int rootNum, int node, int position){
	int i,j;
	int isPolyA = 0;
	for(i=0; i<position; i++){
		if ( treeArr[rootNum][node].posteriornc[i][0] == 1 ){
			isPolyA = 1;
		}else{
			isPolyA = 0;
			break;
		}
	}
	if ( isPolyA == 1 ){
		return 1;
	}
	for(i=numbaseArr[rootNum]-1; i>=position; i--){
		if ( treeArr[rootNum][node].posteriornc[i][0] == 1 ){
			isPolyA = 1;
		}else{
			isPolyA = 0;
			break;
		}
	}
	return isPolyA;
}
type_of_PP getscore_Arr(int alength, int node, int rootNum, char *locQuery, int *positions, int print_all_nodes, FILE* site_scores_file){
	type_of_PP score;
	int pos, i;
	i=0;
	score=0;
	type_of_PP overhang=0;
	int isPolyA=0;
	if (positions[i]==-1){
		score=9999999999;
		return score;
	}
	for (i=0; i<alength; i++){//We assume that we have already ensured that the sequence only contains a, c, t, and g.  Missing data is not aligned
		//pos=positions[i];//test if this is faster rather than referencing multiple times.  Compiler optimization may make the latter faster.
			if ( print_all_nodes == 1){
				fprintf(site_scores_file,"%d\t",positions[i]);
			}
			if (positions[i]==-1){
				score=score+log(0.01);
				if ( print_all_nodes == 1 ){
					fprintf(site_scores_file,"%lf\n",log(0.01));
				}
			}else{
				if ( treeArr[rootNum][node].posteriornc[positions[i]][0]==1 && locQuery[i]=='-' ){
					score=score+0; 
					if ( print_all_nodes == 1){
						fprintf(site_scores_file,"%lf\n",0);
					}
				}else{
					if( treeArr[rootNum][node].posteriornc[positions[i]][0] == 1){
						score = score + log(0.01);
						if ( print_all_nodes == 1){
							fprintf(site_scores_file,"%lf\n",log(0.01));
						}
					}else{
						if (locQuery[i]=='a' || locQuery[i]=='A'){
							score += treeArr[rootNum][node].posteriornc[positions[i]][0];
							if ( print_all_nodes == 1){
								fprintf(site_scores_file,"%lf\n",treeArr[rootNum][node].posteriornc[positions[i]][0]);
							}
						}else if (locQuery[i]=='c' || locQuery[i]=='C'){
							score += treeArr[rootNum][node].posteriornc[positions[i]][1];
							if ( print_all_nodes == 1){
								fprintf(site_scores_file,"%lf\n",treeArr[rootNum][node].posteriornc[positions[i]][1]);
							}
						}else if (locQuery[i]=='g' || locQuery[i]=='G'){
							score += treeArr[rootNum][node].posteriornc[positions[i]][2];
							if ( print_all_nodes == 1){
								fprintf(site_scores_file,"%lf\n",treeArr[rootNum][node].posteriornc[positions[i]][2]);
							}
						}else if (locQuery[i]=='t' || locQuery[i]=='T'){
							score += treeArr[rootNum][node].posteriornc[positions[i]][3];
							if ( print_all_nodes == 1){
								fprintf(site_scores_file,"%lf\n",treeArr[rootNum][node].posteriornc[positions[i]][3]);
							}
						}else if (locQuery[i]=='-'){
							score += log(0.25);
							if ( print_all_nodes == 1){
								fprintf(site_scores_file,"%lf\n",log(0.25));
							}
						}
					}
				}
			}
			//printf("Node: %d Position: %d Score: %lf\n",node,pos,score);
	}
	
	return score;
}
