#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include "bwa.h"
#include "bwamem.h"
#include "kvec.h"
#include "utils.h"
#include "bntseq.h"
#include "../hashmap.h"
#include "../hashmap_base.h"
#include "kseq.h"
#include "../global.h"
KSEQ_DECLARE(gzFile)

extern unsigned char nst_nt4_table[256];

void *kopen(const char *fn, int *_fd);
int kclose(void *a);
void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

typedef struct {
	kseq_t *ks, *ks2;
	mem_opt_t *opt;
	mem_pestat_t *pes0;
	int64_t n_processed;
	int copy_comment, actual_chunk_size;
	bwaidx_t *idx;
	/*char **names;
	char **read1;
	char **read2;*/
	int number_of_seqs;
	bwaMatches *results;
	int concordant;
	int ntree;
	int startline;
	int paired;
	int start;
	int end;
	int max_query_length;
	int max_readname_length;
	int max_acc_name;
} ktp_aux_t;

typedef struct {
	ktp_aux_t *aux;
	int n_seqs;
	bseq1_t *seqs;
} ktp_data_t;

static void *process(void *shared, int step, void *_data)
{
	ktp_aux_t *aux = (ktp_aux_t*)shared;
	ktp_data_t *data = (ktp_data_t*)_data;
	int i;
	if (step == 0) {
		ktp_data_t *ret;
		int64_t size = 0;
		ret = calloc(1, sizeof(ktp_data_t));
		ret->seqs = bseq_read(aux->actual_chunk_size, &ret->n_seqs, aux->ks, aux->ks2, aux->number_of_seqs, aux->paired, aux->start, aux->end, aux->max_readname_length);
		if (ret->seqs == 0) {
			free(ret);
			return 0;
		}
		//if (!aux->copy_comment)
		//	for (i = 0; i < ret->n_seqs; ++i) {
		//		free(ret->seqs[i].comment);
		//		ret->seqs[i].comment = 0;
		//	}
		for (i = 0; i < ret->n_seqs; ++i) size += ret->seqs[i].l_seq;
		if (bwa_verbose >= 3)
			fprintf(stderr, "[M::%s] read %d sequences (%ld bp)...\n", __func__, ret->n_seqs, (long)size);
		return ret;
	} else if (step == 1) {
		const mem_opt_t *opt = aux->opt;
		const bwaidx_t *idx = aux->idx;
		if (opt->flag & MEM_F_SMARTPE) {
			bseq1_t *sep[2];
			int n_sep[2];
			mem_opt_t tmp_opt = *opt;
			bseq_classify(data->n_seqs, data->seqs, n_sep, sep);
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] %d single-end sequences; %d paired-end sequences\n", __func__, n_sep[0], n_sep[1]);
			if (n_sep[0]) {
				tmp_opt.flag &= ~MEM_F_PE;
				mem_process_seqs(&tmp_opt, idx->bwt, idx->bns, idx->pac, aux->n_processed, n_sep[0], sep[0], 0, aux->concordant, aux->startline);
				for (i = 0; i < n_sep[0]; ++i)
					data->seqs[sep[0][i].id].sam = sep[0][i].sam;
			}
			if (n_sep[1]) {
				tmp_opt.flag |= MEM_F_PE;
				mem_process_seqs(&tmp_opt, idx->bwt, idx->bns, idx->pac, aux->n_processed + n_sep[0], n_sep[1], sep[1], aux->pes0, aux->concordant, aux->startline);
				for (i = 0; i < n_sep[1]; ++i)
					data->seqs[sep[1][i].id].sam = sep[1][i].sam;
			}
			free(sep[0]); free(sep[1]);
		} else mem_process_seqs(opt, idx->bwt, idx->bns, idx->pac, aux->n_processed, data->n_seqs, data->seqs, aux->pes0, aux->concordant,aux->startline);
		aux->n_processed += data->n_seqs;
		return data;
	} else if (step == 2) {
		int j=0;
		int l=0;
		int k=0;
		int no_add=0;
		int success=1;
	HASHMAP(char, leafMap) map;
	hashmap_init(&map, hashmap_hash_string, strcmp);
	for(i=0; i<aux->ntree; i++){
		for(j=numspecArr[i]-1; j<2*numspecArr[i]-1; j++){
			struct leafMap *l;
			l = malloc(sizeof(*l));
			l->name = treeArr[i][j].name;
			l->root = i;
			l->node = j;
			hashmap_put(&map,l->name,l);
		}
	}
	j=0;
		for (i = 0; i < data->n_seqs; ++i) {
			/*if (data->seqs[i].sam && aux->concordant==1){
				err_fputs(data->seqs[i].sam, stdout);
				char* token;
				char readname[MAXREADNAME];
				char read1[MAX_NODENAME];
				char read2[MAX_NODENAME];
				token=strtok(data->seqs[i].sam,"\n");
				while(success !=0 && token != NULL){
					success = sscanf(token, "%s %*d %s %*d %*d %*s %s %*d %*d %*s %*s %*s %*s %*s %*s",readname,read1,read2);
					if ( j > 0 && strcmp(readname,aux->results[j-1].readname)==0 ){
						for ( k=0; k<aux->ntree; k++){
							if ( aux->results[j-1].matches[k]==-1 ){
								break;
							}
						}
						for (l=0; l<k; l++){
							if (aux->results[j-1].matches[l]==hashmap_get(&map,read1)){
								no_add=1;
							}
						}
						if ( no_add==0 ){
							aux->results[j-1].matches[k] = hashmap_get(&map,read1);
							k++;
						}
						no_add=0;
						for(l=0; l<k; l++){
							if ( aux->results[j-1].matches[l] == hashmap_get(&map,read2) && strcmp(read2,"=") != 0 && strcmp(read2,"*") != 0){
								no_add=1;
							}
						}
						if (no_add==0){
							aux->results[j-1].matches[k] = hashmap_get(&map,read2);
						}
						no_add=0;
				}else{
					aux->results[j].readname=token;
					token=strtok(NULL,"\t");
					token=strtok(NULL,"\t");
					aux->results[j].matches[0] = token;
					j++;
				}*/
			//}else if ( data->seqs[i].sam && aux->concordant==0){
			if ( data->seqs[i].sam ){
				char *token;
				char readname[aux->max_readname_length];
				char read1[aux->max_acc_name];
				char read2[aux->max_acc_name];
				char cigar[MAX_CIGAR];
				int start_position;
				char *rest = data->seqs[i].sam;
				int decimal=0;
				int decimal_general=0;
				while(token=strtok_r(rest,"\n",&rest)){
					success = sscanf(token, "%s %d %s %d %*d %s %s %*d %*d %*s %*s %*s %*s %*s %*s %*[^.,;]",readname,&decimal,read1,&start_position,cigar,read2);
					/*decimal_general = dec2bin(decimal);
					if (aux->paired==1){
						decimal = dec2bin(decimal);
					}*/
					//if decimal==1, first in pair if decimal==0, second in pair
				//if ( j > 0 && strcmp(readname,aux->results[j-1].readname)==0 ){
				int whichMat=0;
				if (aux->paired==1 && j>0){
					if (strcmp(readname,pairedQueryMat->forward_name[aux->startline+j-1])==0){
						whichMat=1;
					}
				}else if (j>0){
					if (strcmp(readname,singleQueryMat->name[aux->startline+j-1])==0){
						whichMat=1;
					}
				}
				if ( j > 0 && whichMat==1 ){
					for(k=0; k<MAX_NUM_BWA_MATCHES; k++){
						//if ( aux->results[j-1].concordant_matches[k] == -1 ){
						//	break;
						//}
						//if ( strlen(aux->results[j-1].concordant_leaf_matches[k])==0 ){
						//	break;
						//}
						if (aux->results[j-1].concordant_matches_roots[k] == -1){
							break;
						}
					}
					for(l=0; l<k; l++){
						if ( strcmp(read2,"=")==0 ){
							//if ( aux->results[j-1].concordant_matches[l]==hashmap_get(&map,read1) ){
							//	no_add=1;
							//}
							//if ( strcmp(aux->results[j-1].concordant_leaf_matches[l],read1)==0 ){
							//	no_add=1;
							//}
							struct leafMap *leaf_map;
							leaf_map=hashmap_get(&map,read1);
							if ( aux->results[j-1].concordant_matches_roots[l] == leaf_map->root && aux->results[j-1].concordant_matches_nodes[l] == leaf_map->node){
								no_add=1;
							}
						}
					}
					if (aux->concordant==1 && no_add==0 && strcmp(read2,"=")!=0){ no_add=1;}
					if (no_add==0 && strcmp(read2,"=")==0){
							struct leafMap *leaf_map;
						leaf_map=hashmap_get(&map,read1);
						aux->results[j-1].concordant_matches_roots[k] = leaf_map->root;
						aux->results[j-1].concordant_matches_nodes[k] = leaf_map->node;
						//aux->results[j-1].concordant_matches[k] = hashmap_get(&map,read1);
						//strcpy(aux->results[j-1].concordant_leaf_matches[k],read1);
						if ( aux->results[j-1].use_portion == 1){
							strcpy(aux->results[j-1].cigars_forward[k],cigar);
							aux->results[j-1].starts_forward[k] = start_position;
						}
						k++;
					}
					no_add=0;
					for(l=0; l<k; l++){
						//if ( aux->results[j-1].concordant_matches[l] == hashmap_get(&map,read2) || strcmp(read2,"=") == 0 || strcmp(read2,"*") == 0){
						//	no_add=1;
						//}
						//if ( strcmp(aux->results[j-1].concordant_leaf_matches[l],read2)==0 || strcmp(read2,"=") == 0 || strcmp(read2,"*") == 0){
						//	no_add=1;
						//}
						if ( strcmp(read2,"=")== 0 || strcmp(read2,"*")==0){
							no_add=1;
							break;
						}
						struct leafMap *leaf_map;
						leaf_map=hashmap_get(&map,read2);
						if ( aux->results[j-1].concordant_matches_roots[l] == leaf_map->root && aux->results[j-1].concordant_matches_nodes[l] == leaf_map->node ){
							no_add=1;
						}
					}
					//if (aux->concordant==1){ no_add=1;}
					if (k==0){ no_add=1; }
					if (no_add==0 && strcmp(read2,"*")!=0 && strcmp(read1,"=")!=0){
						struct leafMap *leaf_map;
						leaf_map=hashmap_get(&map,read2);
						aux->results[j-1].concordant_matches_roots[k] = leaf_map->root;
						aux->results[j-1].concordant_matches_nodes[k] = leaf_map->node;
						//aux->results[j-1].concordant_matches[k] = hashmap_get(&map,read2);
						//strcpy(aux->results[j-1].concordant_leaf_matches[k],read2);
						if (aux->results[j-1].use_portion==1){
							strcpy(aux->results[j-1].cigars_reverse[k-1],cigar);
							aux->results[j-1].starts_reverse[k-1] = start_position;
						}
					}
					no_add=0;
					for(k=0; k<MAX_NUM_BWA_MATCHES; k++){
						//if ( aux->results[j-1].discordant_matches[k]== -1 ){
						//	break;
						//}
						//if ( strlen(aux->results[j-1].discordant_leaf_matches[k])==0 ){
						//	break;
						//}
						if ( aux->results[j-1].discordant_matches_roots[k]==-1){
							break;
						}
					}
					for(l=0; l<k; l++){
						if (strcmp(read2,"=") != 0 ){
							//if (aux->results[j-1].discordant_matches[l]==hashmap_get(&map,read1)){
							//	no_add=1;
							//}
							//if (strcmp(aux->results[j-1].discordant_leaf_matches[l],read1)==0){
							//	no_add=1;
							//}
							struct leafMap *leaf_map;
							leaf_map=hashmap_get(&map,read1);
							if ( aux->results[j-1].discordant_matches_roots[l] == leaf_map->root && aux->results[j-1].discordant_matches_nodes[l] == leaf_map->node ){
								no_add=1;
							}
						}
					}
					if (no_add==0 && strcmp(read2,"=") != 0){
							struct leafMap *leaf_map;
						leaf_map=hashmap_get(&map,read1);
						aux->results[j-1].discordant_matches_roots[k] = leaf_map->root;
						aux->results[j-1].discordant_matches_nodes[k] = leaf_map->node;
						//aux->results[j-1].discordant_matches[k]=hashmap_get(&map,read1);
						//strcpy(aux->results[j-1].discordant_leaf_matches[k],read1);
						if (aux->results[j-1].use_portion==1 && decimal == 1){
							strcpy(aux->results[j-1].cigars_forward[k],cigar);
							aux->results[j-1].starts_forward[k] = start_position;
						}
						k++;
					}
					no_add=0;
					for(l=0; l<k; l++){
						//if (aux->results[j-1].discordant_matches[l] == hashmap_get(&map,read2) || strcmp(read2,"=") ==0 || strcmp(read2,"*") == 0){
						//	no_add=1;
						//}
						if (strcmp(read2,"*")==0 || strcmp(read2,"=")==0){
							no_add=1;
							break;
						}
						struct leafMap *leaf_map;
						leaf_map=hashmap_get(&map,read2);
						if ( aux->results[j-1].discordant_matches_roots[l] == leaf_map->root && aux->results[j-1].discordant_matches_nodes[l] == leaf_map->node){
							no_add=1;
						}
						//if ( strcmp(aux->results[j-1].discordant_leaf_matches[l],read2)==0 || strcmp(read2,"=")==0 || strcmp(read2,"*") == 0 ){
						//	no_add=1;
						//}
					}
					if (no_add==0 && strcmp(read2,"=")!=0 && strcmp(read2,"*")!=0){
						struct leafMap *leaf_map;
						leaf_map=hashmap_get(&map,read2);
						aux->results[j-1].discordant_matches_roots[k] = leaf_map->root;
						aux->results[j-1].discordant_matches_nodes[k] = leaf_map->node;
						//aux->results[j-1].discordant_matches[k] = hashmap_get(&map,read2);
						//strcpy(aux->results[j-1].discordant_leaf_matches[k],read2);
						if (aux->results[j-1].use_portion==1){
							strcpy(aux->results[j-1].cigars_reverse[k],cigar);
							aux->results[j-1].starts_reverse[k] = start_position;
						}
					}
					no_add=0;
					if ( decimal == 0 && aux->results[j-1].use_portion==1){
						strcpy(aux->results[j-1].cigars_reverse[k],cigar);
						aux->results[j-1].starts_reverse[k] = start_position;
					}
				}else if (j==0 && strcmp(read2,"=")==0){
					//strcpy(aux->results[j].readname,readname);
							struct leafMap *leaf_map;
					leaf_map=hashmap_get(&map,read1);
					aux->results[j].concordant_matches_roots[0] = leaf_map->root;
					aux->results[j].concordant_matches_nodes[0] = leaf_map->node;
					//aux->results[j].concordant_matches[0] = hashmap_get(&map,read1);
					//strcpy(aux->results[j].concordant_leaf_matches[0],read1);
					if (aux->results[j].use_portion==1){
						strcpy(aux->results[j].cigars_forward[0],cigar);
						aux->results[j].starts_forward[0] = start_position;
					}
					//if (strcmp(read2,"=") != 0 && strcmp(read2,"*") != 0 && aux->results[j].concordant_matches[0] != hashmap_get(&map,read2)){
					//if (strcmp(read2,"=") != 0 && strcmp(read2,"*") != 0 && strcmp(aux->results[j].concordant_leaf_matches[0],read2) != 0){
					if (strcmp(read2,"=") != 0 && strcmp(read2,"*") != 0 && (aux->results[j].concordant_matches_roots[0] != leaf_map->root && aux->results[j].concordant_matches_nodes[0] != leaf_map->node)){
							struct leafMap *leaf_map;
						leaf_map=hashmap_get(&map,read2);
						aux->results[j].concordant_matches_roots[1] = leaf_map->root;
						aux->results[j].concordant_matches_nodes[1] = leaf_map->node;
						//aux->results[j].concordant_matches[1] = hashmap_get(&map,read2);
						//strcpy(aux->results[j].concordant_leaf_matches[1],read2);
						if (aux->results[j].use_portion==1){
							strcpy(aux->results[j].cigars_reverse[0],cigar);
							aux->results[j].starts_reverse[0] = start_position;
						}
					}
					j++;
				}else if (j==0 && strcmp(read2,"=")!=0 && aux->paired != 0 && strcmp(read1,"*")!=0){
					//strcpy(aux->results[j].readname,readname);
					struct leafMap *leaf_map;
					leaf_map=hashmap_get(&map,read1);
					aux->results[j].discordant_matches_roots[0] = leaf_map->root;
					aux->results[j].discordant_matches_nodes[0] = leaf_map->node;
					//aux->results[j].discordant_matches[0] = hashmap_get(&map,read1);
					//strcpy(aux->results[j].discordant_leaf_matches[0],read1); 
					if (aux->results[j].use_portion==1){
						strcpy(aux->results[j].cigars_forward[0],cigar);
						aux->results[j].starts_forward[0] = start_position;
					}
					//if (strcmp(read2,"*")!=0 && aux->results[j].discordant_matches[0] != hashmap_get(&map,read2)){
					//if ( strcmp(read2,"*") != 0 && strcmp(aux->results[j].discordant_leaf_matches[0],read2) != 0 ){
					if ( strcmp(read2,"*") != 0 && (aux->results[j].discordant_matches_roots[0] != leaf_map->root && aux->results[j].discordant_matches_nodes[0] != leaf_map->node)){
							struct leafMap *leaf_map;
						leaf_map=hashmap_get(&map,read2);
						aux->results[j].discordant_matches_roots[1] = leaf_map->root;
						aux->results[j].discordant_matches_nodes[1] = leaf_map->node;
						//aux->results[j].discordant_matches[1] = hashmap_get(&map,read2);
						//strcpy(aux->results[j].discordant_leaf_matches[1],read2);
						if (aux->results[j].use_portion==1){
							strcpy(aux->results[j].cigars_reverse[0],cigar);
							aux->results[j].starts_reverse[0] = start_position;
						}
					}
					j++;
				}else{
					if (aux->paired != 0 && strcmp(pairedQueryMat->forward_name[aux->startline+j],readname)==0 && strcmp(read2,"=")==0){
						if ( decimal == 1 || decimal==2){
							//strcpy(aux->results[j].readname,readname);
							struct leafMap *leaf_map;
							leaf_map=hashmap_get(&map,read1);
							aux->results[j].concordant_matches_roots[0] = leaf_map->root;
							aux->results[j].concordant_matches_nodes[0] = leaf_map->node;
							//aux->results[j].concordant_matches[0] = hashmap_get(&map,read1);
							//strcpy(aux->results[j].concordant_leaf_matches[0],read1);
							if (aux->results[j].use_portion==1){
								strcpy(aux->results[j].cigars_forward[0],cigar);
								aux->results[j].starts_forward[0] = start_position;
							}
						}
						//if (aux->results[j].use_portion==1){
						//	strcpy(aux->results[j].cigars_reverse[0],cigar);
						//	aux->results[j].starts_reverse[0] = start_position;
						//}
						//if (strcmp(read2,"=") != 0 && strcmp(read2,"*") != 0 && aux->results[j].concordant_matches[0] != hashmap_get(&map,read2)){
						//if (strcmp(read2,"=") != 0 && strcmp(read2,"*") != 0 && strcmp(aux->results[j].concordant_leaf_matches[0],read2) != 0 ){
						if (strcmp(read2,"=") != 0 && strcmp(read2,"*") != 0 && strcmp(read1,read2) != 0 ){
							struct leafMap *leaf_map;
							leaf_map=hashmap_get(&map,read2);
							aux->results[j].concordant_matches_roots[1] = leaf_map->root;
							aux->results[j].concordant_matches_nodes[1] = leaf_map->node;
							//aux->results[j].concordant_matches[1] = hashmap_get(&map,read2);
							//strcpy(aux->results[j].concordant_leaf_matches[1],read2);
							if (aux->results[j].use_portion==1){
								strcpy(aux->results[j].cigars_reverse[0],cigar);
								aux->results[j].starts_reverse[0] = start_position;
							}
						}
						j++;
					}else if (aux-> paired != 0 && strcmp(pairedQueryMat->forward_name[aux->startline+j],readname)==0 && strcmp(read2,"=")!=0 && strcmp(read1,"*")!=0){
						//strcpy(aux->results[j].readname,readname);
							struct leafMap *leaf_map;
						leaf_map=hashmap_get(&map,read1);
						aux->results[j].discordant_matches_roots[0] = leaf_map->root;
						aux->results[j].discordant_matches_nodes[0] = leaf_map->node;
						//aux->results[j].discordant_matches[0] = hashmap_get(&map,read1);
						//if (strcmp(read2,"*") != 0 && aux->results[j].discordant_matches[0] != hashmap_get(&map,read2) ){
						//if (strcmp(read2,"*") != 0 && strcmp(aux->results[j].discordant_leaf_matches[0],read2) != 0 ){
						if (strcmp(read2,"*") != 0 && strcmp(read1,read2)!=0 ){
							struct leafMap *leaf_map;
							leaf_map=hashmap_get(&map,read2);
							aux->results[j].discordant_matches_roots[1] = leaf_map->root;
							aux->results[j].discordant_matches_nodes[1] = leaf_map->node;
							//aux->results[j].discordant_matches[1] = hashmap_get(&map,read2);
							//strcpy(aux->results[j].discordant_leaf_matches[1],read2);
							if (aux->results[j].use_portion==1){
								strcpy(aux->results[j].cigars_reverse[0],cigar);
								aux->results[j].starts_reverse[0] = start_position;
							}
						}
						j++;
					}else if (aux->paired==0 && strcmp(singleQueryMat->name[aux->startline+j],readname)==0 && strcmp(read2,"=")==0){
						//strcpy(aux->results[j].readname,readname);
							struct leafMap *leaf_map;
						leaf_map=hashmap_get(&map,read1);
						aux->results[j].concordant_matches_roots[0] = leaf_map->root;
						aux->results[j].concordant_matches_nodes[0] = leaf_map->node;
						//aux->results[j].concordant_matches[0] = hashmap_get(&map,read1);
						//strcpy(aux->results[j].concordant_leaf_matches[0],read1);
						if (aux->results[j].use_portion==1){
							strcpy(aux->results[j].cigars_forward[0],cigar);
							aux->results[j].starts_forward[0] = start_position;
						}
						//if (strcmp(read2,"=") != 0 && strcmp(read2,"*") != 0 && strcmp(aux->results[j].concordant_leaf_matches[0],read2) != 0 ){
						if (strcmp(read2,"=") != 0 && strcmp(read2,"*") != 0 && strcmp(read1,read2)!=0 ){
							struct leafMap *leaf_map;
							leaf_map=hashmap_get(&map,read2);
							aux->results[j].concordant_matches_roots[1] = leaf_map->root;
							aux->results[j].concordant_matches_nodes[1] = leaf_map->node;
							//aux->results[j].concordant_matches[1] = hashmap_get(&map,read2);
							//strcpy(aux->results[j].concordant_leaf_matches[1],read2);
							if (aux->results[j].use_portion==1){
								strcpy(aux->results[j].cigars_reverse[0],cigar);
								aux->results[j].starts_reverse[0] = start_position;
							}
						}
						j++;
					}else if (aux->paired==0 && strcmp(singleQueryMat->name[aux->startline+j],readname)==0 && strcmp(read2,"=")!=0){
						//strcpy(aux->results[j].readname,readname);
						if ( strcmp(read1,"*") != 0 ){
							struct leafMap *leaf_map;
							leaf_map=hashmap_get(&map,read1);
							aux->results[j].discordant_matches_roots[0] = leaf_map->root;
							aux->results[j].discordant_matches_nodes[0] = leaf_map->node;
							//aux->results[j].discordant_matches[0] = hashmap_get(&map,read1);
							//strcpy(aux->results[j].discordant_leaf_matches[0],read1);
							if(aux->results[j].use_portion==1){
								strcpy(aux->results[j].cigars_forward[0],cigar);
								aux->results[j].starts_forward[0]=start_position;
							}
							if (strcmp(read2,"*") != 0 && strcmp(read1,read2) != 0 ){
							struct leafMap *leaf_map;
								leaf_map=hashmap_get(&map,read2);
								aux->results[j].discordant_matches_roots[1] = leaf_map->root;
								aux->results[j].discordant_matches_nodes[1] = leaf_map->node;
								//aux->results[j].discordant_matches[1] = hashmap_get(&map,read2);
								//strcpy(aux->results[j].discordant_leaf_matches[1],read2);
								if (aux->results[j].use_portion==1){
									strcpy(aux->results[j].cigars_reverse[0],cigar);
									aux->results[j].starts_reverse[0] = start_position;
								}
							}
						}else{
							aux->results[j].concordant_matches_roots[0]=-1;
							aux->results[j].concordant_matches_nodes[0]=-1;
							aux->results[j].discordant_matches_nodes[0]=-1;
							aux->results[j].discordant_matches_roots[0]=-1;
							//aux->results[j].concordant_matches[0]=-2;
							//aux->results[j].discordant_matches[0]=-2;
						}
						j++;	
					}else{
						if (aux->paired != 0){
							while(strcmp(pairedQueryMat->forward_name[aux->startline+j],readname)!=0){
								//strcpy(aux->results[j].readname,pairedQueryMat->forward_name[aux->startline+j]);
								aux->results[j].concordant_matches_roots[0]=-1;
								aux->results[j].concordant_matches_nodes[0]=-1;
								aux->results[j].discordant_matches_nodes[0]=-1;
								aux->results[j].discordant_matches_roots[0]=-1;
								//aux->results[j].concordant_matches[0]=-2;
								//aux->results[j].discordant_matches[0]=-2;
								j++;
							}
						}else{
							while(strcmp(singleQueryMat->name[aux->startline+j],readname)!=0){
								//strcpy(aux->results[j].readname,singleQueryMat->name[aux->startline+j]);
								aux->results[j].concordant_matches_roots[0]=-1;
								aux->results[j].concordant_matches_nodes[0]=-1;
								aux->results[j].discordant_matches_nodes[0]=-1;
								aux->results[j].discordant_matches_roots[0]=-1;
								//aux->results[j].concordant_matches[0]=-2;
								//aux->results[j].discordant_matches[0]=-2;
								j++;
							}
						}
					}
				}
				}
			free(data->seqs[i].name);
			free(data->seqs[i].seq);
			free(data->seqs[i].sam);
			}
		}
		free(data->seqs); free(data);
		struct leafMap *blob;
		hashmap_foreach_data(blob,&map){
			free(blob);
		}
		hashmap_cleanup(&map);
		return 0;
	}
	return 0;
}

static void update_a(mem_opt_t *opt, const mem_opt_t *opt0)
{
	if (opt0->a) { // matching score is changed
		if (!opt0->b) opt->b *= opt->a;
		if (!opt0->T) opt->T *= opt->a;
		if (!opt0->o_del) opt->o_del *= opt->a;
		if (!opt0->e_del) opt->e_del *= opt->a;
		if (!opt0->o_ins) opt->o_ins *= opt->a;
		if (!opt0->e_ins) opt->e_ins *= opt->a;
		if (!opt0->zdrop) opt->zdrop *= opt->a;
		if (!opt0->pen_clip5) opt->pen_clip5 *= opt->a;
		if (!opt0->pen_clip3) opt->pen_clip3 *= opt->a;
		if (!opt0->pen_unpaired) opt->pen_unpaired *= opt->a;
	}
}

int main_mem(char* databaseFile, int number_of_seqs, int number_of_threads, bwaMatches* bwa_results, int concordant, int numberOfTrees, int startline, int paired, int start, int end, int max_query_length, int max_readname_length, int max_acc_name)
{
	mem_opt_t *opt, opt0;
	int fd, fd2, i, c, ignore_alt = 0, no_mt_io = 0;
	int fixed_chunk_size = -1;
	gzFile fp, fp2 = 0;
	char *p, *rg_line = 0, *hdr_line = 0;
	const char *mode = 0;
	void *ko = 0, *ko2 = 0;
	mem_pestat_t pes[4];
	ktp_aux_t aux;

	memset(&aux, 0, sizeof(ktp_aux_t));
	memset(pes, 0, 4 * sizeof(mem_pestat_t));
	for (i = 0; i < 4; ++i) pes[i].failed = 1;

	aux.opt = opt = mem_opt_init();
	memset(&opt0, 0, sizeof(mem_opt_t));
	bwa_verbose=1;
	/*while ((c = getopt(argc, argv, "51qpaMCSPVYjuk:c:v:s:r:t:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:o:f:W:x:G:h:y:K:X:H:")) >= 0) {
		if (c == 'k') opt->min_seed_len = atoi(optarg), opt0.min_seed_len = 1;
		else if (c == '1') no_mt_io = 1;
		else if (c == 'x') mode = optarg;
		else if (c == 'w') opt->w = atoi(optarg), opt0.w = 1;
		else if (c == 'A') opt->a = atoi(optarg), opt0.a = 1;
		else if (c == 'B') opt->b = atoi(optarg), opt0.b = 1;
		else if (c == 'T') opt->T = atoi(optarg), opt0.T = 1;
		else if (c == 'U') opt->pen_unpaired = atoi(optarg), opt0.pen_unpaired = 1;
		else if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'P') opt->flag |= MEM_F_NOPAIRING;
		else if (c == 'a') opt->flag |= MEM_F_ALL;
		else if (c == 'p') opt->flag |= MEM_F_PE | MEM_F_SMARTPE;
		else if (c == 'M') opt->flag |= MEM_F_NO_MULTI;
		else if (c == 'S') opt->flag |= MEM_F_NO_RESCUE;
		else if (c == 'Y') opt->flag |= MEM_F_SOFTCLIP;
		else if (c == 'V') opt->flag |= MEM_F_REF_HDR;
		else if (c == '5') opt->flag |= MEM_F_PRIMARY5 | MEM_F_KEEP_SUPP_MAPQ; // always apply MEM_F_KEEP_SUPP_MAPQ with -5
		else if (c == 'q') opt->flag |= MEM_F_KEEP_SUPP_MAPQ;
		else if (c == 'u') opt->flag |= MEM_F_XB;
		else if (c == 'c') opt->max_occ = atoi(optarg), opt0.max_occ = 1;
		else if (c == 'd') opt->zdrop = atoi(optarg), opt0.zdrop = 1;
		else if (c == 'v') bwa_verbose = atoi(optarg);
		else if (c == 'j') ignore_alt = 1;
		else if (c == 'r') opt->split_factor = atof(optarg), opt0.split_factor = 1.;
		else if (c == 'D') opt->drop_ratio = atof(optarg), opt0.drop_ratio = 1.;
		else if (c == 'm') opt->max_matesw = atoi(optarg), opt0.max_matesw = 1;
		else if (c == 's') opt->split_width = atoi(optarg), opt0.split_width = 1;
		else if (c == 'G') opt->max_chain_gap = atoi(optarg), opt0.max_chain_gap = 1;
		else if (c == 'N') opt->max_chain_extend = atoi(optarg), opt0.max_chain_extend = 1;
		else if (c == 'o' || c == 'f') xreopen(optarg, "wb", stdout);
		else if (c == 'W') opt->min_chain_weight = atoi(optarg), opt0.min_chain_weight = 1;
		else if (c == 'y') opt->max_mem_intv = atol(optarg), opt0.max_mem_intv = 1;
		else if (c == 'C') aux.copy_comment = 1;
		else if (c == 'K') fixed_chunk_size = atoi(optarg);
		else if (c == 'X') opt->mask_level = atof(optarg);
		else if (c == 'h') {
			opt0.max_XA_hits = opt0.max_XA_hits_alt = 1;
			opt->max_XA_hits = opt->max_XA_hits_alt = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->max_XA_hits_alt = strtol(p+1, &p, 10);
		}
		else if (c == 'Q') {
			opt0.mapQ_coef_len = 1;
			opt->mapQ_coef_len = atoi(optarg);
			opt->mapQ_coef_fac = opt->mapQ_coef_len > 0? log(opt->mapQ_coef_len) : 0;
		} else if (c == 'O') {
			opt0.o_del = opt0.o_ins = 1;
			opt->o_del = opt->o_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->o_ins = strtol(p+1, &p, 10);
		} else if (c == 'E') {
			opt0.e_del = opt0.e_ins = 1;
			opt->e_del = opt->e_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->e_ins = strtol(p+1, &p, 10);
		} else if (c == 'L') {
			opt0.pen_clip5 = opt0.pen_clip3 = 1;
			opt->pen_clip5 = opt->pen_clip3 = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->pen_clip3 = strtol(p+1, &p, 10);
		} else if (c == 'R') {
			if ((rg_line = bwa_set_rg(optarg)) == 0) return 1; // FIXME: memory leak
		} else if (c == 'H') {
			if (optarg[0] != '@') {
				FILE *fp;
				if ((fp = fopen(optarg, "r")) != 0) {
					char *buf;
					buf = calloc(1, 0x10000);
					while (fgets(buf, 0xffff, fp)) {
						i = strlen(buf);
						assert(buf[i-1] == '\n'); // a long line
						buf[i-1] = 0;
						hdr_line = bwa_insert_header(buf, hdr_line);
					}
					free(buf);
					fclose(fp);
				}
			} else hdr_line = bwa_insert_header(optarg, hdr_line);
		} else if (c == 'I') { // specify the insert size distribution
			aux.pes0 = pes;
			pes[1].failed = 0;
			pes[1].avg = strtod(optarg, &p);
			pes[1].std = pes[1].avg * .1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].std = strtod(p+1, &p);
			pes[1].high = (int)(pes[1].avg + 4. * pes[1].std + .499);
			pes[1].low  = (int)(pes[1].avg - 4. * pes[1].std + .499);
			if (pes[1].low < 1) pes[1].low = 1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].high = (int)(strtod(p+1, &p) + .499);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].low  = (int)(strtod(p+1, &p) + .499);
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] mean insert size: %.3f, stddev: %.3f, max: %d, min: %d\n",
						__func__, pes[1].avg, pes[1].std, pes[1].high, pes[1].low);
		}
		else return 1;
	}
*/
	if (rg_line) {
		hdr_line = bwa_insert_header(rg_line, hdr_line);
		free(rg_line);
	}

//	if (opt->n_threads < 1) opt->n_threads = 1;
	opt->n_threads=1;
/*	if (optind + 1 >= argc || optind + 3 < argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]\n\n");
		fprintf(stderr, "Algorithm options:\n\n");
		fprintf(stderr, "       -t INT        number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "       -k INT        minimum seed length [%d]\n", opt->min_seed_len);
		fprintf(stderr, "       -w INT        band width for banded alignment [%d]\n", opt->w);
		fprintf(stderr, "       -d INT        off-diagonal X-dropoff [%d]\n", opt->zdrop);
		fprintf(stderr, "       -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [%g]\n", opt->split_factor);
		fprintf(stderr, "       -y INT        seed occurrence for the 3rd round seeding [%ld]\n", (long)opt->max_mem_intv);
//		fprintf(stderr, "       -s INT        look for internal seeds inside a seed with less than INT occ [%d]\n", opt->split_width);
		fprintf(stderr, "       -c INT        skip seeds with more than INT occurrences [%d]\n", opt->max_occ);
		fprintf(stderr, "       -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [%.2f]\n", opt->drop_ratio);
		fprintf(stderr, "       -W INT        discard a chain if seeded bases shorter than INT [0]\n");
		fprintf(stderr, "       -m INT        perform at most INT rounds of mate rescues for each read [%d]\n", opt->max_matesw);
		fprintf(stderr, "       -S            skip mate rescue\n");
		fprintf(stderr, "       -P            skip pairing; mate rescue performed unless -S also in use\n");
		fprintf(stderr, "\nScoring options:\n\n");
		fprintf(stderr, "       -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [%d]\n", opt->a);
		fprintf(stderr, "       -B INT        penalty for a mismatch [%d]\n", opt->b);
		fprintf(stderr, "       -O INT[,INT]  gap open penalties for deletions and insertions [%d,%d]\n", opt->o_del, opt->o_ins);
		fprintf(stderr, "       -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [%d,%d]\n", opt->e_del, opt->e_ins);
		fprintf(stderr, "       -L INT[,INT]  penalty for 5'- and 3'-end clipping [%d,%d]\n", opt->pen_clip5, opt->pen_clip3);
		fprintf(stderr, "       -U INT        penalty for an unpaired read pair [%d]\n\n", opt->pen_unpaired);
		fprintf(stderr, "       -x STR        read type. Setting -x changes multiple parameters unless overridden [null]\n");
		fprintf(stderr, "                     pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)\n");
		fprintf(stderr, "                     ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)\n");
		fprintf(stderr, "                     intractg: -B9 -O16 -L5  (intra-species contigs to ref)\n");
		fprintf(stderr, "\nInput/output options:\n\n");
		fprintf(stderr, "       -p            smart pairing (ignoring in2.fq)\n");
		fprintf(stderr, "       -R STR        read group header line such as '@RG\\tID:foo\\tSM:bar' [null]\n");
		fprintf(stderr, "       -H STR/FILE   insert STR to header if it starts with @; or insert lines in FILE [null]\n");
		fprintf(stderr, "       -o FILE       sam file to output results to [stdout]\n");
		fprintf(stderr, "       -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)\n");
		fprintf(stderr, "       -5            for split alignment, take the alignment with the smallest coordinate as primary\n");
		fprintf(stderr, "       -q            don't modify mapQ of supplementary alignments\n");
		fprintf(stderr, "       -K INT        process INT input bases in each batch regardless of nThreads (for reproducibility) []\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -v INT        verbosity level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
		fprintf(stderr, "       -T INT        minimum score to output [%d]\n", opt->T);
		fprintf(stderr, "       -h INT[,INT]  if there are <INT hits with score >80%% of the max score, output all in XA [%d,%d]\n", opt->max_XA_hits, opt->max_XA_hits_alt);
		fprintf(stderr, "       -a            output all alignments for SE or unpaired PE\n");
		fprintf(stderr, "       -C            append FASTA/FASTQ comment to SAM output\n");
		fprintf(stderr, "       -V            output the reference FASTA header in the XR tag\n");
		fprintf(stderr, "       -Y            use soft clipping for supplementary alignments\n");
		fprintf(stderr, "       -M            mark shorter split hits as secondary\n\n");
		fprintf(stderr, "       -I FLOAT[,FLOAT[,INT[,INT]]]\n");
		fprintf(stderr, "                     specify the mean, standard deviation (10%% of the mean if absent), max\n");
		fprintf(stderr, "                     (4 sigma from the mean if absent) and min of the insert size distribution.\n");
		fprintf(stderr, "                     FR orientation only. [inferred]\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "Note: Please read the man page for detailed description of the command line and options.\n");
		fprintf(stderr, "\n");
		free(opt);
		return 1;
	}*/
	//opt->T=0; //have no minimum alignment score
	if (mode) {
		if (strcmp(mode, "intractg") == 0) {
			if (!opt0.o_del) opt->o_del = 16;
			if (!opt0.o_ins) opt->o_ins = 16;
			if (!opt0.b) opt->b = 9;
			if (!opt0.pen_clip5) opt->pen_clip5 = 5;
			if (!opt0.pen_clip3) opt->pen_clip3 = 5;
		} else if (strcmp(mode, "pacbio") == 0 || strcmp(mode, "pbref") == 0 || strcmp(mode, "ont2d") == 0) {
			if (!opt0.o_del) opt->o_del = 1;
			if (!opt0.e_del) opt->e_del = 1;
			if (!opt0.o_ins) opt->o_ins = 1;
			if (!opt0.e_ins) opt->e_ins = 1;
			if (!opt0.b) opt->b = 1;
			if (opt0.split_factor == 0.) opt->split_factor = 10.;
			if (strcmp(mode, "ont2d") == 0) {
				if (!opt0.min_chain_weight) opt->min_chain_weight = 20;
				if (!opt0.min_seed_len) opt->min_seed_len = 14;
				if (!opt0.pen_clip5) opt->pen_clip5 = 0;
				if (!opt0.pen_clip3) opt->pen_clip3 = 0;
			} else {
				if (!opt0.min_chain_weight) opt->min_chain_weight = 40;
				if (!opt0.min_seed_len) opt->min_seed_len = 17;
				if (!opt0.pen_clip5) opt->pen_clip5 = 0;
				if (!opt0.pen_clip3) opt->pen_clip3 = 0;
			}
		} else {
			fprintf(stderr, "[E::%s] unknown read type '%s'\n", __func__, mode);
			return 1; // FIXME memory leak
		}
	} else update_a(opt, &opt0);
	bwa_fill_scmat(opt->a, opt->b, opt->mat);

	aux.idx = bwa_idx_load_from_shm(databaseFile);
	if (aux.idx == 0) {
		if ((aux.idx = bwa_idx_load(databaseFile, BWA_IDX_ALL)) == 0) return 1; // FIXME: memory leak
	} else if (bwa_verbose >= 3)
		fprintf(stderr, "[M::%s] load the bwa index from shared memory\n", __func__);
	if (ignore_alt)
		for (i = 0; i < aux.idx->bns->n_seqs; ++i)
			aux.idx->bns->anns[i].is_alt = 0;

	//ko = kopen(read1, &fd);
	//if (ko == 0) {
	//	if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, read1);
	//	return 1;
	//}
	//fp = gzdopen(fd, "r");
	//aux.ks = kseq_init(fp);
	//if (optind + 2 < argc) {
		//if (opt->flag&MEM_F_PE) {
		//	if (bwa_verbose >= 2)
		//		fprintf(stderr, "[W::%s] when '-p' is in use, the second query file is ignored.\n", __func__);
		//} else {
			//ko2 = kopen(read2, &fd2);
			//if (ko2 == 0) {
			//	if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__,read2);
			//	return 1;
			//}
			//fp2 = gzdopen(fd2, "r");
			//aux.ks2 = kseq_init(fp2);
			//opt->flag |= MEM_F_PE;
	//	}
	//}
	if (paired == 0 ){
		opt->flag |= MEM_F_NOPAIRING;
	}else{
		opt->flag |= MEM_F_PE;
	}
	//opt->flag |= MEM_F_ALL;
	//aux.names = names;
	//aux.read1 = read1;
	//aux.read2 = read2;
	aux.max_query_length = max_query_length;
	aux.max_readname_length = max_readname_length;
	aux.max_acc_name = max_acc_name;
	aux.number_of_seqs = number_of_seqs;
	aux.results = bwa_results;
	aux.concordant = concordant;
	aux.ntree = numberOfTrees;
	aux.startline = startline;
	aux.paired = paired;
	aux.start = start;
	aux.end = end;
	//bwa_print_sam_hdr(aux.idx->bns, hdr_line);
	aux.actual_chunk_size = fixed_chunk_size > 0? fixed_chunk_size : opt->chunk_size * opt->n_threads;
	//aux.actual_chunk_size = number_of_seqs;
	//kt_pipeline(no_mt_io? 1 : 2, process, &aux, 3);
	kt_pipeline(1, process, &aux, 3);
	bwa_results = aux.results;
	free(hdr_line);
	free(opt);
	bwa_idx_destroy(aux.idx);
	kseq_destroy(aux.ks);
	//err_gzclose(fp); kclose(ko);
	if (aux.ks2) {
		kseq_destroy(aux.ks2);
		//err_gzclose(fp2); kclose(ko2);
	}
	return 0;
}
/*
int main_fastmap(int argc, char *argv[])
{
	int c, i, min_iwidth = 20, min_len = 17, print_seq = 0, min_intv = 1, max_len = INT_MAX;
	uint64_t max_intv = 0;
	kseq_t *seq;
	bwtint_t k;
	gzFile fp;
	smem_i *itr;
	const bwtintv_v *a;
	bwaidx_t *idx;

	while ((c = getopt(argc, argv, "w:l:pi:I:L:")) >= 0) {
		switch (c) {
			case 'p': print_seq = 1; break;
			case 'w': min_iwidth = atoi(optarg); break;
			case 'l': min_len = atoi(optarg); break;
			case 'i': min_intv = atoi(optarg); break;
			case 'I': max_intv = atol(optarg); break;
			case 'L': max_len  = atoi(optarg); break;
		    default: return 1;
		}
	}
	if (optind + 1 >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa fastmap [options] <idxbase> <in.fq>\n\n");
		fprintf(stderr, "Options: -l INT    min SMEM length to output [%d]\n", min_len);
		fprintf(stderr, "         -w INT    max interval size to find coordiantes [%d]\n", min_iwidth);
		fprintf(stderr, "         -i INT    min SMEM interval size [%d]\n", min_intv);
		fprintf(stderr, "         -L INT    max MEM length [%d]\n", max_len);
		fprintf(stderr, "         -I INT    stop if MEM is longer than -l with a size less than INT [%ld]\n", (long)max_intv);
		fprintf(stderr, "\n");
		return 1;
	}

	fp = xzopen(argv[optind + 1], "r");
	seq = kseq_init(fp);
	if ((idx = bwa_idx_load(argv[optind], BWA_IDX_BWT|BWA_IDX_BNS)) == 0) return 1;
	itr = smem_itr_init(idx->bwt);
	smem_config(itr, min_intv, max_len, max_intv);
	while (kseq_read(seq) >= 0) {
		err_printf("SQ\t%s\t%ld", seq->name.s, seq->seq.l);
		if (print_seq) {
			err_putchar('\t');
			err_puts(seq->seq.s);
		} else err_putchar('\n');
		for (i = 0; i < seq->seq.l; ++i)
			seq->seq.s[i] = nst_nt4_table[(int)seq->seq.s[i]];
		smem_set_query(itr, seq->seq.l, (uint8_t*)seq->seq.s);
		while ((a = smem_next(itr)) != 0) {
			for (i = 0; i < a->n; ++i) {
				bwtintv_t *p = &a->a[i];
				if ((uint32_t)p->info - (p->info>>32) < min_len) continue;
				err_printf("EM\t%d\t%d\t%ld", (uint32_t)(p->info>>32), (uint32_t)p->info, (long)p->x[2]);
				if (p->x[2] <= min_iwidth) {
					for (k = 0; k < p->x[2]; ++k) {
						bwtint_t pos;
						int len, is_rev, ref_id;
						len  = (uint32_t)p->info - (p->info>>32);
						pos = bns_depos(idx->bns, bwt_sa(idx->bwt, p->x[0] + k), &is_rev);
						if (is_rev) pos -= len - 1;
						bns_cnt_ambi(idx->bns, pos, len, &ref_id);
						err_printf("\t%s:%c%ld", idx->bns->anns[ref_id].name, "+-"[is_rev], (long)(pos - idx->bns->anns[ref_id].offset) + 1);
					}
				} else err_puts("\t*");
				err_putchar('\n');
			}
		}
		err_puts("//");
	}

	smem_itr_destroy(itr);
	bwa_idx_destroy(idx);
	kseq_destroy(seq);
	err_gzclose(fp);
	return 0;
}*/
