/*
 * global.h
 */
#ifndef _GLOBAL_
#define _GLOBAL_
#define STATESPACE 20 /*number of categories in approximation of gamma distribution for Ne. Must be at least 4 because some of the memory is used for the nucleotide model*/
#define MAXNUMBEROFINDINSPECIES 500 /*maximum number of individuals belonging to a species*/
#define MAXQUERYLENGTH 30000
#define MAX_CIGAR 100
#define NUMCAT 1/*number of categories in the discretization of the gamma for the nucleotide substituion model*/
#define MINBL 0.000001
#define MAXBL 2.0
#define MAXTIEBREAK 64
#define type_of_PP double
#define MAX_NODENAME 30
#define MAX_NUMBEROFROOTS 20000
#define MAXRESULTSNAME 2000
#define MAXREADNAME 300
#define MAX_NUM_BWA_MATCHES 10
#define SP_SCORE_MIN 0.8
#define FASTA_MAXLINE 40000
#define MAXTAXLENGTHNAME 256
#define MAXFILENAME 1000
#define BUFFER_SIZE 1000
#define BIG_BUFFER 409600
//extern FILE *infile, *outfile, *treefile;
//extern int numspec, numbase, root/*, **seq, numundspec[MAXNUMBEROFINDINSPECIES+1]*/;
extern int *rootArr, *numspecArr, *numbaseArr;
//extern int root,tip,comma; /*globals used to read in the tree*/
//extern double Logfactorial[MAXNUMBEROFINDINSPECIES];
//extern double LRVEC[STATESPACE][STATESPACE], RRVEC[STATESPACE][STATESPACE], RRVAL[STATESPACE], PMAT1[STATESPACE][STATESPACE], PMAT2[STATESPACE][STATESPACE];
//extern double LRVECnc[4][4], RRVECnc[4][4], RRVALnc[4], PMATnc[2][4][5];
//extern double *statevector, UFC, *UFCnc, **templike_nc;
//extern int ***PP; // WE SHOULD TEST WHAT IS THE FASTEST. MAYBE MOVE TO UNSIGNED SHORT?
//extern type_of_PP ***PP;
//extern type_of_PP ***PPcopy;
//extern type_of_PP ****PP_Arr;
//extern char ***taxonomy;
extern char ****taxonomyArr;

#include "needleman_wunsch.h"
#include "hashmap.h"
#include "hashmap_base.h"
typedef struct node{
	int up[2];
	int down;
	int nd;
	int depth;
	double **posteriornc;
	char *name;
	int taxIndex[2];
}node;

typedef struct leafMap{
	char* name;
	int root;
	int node;
}leafMap;
typedef struct ParserState{
	int treeNumber;
	int nodeNumber;
	int i; // position in numbase
	int isIncompleteLine;
} ParserState;

typedef struct queryMatPaired{
	char **query1Mat;
	char **query2Mat;
	char **forward_name;
	char **reverse_name;
}queryMatPaired;

typedef struct queryMatSingle{
	char **queryMat;
	char **name;
}queryMatSingle;

extern queryMatPaired *pairedQueryMat;
extern queryMatSingle *singleQueryMat;

typedef struct resultsStruct{
	nw_aligner_t *nw;
	alignment_t *aln;
	scoring_t *scoring;
	type_of_PP ***nodeScores;
	int **voteRoot;
	int *positions;
	char *locQuery;
	char **taxonPath;
	char **LCAnames;
	int *minNodes;
	int **leaf_coordinates;
	type_of_PP *minimum;
	int print_alignments;
	int *starts_forward;
	int *starts_reverse;
	char **cigars_forward;
	char **cigars_reverse;
}resultsStruct;

typedef struct mystruct{
	char **rootSeqs;
	int ntree;
	int start;
	int end;
	int paired;
	resultsStruct *str;
	int concordant;
	char* databasefile;
	char* alignmentsdir;
	int maxNumSpec;
	int numspec_total;
	int use_nw;
	int print_unassigned;
	int print_alignments_to_file;
	int use_leaf_portion;
	int padding;
	int max_query_length;
	int max_readname_length;
	int max_acc_name;
	int max_numbase;
	int max_lineTaxonomy;
	int number_of_total_nodes;
	int print_all_nodes;
}mystruct;

typedef struct bwaMatches{
	//char *readname;
	int *concordant_matches_roots;
	int *concordant_matches_nodes;
	int *discordant_matches_roots;
	int *discordant_matches_nodes;
	//char **concordant_leaf_matches;
	//char **discordant_leaf_matches;
	int n_matches;
	char **cigars_forward;
	int *starts_forward;
	char **cigars_reverse;
	int *starts_reverse;
	int use_portion;
}bwaMatches;

typedef struct scoresStruct{
	type_of_PP score1;
	type_of_PP score2;
	int	rootNumber;
	int nodeNumber;
}scoresStruct;

typedef struct Options{
	char msa_file[200];
	char tree_file[200];
	char taxonomy_file[200];
	int number_of_trees;
	int reference_mode;
	int use_partitions;
	char reference_file[200];
	char paired_or_single[7];
	char read1_file[2000];
	char read2_file[2000];
	char partitions_directory[200];
	char results_file[200];
	double sp_score;
	char fasta_file[200];
	double cinterval;
	char readdir[2000];
	char print_trees_dir[2000];
	int number_of_partitions;
	int reverse_single_read;
	int reverse_second_of_paired_read;
	char print_alignments_dir[BUFFER_SIZE];
	char treedir[BUFFER_SIZE];
	int number_of_cores;
	int number_of_lines_to_read;
	int print_alignments;
	int use_nw;
	int fastq;
	int unassigned;
	int print_alignments_to_file;
	int use_leaf_portion;
	int padding;
	char print_node_info[BUFFER_SIZE];
	int skip_build;
	int print_leave_seqs;
	double score_constant;
	int print_all_nodes;
}Options;

//extern node *tree;
//extern node *treeArr[MAX_NUMBEROFROOTS];
extern node **treeArr;
//extern struct hashmap_base *map;
/*#define HASHMAP(key_type, data_type)                                    \
    struct {                                                            \
        struct hashmap_base map_base;                                   \
        struct {                                                        \
            const key_type *t_key;                                      \
            data_type *t_data;                                          \
            size_t (*t_hash_func)(const key_type *);                    \
            int (*t_compare_func)(const key_type *, const key_type *);  \
            key_type *(*t_key_dup_func)(const key_type *);              \
            void (*t_key_free_func)(key_type *);                        \
            int (*t_foreach_func)(const key_type *, data_type *, void *); \
            struct {                                                    \
                struct hashmap_base *iter_map;                          \
                struct hashmap_entry *iter_pos;                         \
                struct {                                                \
                    const key_type *t_key;                              \
                    data_type *t_data;                                  \
                } iter_types[0];                                        \
            } t_iterator;                                               \
        } map_types[0];                                                 \
    }*/
//extern HASHMAP(char, struct leafMap) map;
//extern struct map;
//extern struct hashmap_base map;
//extern node *rootTree;
//extern int COUNT2;
//extern int COUNT;
//extern double *localpi;
//extern int localnode;
//extern double currentestimate[10];

//extern char *locQuery;
//extern int positions[MAXQUERYLENGTH];
//extern char *seqInRoot;
//extern char **nodeIDs;
//extern char ***nodeIDsArr;
//extern double parameters[10];
//extern double maxpar[4];
//extern double ml;
//extern int open;
//extern int closed;
extern type_of_PP Cinterval;
//extern double *nw_scores;
//extern int *nodesToCut;
//extern double minVariance;
//extern int minVarNode;
//extern int *nodesToCutMinVar;
//extern int *partitionSizes;
//extern int returnNode;
//extern int SPscoreArr[MAX_NUMBEROFROOTS];
//extern leafMap *leaf_map;
//extern struct hashmap map;
//extern struct hashmap leaf_map;
#endif
