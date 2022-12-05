/*
 * global.h
 */
#ifndef _GLOBAL_
#define _GLOBAL_

#define STATESPACE 20 /*number of categories in approximation of gamma distribution for Ne. Must be at least 4 because some of the memory is used for the nucleotide model*/
#define MAXNUMBEROFINDINSPECIES 500 /*maximum number of individuals belonging to a species*/
#define MAXQUERYLENGTH 1000 /*the maximum length of the query seqeunce*/
#define NUMCAT 1/*number of categories in the discretization of the gamma for the nucleotide substituion model*/
#define MINBL 0.000001
#define MAXBL 2.0
#define MAXTIEBREAK 64
#define type_of_PP double
#define MAX_NODENAME 30
#define MAX_NUMBEROFROOTS 10000
#define MAXRESULTSNAME 2000
#define MAXREADNAME 300
#define MAX_NUM_BWA_MATCHES 5000
#define SP_SCORE_MIN 0.8
#define FASTA_MAXLINE 50000
#define NUM_THREADS 1
#define MAXTAXLENGTHNAME 256
#define MAXFILENAME 1000
#define BUFFER_SIZE 1000
extern FILE *infile, *outfile, *treefile;
extern int numspec, numbase, **seq, numundspec[MAXNUMBEROFINDINSPECIES+1];
extern int *rootArr, *numspecArr, *numbaseArr, ***seqArr;
extern int root,tip,comma; /*globals used to read in the tree*/
extern double Logfactorial[MAXNUMBEROFINDINSPECIES];
extern double LRVEC[STATESPACE][STATESPACE], RRVEC[STATESPACE][STATESPACE], RRVAL[STATESPACE], PMAT1[STATESPACE][STATESPACE], PMAT2[STATESPACE][STATESPACE];
extern double LRVECnc[4][4], RRVECnc[4][4], RRVALnc[4], PMATnc[2][4][5];
extern double *statevector, UFC, *UFCnc, **templike_nc;
//extern int ***PP; // WE SHOULD TEST WHAT IS THE FASTEST. MAYBE MOVE TO UNSIGNED SHORT?
extern type_of_PP ***PP;
extern type_of_PP ***PPcopy;
extern type_of_PP ****PP_Arr;
extern char ***taxonomy;
extern char ****taxonomyArr;

typedef struct node{
	int up[2];
	int down;
	int nd;
	int depth;
	double bl;
	double *like;
	double **likenc;
	double *posterior;
	double **posteriornc;
	int s; /*number of segregating sites*/
	int numsites; /*total number of sites summed over all individuals in specie. numsites/numbase is number of individuals if there are no missing data*/
	int spec; /*indicates which species the node belongs to. -1 if no single species*/
	int mrca; /*is 1 if node is a MFRCA for a spcies*/
	char *name;
	int taxIndex[2];
	type_of_PP score;
	//int nw_fail;
	int SP_fail;
}node;

typedef struct queryMatPaired{
	char **query1Mat;
	char **query2Mat;
	char **name;
}queryMatPaired;

typedef struct queryMatSingle{
	char **queryMat;
	char **name;
}queryMatSingle;

extern queryMatPaired *pairedQueryMat;
extern queryMatSingle *singleQueryMat;

typedef struct bwaMatches{
	char *readname;
	int *concordant_matches;
	int *discordant_matches;
	char **concordant_leaf_matches;
	char **discordant_leaf_matches;
	int n_matches;
}bwaMatches;

typedef struct scoresStruct{
	type_of_PP score1;
	type_of_PP score2;
	int	rootNumber;
	int nodeNumber;
}scoresStruct;

typedef struct partition_files{
	char **tree_files;
	char **msa_files;
	char **tax_files;
}partition_files;

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
	int use_spscore;
	int use_min_leaves;
	int min_leaves;
	double sp_score;
	char fasta_file[200];
	double cinterval;
	char readdir[2000];
	int number_of_partitions;
	int restart;
	char rmrefdir[2000];
	int missing_data;
}Options;

typedef struct masterArr{
	char index[10];
	node *tree;
	int **msa;
	char ***taxonomy;
	int numspec;
	int root;
	int numbase;
	char **names;
}masterArr;
//extern node *treeArr[MAX_NUMBEROFROOTS];
extern node **treeArr;
extern int COUNT2;
extern int COUNT;
extern double *localpi;
extern int localnode;
extern double currentestimate[10];

//extern char *locQuery;
//extern int positions[MAXQUERYLENGTH];
extern char *seqInRoot;
extern char **nodeIDs;
extern char ***nodeIDsArr;
extern double parameters[10];
extern double maxpar[4];
extern double ml;
//extern int open;
extern int closed;
extern type_of_PP Cinterval;
extern double *nw_scores;
extern int *nodesToCut;
extern double minVariance;
extern int minVarNode;
extern int *nodesToCutMinVar;
extern int *partitionSizes;
extern int returnNode;
//extern int SPscoreArr[MAX_NUMBEROFROOTS];
extern int *SPscoreArr;
//extern struct hashmap map;
//extern struct hashmap leaf_map;
#endif
