/*
 * global.h
 */
#ifndef _GLOBAL_
#define _GLOBAL_

#include <stdio.h>

#define STATESPACE                                                             \
    20 /*number of categories in approximation of gamma distribution for Ne.   \
          Must be at least 4 because some of the memory is used for the        \
          nucleotide model*/
#define MAXNUMBEROFINDINSPECIES                                                \
    500 /*maximum number of individuals belonging to a species*/
#define MAXQUERYLENGTH 1000 /*the maximum length of the query seqeunce*/
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
extern int *rootArr, *numspecArr, *numbaseArr, ***seqArr;
extern int tip, comma; /*globals used to read in the tree*/
extern double Logfactorial[MAXNUMBEROFINDINSPECIES];
extern double LRVEC[STATESPACE][STATESPACE], RRVEC[STATESPACE][STATESPACE],
    RRVAL[STATESPACE], PMAT1[STATESPACE][STATESPACE],
    PMAT2[STATESPACE][STATESPACE];
extern double **templike_nc;
extern char ***taxonomy;
extern char ****taxonomyArr;

struct node {
    int up[1000];
    int down;
    int nd;
    int depth;
    double bl;
    double *like;
    double **likenc;
    double *posterior;
    double **posteriornc;
    int s;        /*number of segregating sites*/
    int numsites; /*total number of sites summed over all individuals in specie.
                     numsites/numbase is number of individuals if there are no
                     missing data*/
    int spec;     /* Indicates which species the node belongs to. Species are
                     numbered 1, 2, ...; the value 0 is unused. If both children of
                     a node are of species X, then the parent node will also be
                     species X. This field will be -1 if some descendants of
                     the node belong to different species. The value -2 means
                     "uninitialized". */
    int mrca;     /*is 1 if node is a MRCA for a spcies*/
    char *name;
    int taxIndex[2];
    type_of_PP score;
    // int nw_fail;
    int SP_fail;
};

typedef struct queryMatPaired {
    char **query1Mat;
    char **query2Mat;
    char **name;
} queryMatPaired;

typedef struct queryMatSingle {
    char **queryMat;
    char **name;
} queryMatSingle;

extern queryMatPaired *pairedQueryMat;
extern queryMatSingle *singleQueryMat;

typedef struct bwaMatches {
    char *readname;
    int *concordant_matches;
    int *discordant_matches;
    char **concordant_leaf_matches;
    char **discordant_leaf_matches;
    int n_matches;
} bwaMatches;

typedef struct scoresStruct {
    type_of_PP score1;
    type_of_PP score2;
    int rootNumber;
    int nodeNumber;
} scoresStruct;

typedef struct partition_files {
    char **tree_files;
    char **msa_files;
    char **tax_files;
} partition_files;

typedef struct Options {
    char msa_file[200];
    char tree_file[200];
    char taxonomy_file[200];
    int number_of_trees;
    int two_step;
    int remove_unused;
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
    char prefix[2000];
    int missing_data;
    int famsa_threads;
    int fasttree;
} Options;

typedef struct masterArr {
    char index[10];
    struct node *tree;
    int **msa;
    char ***taxonomy;
    int numspec;
    int numNodes;
    int root;
    int numbase;
    char **names;
    char *filename;
} masterArr;


extern struct node **treeArr;
extern double currentestimate[10];

extern double minVariance;
extern int minVarNode;
extern int *SPscoreArr;
#endif
