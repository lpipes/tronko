#ifndef _READ_REF_
#define _READ_REF_

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>
#include "global.h"
#include "opt.h"
#include "allocatetreememory.h"
#include <string.h>
#include <dirent.h>
#include <regex.h>
int compare_strings(const void* a, const void* b);
int compare_natural(const void *a, const void *b);
void readFilesInDir(char *directory, int number_of_partitions, partition_files *pf);
int readInXNumberOfLines(int numberOfLinesToRead, gzFile query_reads, int whichPair);
int countLinesInFile(FILE *queryFile);
void readSeqArrFull(FILE *seqArrFile);
void readPP_Arr(FILE *PP_ArrFile);
int readReferenceTree( gzFile referenceTree);
void replaceTrees(int*** seqArr_heap, char**** taxonomyArr_heap, char*** nodeIDsArr_heap,Options opt);
#endif /* _READ_REF_ */
