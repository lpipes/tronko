#ifndef _READ_REF_
#define _READ_REF_

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>
#include "global.h"
#include <string.h>
#include <dirent.h>
#include <regex.h>
//void readFilesInDir(char *directory, int number_of_partitions, partition_files *pf);
int readInXNumberOfLines_fastq(int numberOfLinesToRead, gzFile query_reads, int whichPair, Options opt, int max_query_length, int max_readname_length,int first_iter);
int readInXNumberOfLines(int numberOfLinesToRead, gzFile query_reads, int whichPair, Options opt, int max_query_length, int max_readname_length);
int readReferenceTree( gzFile referenceTree,int* name_specs);
int setNumbase_setNumspec(int numberOfPartitions,int* specs);
#endif /* _READ_REF_ */
