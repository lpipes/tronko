#ifndef _READ_FASTA_
#define _READ_FASTA_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "global.h"
#include "opt.h"
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <errno.h>                      
#include <fcntl.h>
#include <string.h>
#include <zlib.h>

void setNumspec(gzFile infile, int* specs);
int setNumspecArr(gzFile partitionsFile);
//void readSeqArr(FILE *partitionsFile, int whichPartitions, int maxname);
//void readSeqArr_UsePartitions(FILE *partitionsFile, int whichPartitions,int*** seqArr_heap, char*** nodeIDsArr_heap);
FILE *openFasta(char *infile);
void readfasta(FILE *infile);
void readseq(gzFile infile, int max_nodename, struct masterArr *m);

#endif /* _READ_FASTA */
