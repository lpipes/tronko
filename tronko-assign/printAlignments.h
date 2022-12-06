#ifndef _PRINTALIGN_
#define _PRINTALIGN_

#include <stdio.h>
#include "global.h"
#include <unistd.h>

void createNewFile( char *referenceName, char *alignments_directory, alignment_t *aln, char *readname);
void printToFile (char *referenceName, char *alignments_directory, alignment_t *aln, char *readname);
void createNewFile2( char *referenceName, char *alignments_directory, alignment_t *aln, char *readname);
void printToFile2 (char *referenceName, char *alignments_directory, alignment_t *aln, char *readname, int length_of_read, int* positionsInRoot);
void printLeaveSeqsToFile(FILE *output, int node, int whichTree, int numbase);
void printTreesNumbersToFile(FILE *output, int node, int whichTree);
#endif /* _PRINTALIGN_ */
