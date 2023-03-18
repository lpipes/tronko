#ifndef _GETSEQINROOT_
#define _GETSEQINROOT_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "global.h"

int getStartPosition(int start, int rootNum, int node, int padding);
int getEndPosition(char* s, int rootNum, int node, int start, int padding);
void getSequenceInNode(int rootNum, int node, char *sequenceInNode);
void getSequenceInNodeWithoutNs(int rootNum, int node, char *sequenceInNode, int *positionsInRoot, int start_position, int end_position);
#endif /* _GETSEQINROOT_ */
