#ifndef _PRINTTREE_
#define _PRINTTREE_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "global.h"
#include "opt.h"

void printtree();
void printtreeArr();
void printTreeFile(int numberOfTrees, int max_nodename, int max_tax_name, int max_lineTaxonomy, Options opt);
void printTaxonomyArrToFile(int numberOfTrees);
void printTaxonomyToFile();
void printTreeToFile();
void printSeqArr(int numberOfTrees);
void printSeqToFile();

#endif /* _PRINTTREE_ */
