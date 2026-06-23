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

int compare_natural(const void *a, const void *b);
void readFilesInDir(char *directory, int number_of_partitions,
                    partition_files *pf);

#endif /* _READ_REF_ */
