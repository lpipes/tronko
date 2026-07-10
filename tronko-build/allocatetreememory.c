#include <string.h>
#include "allocatetreememory.h"

void allocatetreememory_for_nucleotide_Arr(int numberOfTrees) {
    int i, j, k;
    for (k = 0; k < numberOfTrees; k++) {
        for (i = 0; i < (numspecArr[k] * 2 - 1); i++) {
            treeArr[k][i].likenc = malloc(numbaseArr[k] * (sizeof(double *)));
            treeArr[k][i].posteriornc =
                malloc(numbaseArr[k] * (sizeof(double *)));
            for (j = 0; j < numbaseArr[k]; j++) {
                treeArr[k][i].likenc[j] = malloc(4 * (sizeof(double)));
                treeArr[k][i].posteriornc[j] = malloc(4 * (sizeof(double)));
            }
        }
    }
}

void allocateTreeArrMemory(struct masterArr *m, int max_nodename) {
    int i;
    for (i = 0; i < m->numspec * 2 - 1; i++) {
        m->tree[i].like = malloc(STATESPACE * (sizeof(double)));
        m->tree[i].posterior = malloc(STATESPACE * (sizeof(double)));
        m->tree[i].name = malloc((max_nodename + 1) * sizeof(char));
        m->tree[i].up[0] = -2; // -2 is uninitialized / NULL value
        m->tree[i].up[1] = -2;
        m->tree[i].down = -2;
        m->tree[i].nd = -2;
        m->tree[i].depth = -2;
        m->tree[i].bl = -2.0;
        m->tree[i].likenc = NULL;
        m->tree[i].posteriornc = NULL;
        m->tree[i].s = -2;
        m->tree[i].numsites = -2;
        m->tree[i].spec = -2;
        m->tree[i].mrca = -2;
        memset(m->tree[i].name, '\0', max_nodename + 1); // name is 20 bytes
    }
}

void freeTreeMemory(int whichPartition) {
    int i, j;
    for (i = 0; i < 2 * numspecArr[whichPartition] - 1; i++) {
        for (j = 0; j < numbaseArr[whichPartition]; j++) {
            free(treeArr[whichPartition][i].likenc[j]);
            free(treeArr[whichPartition][i].posteriornc[j]);
        }
        free(treeArr[whichPartition][i].like);
        free(treeArr[whichPartition][i].posterior);
        free(treeArr[whichPartition][i].likenc);
        free(treeArr[whichPartition][i].posteriornc);
        free(treeArr[whichPartition][i].name);
    }
    free(treeArr[whichPartition]);
}

// by default calloc returns a pointer to allocated memory that has been
// initialized to zero (vs malloc, which just returns a pointer to allocated
// memory and doesn't initialize it), but it's good practice to check that the
// returned pointer is valid, since calloc (and malloc) will return NULL if
// the system is out of memory.
void *calloc_check(size_t nmemb, size_t size) {
    void *ptr = calloc(nmemb, size);
    if (!ptr) { // if calloc returns NULL, we've run out of memory
        fprintf(stderr, "Out of memory!\n");
        exit(1);
    }
    return ptr;
}
