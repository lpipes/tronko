#include "printtree.h"

// Prints a tree specified by `whichRoot` in human-readable format (for
// debugging)
void printtreeArr(int whichRoot) {
    int i, j;
    printf("\n");
    j = whichRoot;
    for (i = 0; i < 2 * numspecArr[j] - 1; i++) {
        struct node *cnode = &(treeArr[j][i]);
        printf("Node %5i: up: (%5i, %5i), down: %5i, spec=%4d, mrca=%4d, "
               "s=%2d, taxIndex[0]: %5d, taxIndex[1]: %1d, depth: "
               "%2d, nd: %d, name: %s, ",
               i, cnode->up[0], cnode->up[1], cnode->down, cnode->spec,
               cnode->mrca, cnode->s, cnode->taxIndex[0], cnode->taxIndex[1],
               cnode->depth, cnode->nd, cnode->name);
        if (cnode->taxIndex[0] != -1 && cnode->taxIndex[1] != -1)
            printf("taxonomy: %s",
                   taxonomyArr[j][cnode->taxIndex[0]][cnode->taxIndex[1]]);
        else printf("taxonomy: Domain");

        if (i != rootArr[j]) printf(" (bl: %f)", cnode->bl);
        printf("\n");
    }
}

void printTreeFile(int numberOfTrees, int max_nodename, int max_tax_name,
                   int max_lineTaxonomy, Options opt) {
    int i, j, k, l;
    char buf[BUFFER_SIZE];
    struct stat st = {0};
    if (stat(opt.partitions_directory, &st) == -1) {
        mkdir(opt.partitions_directory, 0700);
    }
    snprintf(buf, BUFFER_SIZE, "%s/reference_tree.txt",
             opt.partitions_directory);
    FILE *outputTree = fopen(buf, "w");
    if (outputTree == NULL) {
        printf("Error opening reference tree file!\n");
        exit(-1);
    }
    fprintf(outputTree, "%d\n", numberOfTrees);
    fprintf(outputTree, "%d\n", max_nodename);
    fprintf(outputTree, "%d\n", max_tax_name);
    fprintf(outputTree, "%d\n", max_lineTaxonomy);
    for (i = 0; i < numberOfTrees; i++) {
        fprintf(outputTree, "%d\t%d\t%d\n", numbaseArr[i], rootArr[i],
                numspecArr[i]);
    }
    for (i = 0; i < numberOfTrees; i++) {
        for (j = 0; j < numspecArr[i]; j++) {
            for (k = 0; k < 7; k++) {
                if (k == 6) {
                    fprintf(outputTree, "%s\n", taxonomyArr[i][j][k]);
                } else {
                    fprintf(outputTree, "%s;", taxonomyArr[i][j][k]);
                }
            }
        }
    }
    for (i = 0; i < numberOfTrees; i++) {
        for (j = 0; j < 2 * numspecArr[i] - 1; j++) {
            fprintf(outputTree, "%d\t", i);
            fprintf(outputTree, "%d\t", j);
            fprintf(outputTree, "%d\t", treeArr[i][j].up[0]);
            fprintf(outputTree, "%d\t", treeArr[i][j].up[1]);
            fprintf(outputTree, "%d\t", treeArr[i][j].down);
            fprintf(outputTree, "%d\t", treeArr[i][j].depth);
            fprintf(outputTree, "%d\t", treeArr[i][j].taxIndex[0]);
            fprintf(outputTree, "%d\t", treeArr[i][j].taxIndex[1]);
            if (treeArr[i][j].up[0] != -1 && treeArr[i][j].up[1] != -1) {
                fprintf(outputTree, "\n");
            } else {
                fprintf(outputTree, "%s\n", treeArr[i][j].name);
            }
            for (k = 0; k < numbaseArr[i] - 1; k++) {
                for (l = 0; l < 3; l++) {
                    fprintf(outputTree, "%.17g\t",
                            treeArr[i][j].posteriornc[k][l]);
                }
                fprintf(outputTree, "%.17g\n", treeArr[i][j].posteriornc[k][3]);
            }
            fprintf(outputTree, "%.17g\t%.17g\t%.17g\t%.17g\n",
                    treeArr[i][j].posteriornc[numbaseArr[i] - 1][0],
                    treeArr[i][j].posteriornc[numbaseArr[i] - 1][1],
                    treeArr[i][j].posteriornc[numbaseArr[i] - 1][2],
                    treeArr[i][j].posteriornc[numbaseArr[i] - 1][3]);
        }
    }
    fclose(outputTree);
}
