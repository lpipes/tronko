#include "readreference.h"

#include <ctype.h>
#include <string.h>

int compare_natural(const void *a, const void *b) {
    const char *ia = *(const char **)a;
    const char *ib = *(const char **)b;
    while (*ia && *ib) {
        if (isdigit(*ia) && isdigit(*ib)) {
            // Compare the numbers
            long na = strtol(ia, (char **)&ia, 10);
            long nb = strtol(ib, (char **)&ib, 10);
            if (na < nb) {
                return -1;
            } else if (na > nb) {
                return 1;
            }
        } else {
            // Compare the characters
            if (*ia < *ib) {
                return -1;
            } else if (*ia > *ib) {
                return 1;
            }
            ++ia;
            ++ib;
        }
    }
    if (*ia) {
        return 1;
    } else if (*ib) {
        return -1;
    } else {
        return 0;
    }
}

void readFilesInDir(char *directory, int number_of_partitions,
                    partition_files *pf) {
    struct dirent *de;
    DIR *dr = opendir(directory);
    regex_t regex1;
    char msgbuf[100];
    int i, j = 0;
    char name[300];
    char name2[300];
    char name3[300];
    char **msa_files =
        (malloc)(number_of_partitions *
                 sizeof(char *)); // intermediate array to store filenames
    for (i = 0; i < number_of_partitions; i++) {
        msa_files[i] = (malloc)(300 * sizeof(char));
        for (j = 0; j < 300; j++) {
            msa_files[i][j] = '\0';
        }
    }
    int reti = regcomp(&regex1, "_MSA\\.fasta$", 0);
    if (reti) {
        fprintf(stderr, "Could not compile regex\n");
        exit(-1);
    }

    if (dr == NULL) {
        printf("Could not open directory for reads");
        exit(-1);
    }

    int file_count = 0;
    while ((de = readdir(dr)) != NULL) {
        reti = regexec(&regex1, de->d_name, 0, NULL, 0);
        if (!reti) {
            msa_files[file_count] =
                strdup(de->d_name); // store filenames in the array
            file_count++;
        } else if (reti != REG_NOMATCH) {
            regerror(reti, &regex1, msgbuf, sizeof(msgbuf));
            fprintf(stderr, "Regex match failed: %s\n", msgbuf);
            exit(1);
        }
    }

    closedir(dr);

    // Sort the msa_files
    qsort(msa_files, file_count, sizeof(char *), compare_natural);
    j = 0;
    // Now loop over the sorted list and process each file
    for (int k = 0; k < file_count; k++) {
        strcpy(pf->msa_files[j], msa_files[k]);

        for (i = 0; i < 300; i++) {
            name[i] = '\0';
        }

        for (i = 0; i < (int)strlen(msa_files[k]) - 10; i++) {
            name[i] = msa_files[k][i];
        }

        name[i] = '\0';
        strcpy(name2, name);
        strncat(name2, "_taxonomy.txt", 250);
        strcpy(pf->tax_files[j], name2);

        strcpy(name3, "RAxML_bestTree.");
        strcat(name3, name);
        strcat(name3, ".reroot");
        strcpy(pf->tree_files[j], name3);

        j++;
        free(msa_files[k]);
    }
    free(msa_files);
    regfree(&regex1);
}
