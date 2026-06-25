/*
 * readsam.h
 *
 * Read a user-supplied SAM/BAM alignment file in place of running bwa-mem
 * internally. The SAM records map query reads to reference accessions; this
 * module resolves those accessions to tree (root,node) coordinates and fills
 * the same bwaMatches structure that the bwa pipeline produces, so the
 * downstream placement code is unchanged.
 *
 * To scale to arbitrarily many reads, the file is *streamed* in lockstep with
 * tronko's read batches: only the records for the reads in the current batch
 * are held in memory at once. This requires the alignment file to be in the
 * same order as the read file(s) -- i.e. the natural, unsorted output of the
 * aligner (e.g. `bwa mem reads.fq`), NOT a coordinate- or name-sorted BAM. The
 * order is validated while streaming.
 */
#ifndef _READSAM_
#define _READSAM_
#include "global.h"

/* Open the alignment file once. A path ending in ".bam" is streamed through
 * `samtools view`; anything else is read with zlib (plain or gzipped SAM). */
void sam_stream_open(const char *path);

/* Build the accession-name -> (root,node) map a single time, reused by every
 * batch and worker thread. Call after the reference tree has been loaded. */
void build_leaf_map_global(int numberOfTrees);

/* Load the alignment records for the current batch of n reads (query-matrix
 * indices [0,n)) into memory. `paired` selects forward_name vs single name. */
void sam_load_batch(int n, int paired);

/* Release the current batch's records. Call after the batch's worker threads
 * have joined and before the next batch is read. */
void sam_free_batch(void);

/* Close the alignment file and free the global leaf map. */
void sam_stream_close(void);

/*
 * Fill bwa_results[0 .. end-start) for the chunk of reads [start,end) using the
 * records loaded for the current batch, mirroring the contract of main_mem().
 */
int read_sam(int start, int end, bwaMatches *bwa_results, int concordant,
             int numberOfTrees, int paired, int max_query_length,
             int max_readname_length, int max_acc_name);

#endif
