/*
 * readsam.c  -- see readsam.h
 *
 * The SAM classification logic in read_sam() (the per-line concordant/
 * discordant resolution) is transcribed verbatim from the step-2 parser in
 * bwa_source_files/fastmap.c so that reading a SAM file produces results
 * identical to the internal bwa-mem path given the same alignments.
 *
 * To keep memory flat regardless of read count, the alignment file is streamed
 * in lockstep with tronko's read batches: sam_load_batch() pulls only the
 * records for the current batch's reads, and sam_free_batch() releases them
 * before the next batch. This requires the file to be in read-file order
 * (unsorted aligner output); the order is validated as we stream.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/types.h>
#include <zlib.h>
#include "global.h"
#include "hashmap.h"
#include "hashmap_base.h"
#include "readsam.h"

/* One read's worth of SAM text plus its position within the current batch. */
typedef struct samRec {
	char *lines;     /* all alignment lines for this read, '\n'-separated */
	size_t len;
	size_t cap;
	int order;       /* index of this read within the current batch */
} samRec;

/* Records for the current batch only, keyed by read name. The key pointers are
 * borrowed from the query matrices and are valid for the batch's lifetime. */
static HASHMAP(char, samRec) batchMap;
static int batch_inited = 0;

/* Accession-name -> (root,node) map, built once and reused by every batch. */
static HASHMAP(char, leafMap) g_leafMap;
static int leaf_inited = 0;

/* Persistent alignment-file stream. */
static gzFile g_gz = Z_NULL;   /* SAM (plain or gzipped) */
static FILE  *g_fp = NULL;     /* `samtools view` pipe for BAM */
static int    g_is_bam = 0;

/* Reusable line buffer and a one-record lookahead ("pushback") that belongs to
 * a later batch than the one currently being filled. */
static char  *g_line = NULL;
static size_t g_line_cap = 0;
static char  *g_pb = NULL;     /* pushed-back raw line, or NULL */
static size_t g_pb_cap = 0;
static int    g_pb_valid = 0;
static long   g_total_records = 0;

static void append_line(samRec *rec, const char *line, size_t llen)
{
	int need_nl = (llen == 0 || line[llen - 1] != '\n');
	size_t add = llen + (need_nl ? 1 : 0);
	if (rec->len + add + 1 > rec->cap) {
		while (rec->len + add + 1 > rec->cap)
			rec->cap = rec->cap ? rec->cap * 2 : 512;
		rec->lines = realloc(rec->lines, rec->cap);
	}
	memcpy(rec->lines + rec->len, line, llen);
	rec->len += llen;
	if (need_nl)
		rec->lines[rec->len++] = '\n';
	rec->lines[rec->len] = '\0';
}

void build_leaf_map_global(int numberOfTrees)
{
	int i, j;
	hashmap_init(&g_leafMap, hashmap_hash_string, strcmp);
	leaf_inited = 1;
	for (i = 0; i < numberOfTrees; i++) {
		for (j = numspecArr[i] - 1; j < 2 * numspecArr[i] - 1; j++) {
			struct leafMap *lm = malloc(sizeof(*lm));
			lm->name = treeArr[i][j].name;
			lm->root = i;
			lm->node = j;
			hashmap_put(&g_leafMap, lm->name, lm);
		}
	}
}

/* Read one line (including trailing newline) from the open stream into g_line.
 * Returns the length, or -1 at end of file. */
static long stream_getline(void)
{
	if (g_line_cap == 0) {
		g_line_cap = 65536;
		g_line = (char *)malloc(g_line_cap);
	}
	if (g_is_bam) {
		size_t n = g_line_cap;
		ssize_t got = getline(&g_line, &n, g_fp);
		g_line_cap = n;
		return got; /* -1 at EOF */
	}
	size_t len = 0;
	for (;;) {
		if (len + 1 >= g_line_cap) {
			g_line_cap *= 2;
			g_line = (char *)realloc(g_line, g_line_cap);
		}
		if (gzgets(g_gz, g_line + len, (int)(g_line_cap - len)) == NULL)
			break;
		size_t got = strlen(g_line + len);
		len += got;
		if (got == 0)
			break;
		if (g_line[len - 1] == '\n')
			break;
	}
	return len == 0 ? -1 : (long)len;
}

void sam_stream_open(const char *path)
{
	size_t plen = strlen(path);
	g_is_bam = (plen >= 4 && strcasecmp(path + plen - 4, ".bam") == 0);
	if (g_is_bam) {
		size_t clen = strlen(path) + 64;
		char *cmd = (char *)malloc(clen);
		snprintf(cmd, clen, "samtools view \"%s\"", path);
		g_fp = popen(cmd, "r");
		free(cmd);
		if (g_fp == NULL) {
			fprintf(stderr, "Error: could not run samtools to read BAM file '%s'. Is samtools on your PATH?\n", path);
			exit(1);
		}
	} else {
		g_gz = gzopen(path, "r");
		if (g_gz == Z_NULL) {
			fprintf(stderr, "Error: could not open SAM file '%s'.\n", path);
			exit(1);
		}
	}
}

/* Extract the QNAME (first tab-delimited field) into qbuf. Returns 0 on a
 * header/blank line that should be skipped, 1 otherwise. */
static int parse_qname(const char *line, char *qbuf, size_t qbuf_sz)
{
	if (line[0] == '@' || line[0] == '\n' || line[0] == '\0')
		return 0;
	const char *tab = strchr(line, '\t');
	if (tab == NULL)
		return 0;
	size_t qlen = (size_t)(tab - line);
	if (qlen >= qbuf_sz)
		qlen = qbuf_sz - 1;
	memcpy(qbuf, line, qlen);
	qbuf[qlen] = '\0';
	return 1;
}

/* Try to attach a raw line to the current batch. Returns 1 if the line's read
 * belongs to this batch (and was attached), 0 if it belongs to a later batch. */
static int attach_line(const char *line, size_t llen, int *last_order)
{
	char qbuf[MAXREADNAME + 1];
	if (!parse_qname(line, qbuf, sizeof(qbuf)))
		return 1; /* header/blank: consumed, not for a later batch */
	samRec *rec = hashmap_get(&batchMap, qbuf);
	if (rec == NULL)
		return 0; /* belongs to a later batch */
	if (rec->order < *last_order) {
		fprintf(stderr, "Error: alignment file is not in read order (read '%s' appears out of order). tronko-assign -b requires the aligner's unsorted output, not a coordinate-/name-sorted BAM. Exiting...\n", qbuf);
		exit(1);
	}
	*last_order = rec->order;
	append_line(rec, line, llen);
	g_total_records++;
	return 1;
}

void sam_load_batch(int n, int paired)
{
	int i;
	if (batch_inited)
		sam_free_batch();
	hashmap_init(&batchMap, hashmap_hash_string, strcmp);
	batch_inited = 1;

	/* Seed the map with this batch's read names so membership can be tested. */
	for (i = 0; i < n; i++) {
		char *name = paired ? pairedQueryMat->forward_name[i]
		                    : singleQueryMat->name[i];
		if (hashmap_get(&batchMap, name) != NULL)
			continue; /* duplicate name within batch: keep first */
		samRec *rec = (samRec *)malloc(sizeof(samRec));
		rec->lines = NULL;
		rec->len = 0;
		rec->cap = 0;
		rec->order = i;
		hashmap_put(&batchMap, name, rec);
	}

	int last_order = -1;
	/* A record buffered from the previous batch may belong to this one. */
	if (g_pb_valid) {
		if (attach_line(g_pb, strlen(g_pb), &last_order))
			g_pb_valid = 0; /* consumed */
		else
			return; /* still belongs to a later batch: this batch has none */
	}
	long got;
	while ((got = stream_getline()) != -1) {
		if (attach_line(g_line, (size_t)got, &last_order))
			continue;
		/* Belongs to a later batch: stash it and stop. */
		if ((size_t)got + 1 > g_pb_cap) {
			g_pb_cap = (size_t)got + 1;
			g_pb = (char *)realloc(g_pb, g_pb_cap);
		}
		memcpy(g_pb, g_line, (size_t)got);
		g_pb[got] = '\0';
		g_pb_valid = 1;
		break;
	}
}

void sam_free_batch(void)
{
	if (!batch_inited)
		return;
	samRec *rec;
	hashmap_foreach_data(rec, &batchMap) {
		free(rec->lines);
		free(rec);
	}
	hashmap_cleanup(&batchMap);
	batch_inited = 0;
}

void sam_stream_close(void)
{
	if (g_is_bam) {
		if (g_fp) { pclose(g_fp); g_fp = NULL; }
	} else {
		if (g_gz != Z_NULL) { gzclose(g_gz); g_gz = Z_NULL; }
	}
	free(g_line); g_line = NULL; g_line_cap = 0;
	free(g_pb); g_pb = NULL; g_pb_cap = 0; g_pb_valid = 0;
	if (leaf_inited) {
		struct leafMap *blob;
		hashmap_foreach_data(blob, &g_leafMap) { free(blob); }
		hashmap_cleanup(&g_leafMap);
		leaf_inited = 0;
	}
	fprintf(stderr, "[readsam] Processed %ld alignment records total.\n", g_total_records);
}

int read_sam(int start, int end, bwaMatches *bwa_results, int concordant,
             int numberOfTrees, int paired, int max_query_length,
             int max_readname_length, int max_acc_name)
{
	int i, j, k, l;
	int no_add = 0;
	int success = 1;

	/* Reassemble this chunk's SAM text in read-file order from the batch's
	 * records. Reads with no record contribute nothing (they stay unmapped via
	 * the caller's -1 initialization). For single-end reads we synthesize a
	 * minimal unmapped record so the parser's read index stays aligned without
	 * losing the following read's first line. */
	char *chunk = NULL;
	size_t chunk_len = 0, chunk_cap = 0;
	for (i = start; i < end; i++) {
		const char *name = paired ? pairedQueryMat->forward_name[i]
		                          : singleQueryMat->name[i];
		samRec *rec = hashmap_get(&batchMap, (char *)name);
		if (rec != NULL && rec->len == 0)
			rec = NULL; /* seeded but no alignment lines: treat as missing */
		const char *add = NULL;
		size_t addlen = 0;
		char synth[MAXREADNAME + 32];
		if (rec != NULL) {
			add = rec->lines;
			addlen = rec->len;
		} else if (!paired) {
			snprintf(synth, sizeof(synth), "%s\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n", name);
			add = synth;
			addlen = strlen(synth);
		}
		if (add == NULL)
			continue;
		if (chunk_len + addlen + 1 > chunk_cap) {
			while (chunk_len + addlen + 1 > chunk_cap)
				chunk_cap = chunk_cap ? chunk_cap * 2 : 65536;
			chunk = (char *)realloc(chunk, chunk_cap);
		}
		memcpy(chunk + chunk_len, add, addlen);
		chunk_len += addlen;
		chunk[chunk_len] = '\0';
	}

	if (chunk != NULL && chunk_len > 0) {
		bwaMatches *results = bwa_results; /* alias to match transcribed code */
		char readname[max_readname_length];
		char read1[max_acc_name];
		char read2[max_acc_name];
		char cigar[MAX_CIGAR];
		int start_position;
		char *rest = chunk;
		int decimal = 0;
		char *token;
		j = 0;
		/* ---- begin verbatim transcription of fastmap.c process() step 2 ---- */
		while (token = strtok_r(rest, "\n", &rest)) {
			success = sscanf(token, "%s %d %s %d %*d %s %s %*d %*d %*s %*s %*s %*s %*s %*s %*[^.,;]", readname, &decimal, read1, &start_position, cigar, read2);
			int whichMat = 0;
			if (paired == 1 && j > 0) {
				if (strcmp(readname, pairedQueryMat->forward_name[start + j - 1]) == 0) {
					whichMat = 1;
				}
			} else if (j > 0) {
				if (strcmp(readname, singleQueryMat->name[start + j - 1]) == 0) {
					whichMat = 1;
				}
			}
			if (j > 0 && whichMat == 1) {
				for (k = 0; k < MAX_NUM_BWA_MATCHES; k++) {
					if (results[j - 1].concordant_matches_roots[k] == -1) {
						break;
					}
				}
				for (l = 0; l < k; l++) {
					if (strcmp(read2, "=") == 0) {
						struct leafMap *leaf_map;
						leaf_map = hashmap_get(&g_leafMap, read1);
						if (results[j - 1].concordant_matches_roots[l] == leaf_map->root && results[j - 1].concordant_matches_nodes[l] == leaf_map->node) {
							no_add = 1;
						}
					}
				}
				if (concordant == 1 && no_add == 0 && strcmp(read2, "=") != 0) { no_add = 1; }
				if (no_add == 0 && strcmp(read2, "=") == 0) {
					struct leafMap *leaf_map;
					leaf_map = hashmap_get(&g_leafMap, read1);
					results[j - 1].concordant_matches_roots[k] = leaf_map->root;
					results[j - 1].concordant_matches_nodes[k] = leaf_map->node;
					if (results[j - 1].use_portion == 1) {
						strcpy(results[j - 1].cigars_forward[k], cigar);
						results[j - 1].starts_forward[k] = start_position;
					}
					k++;
				}
				no_add = 0;
				for (l = 0; l < k; l++) {
					if (strcmp(read2, "=") == 0 || strcmp(read2, "*") == 0) {
						no_add = 1;
						break;
					}
					struct leafMap *leaf_map;
					leaf_map = hashmap_get(&g_leafMap, read2);
					if (results[j - 1].concordant_matches_roots[l] == leaf_map->root && results[j - 1].concordant_matches_nodes[l] == leaf_map->node) {
						no_add = 1;
					}
				}
				if (k == 0) { no_add = 1; }
				if (no_add == 0 && strcmp(read2, "*") != 0 && strcmp(read1, "=") != 0) {
					struct leafMap *leaf_map;
					leaf_map = hashmap_get(&g_leafMap, read2);
					results[j - 1].concordant_matches_roots[k] = leaf_map->root;
					results[j - 1].concordant_matches_nodes[k] = leaf_map->node;
					if (results[j - 1].use_portion == 1) {
						strcpy(results[j - 1].cigars_reverse[k - 1], cigar);
						results[j - 1].starts_reverse[k - 1] = start_position;
					}
				}
				no_add = 0;
				for (k = 0; k < MAX_NUM_BWA_MATCHES; k++) {
					if (results[j - 1].discordant_matches_roots[k] == -1) {
						break;
					}
				}
				for (l = 0; l < k; l++) {
					if (strcmp(read2, "=") != 0) {
						struct leafMap *leaf_map;
						leaf_map = hashmap_get(&g_leafMap, read1);
						if (results[j - 1].discordant_matches_roots[l] == leaf_map->root && results[j - 1].discordant_matches_nodes[l] == leaf_map->node) {
							no_add = 1;
						}
					}
				}
				if (no_add == 0 && strcmp(read2, "=") != 0) {
					struct leafMap *leaf_map;
					leaf_map = hashmap_get(&g_leafMap, read1);
					results[j - 1].discordant_matches_roots[k] = leaf_map->root;
					results[j - 1].discordant_matches_nodes[k] = leaf_map->node;
					if (results[j - 1].use_portion == 1 && decimal == 1) {
						strcpy(results[j - 1].cigars_forward[k], cigar);
						results[j - 1].starts_forward[k] = start_position;
					}
					k++;
				}
				no_add = 0;
				for (l = 0; l < k; l++) {
					if (strcmp(read2, "*") == 0 || strcmp(read2, "=") == 0) {
						no_add = 1;
						break;
					}
					struct leafMap *leaf_map;
					leaf_map = hashmap_get(&g_leafMap, read2);
					if (results[j - 1].discordant_matches_roots[l] == leaf_map->root && results[j - 1].discordant_matches_nodes[l] == leaf_map->node) {
						no_add = 1;
					}
				}
				if (no_add == 0 && strcmp(read2, "=") != 0 && strcmp(read2, "*") != 0) {
					struct leafMap *leaf_map;
					leaf_map = hashmap_get(&g_leafMap, read2);
					results[j - 1].discordant_matches_roots[k] = leaf_map->root;
					results[j - 1].discordant_matches_nodes[k] = leaf_map->node;
					if (results[j - 1].use_portion == 1) {
						strcpy(results[j - 1].cigars_reverse[k], cigar);
						results[j - 1].starts_reverse[k] = start_position;
					}
				}
				no_add = 0;
				if (decimal == 0 && results[j - 1].use_portion == 1) {
					strcpy(results[j - 1].cigars_reverse[k], cigar);
					results[j - 1].starts_reverse[k] = start_position;
				}
			} else if (j == 0 && strcmp(read2, "=") == 0) {
				struct leafMap *leaf_map;
				leaf_map = hashmap_get(&g_leafMap, read1);
				results[j].concordant_matches_roots[0] = leaf_map->root;
				results[j].concordant_matches_nodes[0] = leaf_map->node;
				if (results[j].use_portion == 1) {
					strcpy(results[j].cigars_forward[0], cigar);
					results[j].starts_forward[0] = start_position;
				}
				if (strcmp(read2, "=") != 0 && strcmp(read2, "*") != 0 && (results[j].concordant_matches_roots[0] != leaf_map->root && results[j].concordant_matches_nodes[0] != leaf_map->node)) {
					struct leafMap *leaf_map;
					leaf_map = hashmap_get(&g_leafMap, read2);
					results[j].concordant_matches_roots[1] = leaf_map->root;
					results[j].concordant_matches_nodes[1] = leaf_map->node;
					if (results[j].use_portion == 1) {
						strcpy(results[j].cigars_reverse[0], cigar);
						results[j].starts_reverse[0] = start_position;
					}
				}
				j++;
			} else if (j == 0 && strcmp(read2, "=") != 0 && paired != 0 && strcmp(read1, "*") != 0) {
				struct leafMap *leaf_map;
				leaf_map = hashmap_get(&g_leafMap, read1);
				results[j].discordant_matches_roots[0] = leaf_map->root;
				results[j].discordant_matches_nodes[0] = leaf_map->node;
				if (results[j].use_portion == 1) {
					strcpy(results[j].cigars_forward[0], cigar);
					results[j].starts_forward[0] = start_position;
				}
				if (strcmp(read2, "*") != 0 && (results[j].discordant_matches_roots[0] != leaf_map->root && results[j].discordant_matches_nodes[0] != leaf_map->node)) {
					struct leafMap *leaf_map;
					leaf_map = hashmap_get(&g_leafMap, read2);
					results[j].discordant_matches_roots[1] = leaf_map->root;
					results[j].discordant_matches_nodes[1] = leaf_map->node;
					if (results[j].use_portion == 1) {
						strcpy(results[j].cigars_reverse[0], cigar);
						results[j].starts_reverse[0] = start_position;
					}
				}
				j++;
			} else {
				if (paired != 0 && strcmp(pairedQueryMat->forward_name[start + j], readname) == 0 && strcmp(read2, "=") == 0) {
					if (decimal == 1 || decimal == 2) {
						struct leafMap *leaf_map;
						leaf_map = hashmap_get(&g_leafMap, read1);
						results[j].concordant_matches_roots[0] = leaf_map->root;
						results[j].concordant_matches_nodes[0] = leaf_map->node;
						if (results[j].use_portion == 1) {
							strcpy(results[j].cigars_forward[0], cigar);
							results[j].starts_forward[0] = start_position;
						}
					}
					if (strcmp(read2, "=") != 0 && strcmp(read2, "*") != 0 && strcmp(read1, read2) != 0) {
						struct leafMap *leaf_map;
						leaf_map = hashmap_get(&g_leafMap, read2);
						results[j].concordant_matches_roots[1] = leaf_map->root;
						results[j].concordant_matches_nodes[1] = leaf_map->node;
						if (results[j].use_portion == 1) {
							strcpy(results[j].cigars_reverse[0], cigar);
							results[j].starts_reverse[0] = start_position;
						}
					}
					j++;
				} else if (paired != 0 && strcmp(pairedQueryMat->forward_name[start + j], readname) == 0 && strcmp(read2, "=") != 0 && strcmp(read1, "*") != 0) {
					struct leafMap *leaf_map;
					leaf_map = hashmap_get(&g_leafMap, read1);
					results[j].discordant_matches_roots[0] = leaf_map->root;
					results[j].discordant_matches_nodes[0] = leaf_map->node;
					if (strcmp(read2, "*") != 0 && strcmp(read1, read2) != 0) {
						struct leafMap *leaf_map;
						leaf_map = hashmap_get(&g_leafMap, read2);
						results[j].discordant_matches_roots[1] = leaf_map->root;
						results[j].discordant_matches_nodes[1] = leaf_map->node;
						if (results[j].use_portion == 1) {
							strcpy(results[j].cigars_reverse[0], cigar);
							results[j].starts_reverse[0] = start_position;
						}
					}
					j++;
				} else if (paired == 0 && strcmp(singleQueryMat->name[start + j], readname) == 0 && strcmp(read2, "=") == 0) {
					struct leafMap *leaf_map;
					leaf_map = hashmap_get(&g_leafMap, read1);
					results[j].concordant_matches_roots[0] = leaf_map->root;
					results[j].concordant_matches_nodes[0] = leaf_map->node;
					if (results[j].use_portion == 1) {
						strcpy(results[j].cigars_forward[0], cigar);
						results[j].starts_forward[0] = start_position;
					}
					if (strcmp(read2, "=") != 0 && strcmp(read2, "*") != 0 && strcmp(read1, read2) != 0) {
						struct leafMap *leaf_map;
						leaf_map = hashmap_get(&g_leafMap, read2);
						results[j].concordant_matches_roots[1] = leaf_map->root;
						results[j].concordant_matches_nodes[1] = leaf_map->node;
						if (results[j].use_portion == 1) {
							strcpy(results[j].cigars_reverse[0], cigar);
							results[j].starts_reverse[0] = start_position;
						}
					}
					j++;
				} else if (paired == 0 && strcmp(singleQueryMat->name[start + j], readname) == 0 && strcmp(read2, "=") != 0) {
					if (strcmp(read1, "*") != 0) {
						struct leafMap *leaf_map;
						leaf_map = hashmap_get(&g_leafMap, read1);
						results[j].discordant_matches_roots[0] = leaf_map->root;
						results[j].discordant_matches_nodes[0] = leaf_map->node;
						if (results[j].use_portion == 1) {
							strcpy(results[j].cigars_forward[0], cigar);
							results[j].starts_forward[0] = start_position;
						}
						if (strcmp(read2, "*") != 0 && strcmp(read1, read2) != 0) {
							struct leafMap *leaf_map;
							leaf_map = hashmap_get(&g_leafMap, read2);
							results[j].discordant_matches_roots[1] = leaf_map->root;
							results[j].discordant_matches_nodes[1] = leaf_map->node;
							if (results[j].use_portion == 1) {
								strcpy(results[j].cigars_reverse[0], cigar);
								results[j].starts_reverse[0] = start_position;
							}
						}
					} else {
						results[j].concordant_matches_roots[0] = -1;
						results[j].concordant_matches_nodes[0] = -1;
						results[j].discordant_matches_nodes[0] = -1;
						results[j].discordant_matches_roots[0] = -1;
					}
					j++;
				} else {
					if (paired != 0) {
						while (strcmp(pairedQueryMat->forward_name[start + j], readname) != 0) {
							results[j].concordant_matches_roots[0] = -1;
							results[j].concordant_matches_nodes[0] = -1;
							results[j].discordant_matches_nodes[0] = -1;
							results[j].discordant_matches_roots[0] = -1;
							j++;
						}
					} else {
						while (strcmp(singleQueryMat->name[start + j], readname) != 0) {
							results[j].concordant_matches_roots[0] = -1;
							results[j].concordant_matches_nodes[0] = -1;
							results[j].discordant_matches_nodes[0] = -1;
							results[j].discordant_matches_roots[0] = -1;
							j++;
						}
					}
				}
			}
		}
		/* ---- end verbatim transcription ---- */
	}

	free(chunk);
	return 0;
}
