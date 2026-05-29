#ifndef KMER_COMMON_H
#define KMER_COMMON_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ============================================================
 * AZF BED data structures & functions
 * ============================================================ */

typedef struct {
	char *chr;
	int64_t start;  /* 0-based, inclusive */
	int64_t end;    /* 0-based, exclusive */
} bed_interval_t;

typedef struct {
	bed_interval_t *a;
	int n, m;
} bed_file_t;

static bed_file_t *bed_load(const char *fn)
{
	FILE *fp = fopen(fn, "r");
	if (!fp) {
		fprintf(stderr, "WARNING: failed to open BED file '%s'\n", fn);
		return NULL;
	}

	bed_file_t *bf = calloc(1, sizeof(bed_file_t));
	char line[4096];
	while (fgets(line, sizeof(line), fp)) {
		if (bf->n == bf->m) {
			bf->m = bf->m ? bf->m * 2 : 64;
			bf->a = realloc(bf->a, bf->m * sizeof(bed_interval_t));
		}
		char chr[256];
		int64_t start, end;
		if (sscanf(line, "%255s %ld %ld", chr, &start, &end) == 3) {
			bf->a[bf->n].chr  = strdup(chr);
			bf->a[bf->n].start = start;
			bf->a[bf->n].end   = end;
			bf->n++;
		}
	}
	fclose(fp);
	return bf;
}

static void bed_destroy(bed_file_t *bf)
{
	if (!bf) return;
	for (int i = 0; i < bf->n; i++) free(bf->a[i].chr);
	free(bf->a);
	free(bf);
}

/* returns 1 if (chr, pos) falls within any BED interval, else 0 */
static int bed_query(const bed_file_t *bf, const char *chr, int64_t pos)
{
	if (!bf) return 0;
	for (int i = 0; i < bf->n; i++) {
		if (strcmp(bf->a[i].chr, chr) == 0 &&
		    pos >= bf->a[i].start && pos < bf->a[i].end)
			return 1;
	}
	return 0;
}

/* ============================================================
 * Chain / Liftover data structures & functions
 * ============================================================ */

typedef struct {
	int64_t t_pos;  /* target start of this block (accumulated) */
	int64_t q_pos;  /* query start of this block (accumulated) */
	int64_t size;   /* block length */
} chain_block_t;

typedef struct {
	char *tName;
	int64_t tSize;
	char tStrand;
	int64_t tStart;
	int64_t tEnd;
	char *qName;
	int64_t qSize;
	char qStrand;
	int64_t qStart;
	int64_t qEnd;
	int id;
	int score;

	chain_block_t *blocks;
	int n_blocks;
} chain_t;

typedef struct {
	chain_t *a;
	int n, m;
} chain_file_t;

static chain_file_t *chain_load(const char *fn)
{
	FILE *fp = fopen(fn, "r");
	if (!fp) {
		fprintf(stderr, "WARNING: failed to open chain file '%s'\n", fn);
		return NULL;
	}

	chain_file_t *cf = calloc(1, sizeof(chain_file_t));
	chain_t *cur = NULL;
	int64_t t_current = 0, q_current = 0;
	char line[65536];

	while (fgets(line, sizeof(line), fp)) {
		if (line[0] == '#') continue;
		if (strncmp(line, "chain ", 6) == 0) {
			/* header line */
			if (cf->n == cf->m) {
				cf->m = cf->m ? cf->m * 2 : 256;
				cf->a = realloc(cf->a, cf->m * sizeof(chain_t));
			}
			cur = &cf->a[cf->n++];
			memset(cur, 0, sizeof(chain_t));

			char tName[256], qName[256], tStrand[2], qStrand[2];
			long long tmp_v[6];
			sscanf(line, "chain %d %255s %lld %1s %lld %lld %255s %lld %1s %lld %lld %d",
			       &cur->score, tName, &tmp_v[0], tStrand,
			       &tmp_v[1], &tmp_v[2],
			       qName, &tmp_v[3], qStrand,
			       &tmp_v[4], &tmp_v[5], &cur->id);
			cur->tSize  = tmp_v[0]; cur->tStart = tmp_v[1]; cur->tEnd = tmp_v[2];
			cur->qSize  = tmp_v[3]; cur->qStart = tmp_v[4]; cur->qEnd = tmp_v[5];
			cur->tName = strdup(tName);
			cur->qName = strdup(qName);
			cur->tStrand = tStrand[0];
			cur->qStrand = qStrand[0];
			cur->n_blocks = 0;

			/* first block starts at chain start */
			t_current = cur->tStart;
			q_current = cur->qStart;
		} else if (cur) {
			long long size, dt, dq;
			if (sscanf(line, "%lld %lld %lld", &size, &dt, &dq) == 3) {
				cur->n_blocks++;
				cur->blocks = realloc(cur->blocks, cur->n_blocks * sizeof(chain_block_t));
				chain_block_t *blk = &cur->blocks[cur->n_blocks - 1];

				blk->t_pos = t_current;
				blk->q_pos = q_current;
				blk->size  = size;
				/* dt/dq are the gaps between this block and the NEXT block */
				t_current += size + dt;
				q_current += size + dq;
			} else {
				/* blank line or unparseable → end current chain */
				cur = NULL;
			}
		}
	}
	fclose(fp);
	fprintf(stderr, "Loaded %d chains from %s\n", cf->n, fn);
	return cf;
}

static void chain_destroy(chain_file_t *cf)
{
	if (!cf) return;
	for (int i = 0; i < cf->n; i++) {
		free(cf->a[i].tName);
		free(cf->a[i].qName);
		free(cf->a[i].blocks);
	}
	free(cf->a);
	free(cf);
}

/*
 * Liftover from target (BWA-index genome) to query (output genome).
 * Maps from the chain's target to query direction.
 *
 * Returns 1 on success (sets *q_chr and *q_pos, *q_strand),
 * or 0 if no chain covers this position.
 */
static int chain_liftover(const chain_file_t *cf, const char *t_chr, int64_t t_pos,
                          char **q_chr, int64_t *q_pos, char *q_strand)
{
	if (!cf) return 0;

	/* linear scan; a binary-search-index could be added for very large chain sets */
	for (int i = 0; i < cf->n; i++) {
		chain_t *c = &cf->a[i];

		if (strcmp(c->tName, t_chr) != 0) continue;
		if (t_pos < c->tStart || t_pos >= c->tEnd) continue;

		/* walk blocks and map t_pos → q_pos */
		for (int j = 0; j < c->n_blocks; j++) {
			chain_block_t *blk = &c->blocks[j];
			if (t_pos >= blk->t_pos && t_pos < blk->t_pos + blk->size) {
				int64_t offset = t_pos - blk->t_pos;
				*q_chr = c->qName;
				*q_pos = blk->q_pos + offset;
				if (q_strand) {
					/* For reverse strand chains, we keep the mapping as-is;
					 * both target and query are on the same strand for
					 * most chains, but if query strand is '-',
					 * the position is already in query coordinates. */
					*q_strand = c->qStrand;
				}
				return 1;
			}
		}

		/* position is in a gap between blocks → not alignable */
		return 0;
	}
	return 0;
}

/*
 * Helper: build SA coordinate → (chr, pos, strand) using BWA index.
 * Returns 1 on success.
 */
static int bwa_pos2coord(const bwaidx_t *idx, bwtint_t sa_pos, int read_len,
                         char **chr, long long *pos, char *strand)
{
	int is_rev;
	int64_t global_pos = bns_depos(idx->bns, sa_pos, &is_rev);
	if (is_rev) global_pos -= read_len - 1;

	int ref_id = 0;
	while (ref_id < idx->bns->n_seqs - 1 &&
	       global_pos >= idx->bns->anns[ref_id + 1].offset)
		ref_id++;

	int64_t chr_pos = global_pos - idx->bns->anns[ref_id].offset;
	*chr = idx->bns->anns[ref_id].name;
	*pos = chr_pos + 1; /* 1-based */
	*strand = is_rev ? '-' : '+';
	return 1;
}

#endif
