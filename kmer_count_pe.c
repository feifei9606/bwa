// gcc -O3 -march=native -pthread -o kmer_count_pe kmer_count_pe.c bwt.c bntseq.c bwa.c utils.c ksw.c bwashm.c -lz -lrt -lpthread
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <zlib.h>
#include <pthread.h>

#include "bwa.h"
#include "kseq.h"
#include "utils.h"
#include "kmer_common.h"

KSEQ_INIT(gzFile, gzread)

#define MAX_MULTI_ENUM    10000
#define MAX_POS_PER_READ    100
#define MAX_COMBINATIONS    500

typedef struct {
	bwaidx_t *idx;                 // 只读，线程共享
	kseq_t *seq1;                   // R1 或 interleaved 输入
	kseq_t *seq2;                   // R2（双文件模式），若为 NULL 表示 interleaved
	pthread_mutex_t *read_lock;     // 保护读取一对序列
	uint64_t pair_counter;          // 配对序号
	int interleaved;                // 1=interleaved, 0=双文件

	uint64_t unique;                // 所有唯一匹配的 read 数
	uint64_t multi;                 // 所有多匹配的 read 数
	uint64_t zero;                  // 所有零匹配或无效的 read 数
	pthread_mutex_t *count_lock;    // 保护统计计数器

	pthread_mutex_t *print_lock;    // 保护输出

	bed_file_t *azf;                // AZF BED 区域(可为NULL)
	chain_file_t *chain;            // liftover chain (可为NULL)
} shared_t;

/*
 * 收集一个 read 的可用比对位置列表
 * ret 来自 bwt_match_exact, k/l 是 SA 区间
 * 填充 hit_* 数组，返回命中数 (0 = 无可用位置)
 */
static int collect_hits(bwaidx_t *idx, int ret, bwtint_t k, bwtint_t l,
                        int len, bed_file_t *azf,
                        char **hit_chrs, long long *hit_positions,
                        char *hit_strands, int max_hits)
{
	if (ret == 1) {
		/* 唯一匹配 */
		bwtint_t sa_pos = bwt_sa(idx->bwt, k);
		char *chr = NULL;
		long long pos = 0;
		char strand = 0;
		bwa_pos2coord(idx, sa_pos, len, &chr, &pos, &strand);
		hit_chrs[0] = chr;
		hit_positions[0] = pos;
		hit_strands[0] = strand;
		return 1;
	}

	if (ret > 1 && azf) {
		/* 多匹配 + AZF 模式: 检查所有命中是否都在 AZF */
		bwtint_t n_hits = l - k + 1;
		if (n_hits > MAX_MULTI_ENUM || n_hits > (bwtint_t)max_hits) return -1;

		bwtint_t *sas = malloc(n_hits * sizeof(bwtint_t));
		char   **c = n_hits ? malloc(n_hits * sizeof(char*)) : NULL;
		long long *p = n_hits ? malloc(n_hits * sizeof(long long)) : NULL;
		char    *s = n_hits ? malloc(n_hits * sizeof(char)) : NULL;

		int all_in_azf = 1;
		for (bwtint_t i = 0; i < n_hits; i++) {
			sas[i] = bwt_sa(idx->bwt, k + i);
			c[i] = NULL;
			char strand = 0;
			long long pos = 0;
			bwa_pos2coord(idx, sas[i], len, &c[i], &pos, &strand);
			p[i] = pos;
			s[i] = strand;
			if (!bed_query(azf, c[i], pos - 1)) {
				all_in_azf = 0;
			}
		}

		if (all_in_azf) {
			for (bwtint_t i = 0; i < n_hits; i++) {
				hit_chrs[i] = c[i];
				hit_positions[i] = p[i];
				hit_strands[i] = s[i];
			}
			free(sas);
			free(p);
			free(s);
			free(c);
			return (int)n_hits;
		}

		/* c[i] 指向 BNS 索引字符串，不可 free */
		free(sas); free(c); free(p); free(s);
		return -1;
	}

	/* 零匹配 或 多匹配无 AZF */
	return ret == 0 ? 0 : -1;
}

/*
 * 对一个命中执行 liftover（如果 chain 存在）
 */
static void liftover_hit(chain_file_t *chain,
                         char *orig_chr, long long orig_pos, char orig_strand,
                         char **out_chr, long long *out_pos, char *out_strand)
{
	if (!chain) {
		*out_chr = orig_chr;
		*out_pos = orig_pos;
		*out_strand = orig_strand;
		return;
	}

	char *q_chr = NULL;
	int64_t q_pos = 0;
	char q_strand = 0;
	if (chain_liftover(chain, orig_chr, orig_pos - 1, &q_chr, &q_pos, &q_strand)) {
		*out_chr = q_chr;
		*out_pos = q_pos + 1;
		*out_strand = q_strand;
	} else {
		*out_chr = "0";
		*out_pos = 0;
		*out_strand = '0';
	}
}

static void* worker_func(void *arg)
{
	shared_t *s = (shared_t*)arg;

	uint64_t local_unique = 0;
	uint64_t local_multi  = 0;
	uint64_t local_zero   = 0;

	while (1) {
		char *name1 = NULL, *name2 = NULL;
		char *seq_str1 = NULL, *seq_str2 = NULL;
		int len1 = 0, len2 = 0;
		uint64_t pair_idx = 0;

		/* 1. 读取一对序列 */
		pthread_mutex_lock(s->read_lock);

		if (s->interleaved) {
			int ret1 = kseq_read(s->seq1);
			if (ret1 < 0) { pthread_mutex_unlock(s->read_lock); break; }
			name1 = strdup(s->seq1->name.s);
			seq_str1 = strdup(s->seq1->seq.s);
			len1 = s->seq1->seq.l;

			int ret2 = kseq_read(s->seq1);
			if (ret2 < 0) {
				free(name1); free(seq_str1);
				pthread_mutex_unlock(s->read_lock);
				break;
			}
			name2 = strdup(s->seq1->name.s);
			seq_str2 = strdup(s->seq1->seq.s);
			len2 = s->seq1->seq.l;
			pair_idx = ++s->pair_counter;
		} else {
			int ret1 = kseq_read(s->seq1);
			int ret2 = kseq_read(s->seq2);
			if (ret1 < 0 || ret2 < 0) { pthread_mutex_unlock(s->read_lock); break; }
			name1 = strdup(s->seq1->name.s);
			seq_str1 = strdup(s->seq1->seq.s);
			len1 = s->seq1->seq.l;
			name2 = strdup(s->seq2->name.s);
			seq_str2 = strdup(s->seq2->seq.s);
			len2 = s->seq2->seq.l;
			pair_idx = ++s->pair_counter;
		}

		if (!name1 || !seq_str1 || !name2 || !seq_str2) {
			perror("strdup"); exit(1);
		}
		pthread_mutex_unlock(s->read_lock);

		/* 2. 编码 read 1 */
		uint8_t *q1 = (uint8_t*)malloc(len1);
		if (!q1) { perror("malloc"); exit(1); }
		int bad1 = 0, gc1 = 0;
		for (int i = 0; i < len1; ++i) {
			q1[i] = nst_nt4_table[(int)seq_str1[i]];
			if (q1[i] > 3) { bad1 = 1; break; }
			if (q1[i] == 1 || q1[i] == 2) gc1++;
		}
		bwtint_t k1, l1;
		int ret_match1 = bad1 ? 0 : bwt_match_exact(s->idx->bwt, len1, q1, &k1, &l1);

		/* 3. 编码 read 2 */
		uint8_t *q2 = (uint8_t*)malloc(len2);
		if (!q2) { perror("malloc"); exit(1); }
		int bad2 = 0, gc2 = 0;
		for (int i = 0; i < len2; ++i) {
			q2[i] = nst_nt4_table[(int)seq_str2[i]];
			if (q2[i] > 3) { bad2 = 1; break; }
			if (q2[i] == 1 || q2[i] == 2) gc2++;
		}
		bwtint_t k2, l2;
		int ret_match2 = bad2 ? 0 : bwt_match_exact(s->idx->bwt, len2, q2, &k2, &l2);

		/* 4. 收集两个 read 的命中位置 */
		char *hit1_chrs[MAX_POS_PER_READ];
		long long hit1_positions[MAX_POS_PER_READ];
		char hit1_strands[MAX_POS_PER_READ];

		char *hit2_chrs[MAX_POS_PER_READ];
		long long hit2_positions[MAX_POS_PER_READ];
		char hit2_strands[MAX_POS_PER_READ];

		int n1 = collect_hits(s->idx, ret_match1, k1, l1, len1, s->azf,
		                      hit1_chrs, hit1_positions, hit1_strands,
		                      MAX_POS_PER_READ);
		int n2 = collect_hits(s->idx, ret_match2, k2, l2, len2, s->azf,
		                      hit2_chrs, hit2_positions, hit2_strands,
		                      MAX_POS_PER_READ);

		if (n1 > 0 || n2 > 0) {
			/* 至少一个 read 有有效命中 */
			local_unique++;

			int c1 = (n1 > 0) ? n1 : 1;  /* 如果无命中，用1个占位行 */
			int c2 = (n2 > 0) ? n2 : 1;

			int total_lines = 0;
			for (int i = 0; i < c1 && total_lines < MAX_COMBINATIONS; i++) {
				for (int j = 0; j < c2 && total_lines < MAX_COMBINATIONS; j++) {
					total_lines++;

					/* Read 1 字段 */
					char *out_chr1 = NULL;
					long long out_pos1 = 0;
					char out_strand1 = 0;
					if (n1 > 0) {
						liftover_hit(s->chain, hit1_chrs[i], hit1_positions[i],
						             hit1_strands[i],
						             &out_chr1, &out_pos1, &out_strand1);
					} else {
						out_chr1 = "0"; out_pos1 = 0; out_strand1 = '0';
					}

					/* Read 2 字段 */
					char *out_chr2 = NULL;
					long long out_pos2 = 0;
					char out_strand2 = 0;
					if (n2 > 0) {
						liftover_hit(s->chain, hit2_chrs[j], hit2_positions[j],
						             hit2_strands[j],
						             &out_chr2, &out_pos2, &out_strand2);
					} else {
						out_chr2 = "0"; out_pos2 = 0; out_strand2 = '0';
					}

					pthread_mutex_lock(s->print_lock);
					printf("%llu\t%s\t%lld\t%d\t%c\t%s\t%lld\t%d\t%c\n",
					       (unsigned long long)pair_idx,
					       out_chr1, out_pos1, (n1 > 0 ? gc1 : 0), out_strand1,
					       out_chr2, out_pos2, (n2 > 0 ? gc2 : 0), out_strand2);
					fflush(stdout);
					pthread_mutex_unlock(s->print_lock);
				}
			}

		} else if (n1 == 0 || n2 == 0) {
			/* 至少一个零匹配 */
			local_zero++;
		} else {
			/* 两者都是多匹配且不在 AZF 内 (n1<0 且 n2<0) */
			local_multi++;
		}

		/* 释放从 collect_hits 分配的内存 (通过 idx->bns->anns 字段无需释放) */
		if (n1 > 0) {
			/* hit1_chrs 指向 bns->anns[].name，无需 free */
		}
		if (n2 > 0) {
			/* 同上 */
		}

		free(q1); free(q2);
		free(name1); free(name2);
		free(seq_str1); free(seq_str2);
	}

	pthread_mutex_lock(s->count_lock);
	s->unique += local_unique;
	s->multi  += local_multi;
	s->zero   += local_zero;
	pthread_mutex_unlock(s->count_lock);

	return NULL;
}

int main(int argc, char *argv[])
{
	const char *azf_file   = NULL;
	const char *chain_file = NULL;
	int pos_skip = 0;

	/* 解析可选参数 --azf <file> / --chain <file> */
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "--azf") == 0 && i + 1 < argc) {
			azf_file = argv[++i];
			pos_skip += 2;
		} else if (strcmp(argv[i], "--chain") == 0 && i + 1 < argc) {
			chain_file = argv[++i];
			pos_skip += 2;
		}
	}

	int pos_argc = argc - pos_skip;
	if (pos_argc < 4 || pos_argc > 5) {
		fprintf(stderr,
		        "Usage: %s [--azf <bed_file>] [--chain <chain_file>] <idxbase> <in.fq> <threads>     (interleaved)\n"
		        "   or: %s [--azf <bed_file>] [--chain <chain_file>] <idxbase> <r1.fq> <r2.fq> <threads> (two-file)\n",
		        argv[0], argv[0]);
		return 1;
	}

	const char *idxbase  = argv[1 + pos_skip];
	const char *fqname1  = argv[2 + pos_skip];
	const char *fqname2  = NULL;
	int n_threads;

	if (pos_argc == 4) {
		/* interleaved: 程序名 + idxbase + fq + threads = 4 */
		n_threads = atoi(argv[3 + pos_skip]);
	} else {
		/* two-file: 程序名 + idxbase + fq1 + fq2 + threads = 5 */
		fqname2   = argv[3 + pos_skip];
		n_threads = atoi(argv[4 + pos_skip]);
	}

	/* 加载 AZF BED */
	bed_file_t *azf = NULL;
	if (azf_file) {
		azf = bed_load(azf_file);
		if (!azf) return 1;
	}

	/* 加载 chain 文件 */
	chain_file_t *chain = NULL;
	if (chain_file) {
		chain = chain_load(chain_file);
		if (!chain) return 1;
	}

	/* 加载 BWA 索引 */
	bwaidx_t *idx = bwa_idx_load_from_shm(idxbase);
	if (!idx) {
		idx = bwa_idx_load(idxbase, BWA_IDX_ALL);
	}
	if (!idx) {
		fprintf(stderr, "ERROR: Failed to load index '%s'\n", idxbase);
		return 1;
	}
	if (!idx->bwt || !idx->bns) {
		fprintf(stderr, "ERROR: Index loaded but BWT or BNS part is missing\n");
		bwa_idx_destroy(idx);
		return 1;
	}

	/* 打开输入文件 */
	gzFile fp1 = xzopen(fqname1, "r");
	if (!fp1) {
		fprintf(stderr, "ERROR: Failed to open '%s'\n", fqname1);
		bwa_idx_destroy(idx);
		return 1;
	}
	gzFile fp2 = NULL;
	if (fqname2) {
		fp2 = xzopen(fqname2, "r");
		if (!fp2) {
			fprintf(stderr, "ERROR: Failed to open '%s'\n", fqname2);
			gzclose(fp1);
			bwa_idx_destroy(idx);
			return 1;
		}
	}

	kseq_t *seq1 = kseq_init(fp1);
	kseq_t *seq2 = fqname2 ? kseq_init(fp2) : NULL;

	pthread_mutex_t read_lock  = PTHREAD_MUTEX_INITIALIZER;
	pthread_mutex_t count_lock = PTHREAD_MUTEX_INITIALIZER;
	pthread_mutex_t print_lock = PTHREAD_MUTEX_INITIALIZER;

	shared_t shared = {
		.idx = idx,
		.seq1 = seq1,
		.seq2 = seq2,
		.read_lock = &read_lock,
		.pair_counter = 0,
		.interleaved = (seq2 == NULL) ? 1 : 0,
		.unique = 0,
		.multi  = 0,
		.zero   = 0,
		.count_lock = &count_lock,
		.print_lock = &print_lock,
		.azf   = azf,
		.chain = chain,
	};

	pthread_t *threads = malloc(n_threads * sizeof(pthread_t));
	if (!threads) { perror("malloc"); return 1; }
	for (int i = 0; i < n_threads; ++i) {
		if (pthread_create(&threads[i], NULL, worker_func, &shared) != 0) {
			perror("pthread_create");
			free(threads);
			return 1;
		}
	}

	for (int i = 0; i < n_threads; ++i) {
		pthread_join(threads[i], NULL);
	}

	uint64_t total = shared.unique + shared.multi + shared.zero;
	fprintf(stderr, "TOTAL_READS\t%llu\n", (unsigned long long)total);
	fprintf(stderr, "UNIQUE\t%llu\n", (unsigned long long)shared.unique);
	fprintf(stderr, "MULTI\t%llu\n",  (unsigned long long)shared.multi);
	fprintf(stderr, "ZERO\t%llu\n",   (unsigned long long)shared.zero);

	free(threads);
	kseq_destroy(seq1);
	if (seq2) kseq_destroy(seq2);
	gzclose(fp1);
	if (fp2) gzclose(fp2);
	bwa_idx_destroy(idx);
	bed_destroy(azf);
	chain_destroy(chain);

	return 0;
}
