// gcc -O3 -march=native -pthread kmer_count_se.c bwt.c bntseq.c bwa.c utils.c ksw.c bwashm.c -lz -lrt -lpthread -o kmer_count_se
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

#define MAX_MULTI_ENUM 10000

typedef struct {
	bwaidx_t *idx;               // 只读，线程共享
	kseq_t *seq;                  // 共享，需要锁保护
	pthread_mutex_t *read_lock;   // 保护 kseq_read
	uint64_t seq_counter;          // 序列计数器

	uint64_t unique;
	uint64_t multi;
	uint64_t zero;
	pthread_mutex_t *count_lock;  // 保护统计计数器

	pthread_mutex_t *print_lock;  // 保护输出

	bed_file_t *azf;              // AZF BED 区域(可为NULL)
	chain_file_t *chain;          // liftover chain (可为NULL)
} shared_t;

/*
 * 输出一行结果 (对唯一匹配或多匹配AZF区域都适用)
 * 如果 chain 存在，先尝试 liftover；失败则用回原始坐标
 */
static void output_line(shared_t *s, uint64_t seq_idx,
                        char *orig_chr, long long orig_pos, int gc_count, char orig_strand)
{
	char *out_chr = orig_chr;
	long long out_pos = orig_pos;
	char out_strand = orig_strand;

	if (s->chain) {
		char *q_chr = NULL;
		int64_t q_pos = 0;
		char q_strand = 0;
		/* 0-based position for chain_liftover */
		if (chain_liftover(s->chain, orig_chr, orig_pos - 1, &q_chr, &q_pos, &q_strand)) {
			out_chr = q_chr;
			out_pos = q_pos + 1; /* back to 1-based */
			out_strand = q_strand;
		} else {
			/* liftover 失败：输出 0 */
			out_chr = "0";
			out_pos = 0;
			out_strand = '0';
		}
	}

	pthread_mutex_lock(s->print_lock);
	printf("%llu\t%s\t%lld\t%d\t%c\n",
	       (unsigned long long)seq_idx,
	       out_chr,
	       out_pos,
	       gc_count,
	       out_strand);
	fflush(stdout);
	pthread_mutex_unlock(s->print_lock);
}

static void* worker_func(void *arg)
{
	shared_t *s = (shared_t*)arg;

	uint64_t local_unique = 0;
	uint64_t local_multi  = 0;
	uint64_t local_zero   = 0;

	while (1) {
		// 1. 锁定读取一条序列
		pthread_mutex_lock(s->read_lock);
		if (kseq_read(s->seq) < 0) {
			pthread_mutex_unlock(s->read_lock);
			break;
		}

		int len = s->seq->seq.l;
		char *name = strdup(s->seq->name.s);
		char *seq_str = strdup(s->seq->seq.s);
		uint64_t seq_idx = ++s->seq_counter;
		pthread_mutex_unlock(s->read_lock);

		if (!name || !seq_str) {
			perror("strdup");
			exit(1);
		}

		// 2. 编码序列
		uint8_t *q = (uint8_t*)malloc(len);
		if (!q) { perror("malloc"); exit(1); }

		int bad = 0;
		int gc_count = 0;
		for (int i = 0; i < len; ++i) {
			q[i] = nst_nt4_table[(int)seq_str[i]];
			if (q[i] > 3) { bad = 1; break; }
			if (q[i] == 1 || q[i] == 2) gc_count++;
		}

		if (bad) {
			local_zero++;
			free(q); free(name); free(seq_str);
			continue;
		}

		bwtint_t k, l;
		int ret = bwt_match_exact(s->idx->bwt, len, q, &k, &l);

		if (ret == 1) {
			// 唯一匹配
			bwtint_t sa_pos = bwt_sa(s->idx->bwt, k);
			char *chr = NULL;
			long long pos = 0;
			char strand = 0;
			bwa_pos2coord(s->idx, sa_pos, len, &chr, &pos, &strand);
			output_line(s, seq_idx, chr, pos, gc_count, strand);
			local_unique++;

		} else if (ret > 1 && s->azf) {
			// 多匹配 + AZF 模式: 检查是否全部命中都在 AZF 区域
			bwtint_t n_hits = l - k + 1;
			if (n_hits > MAX_MULTI_ENUM) {
				local_multi++;
				free(q); free(name); free(seq_str);
				continue;
			}

			int all_in_azf = 1;
			bwtint_t *sa_positions = malloc(n_hits * sizeof(bwtint_t));
			char   **chrs   = malloc(n_hits * sizeof(char*));
			long long *positions = malloc(n_hits * sizeof(long long));
			char    *strands = n_hits ? malloc(n_hits * sizeof(char)) : NULL;

			for (bwtint_t i = 0; i < n_hits; i++) {
				bwtint_t sa_p = bwt_sa(s->idx->bwt, k + i);
				sa_positions[i] = sa_p;
				char *chr = NULL;
				long long pos = 0;
				char strand = 0;
				bwa_pos2coord(s->idx, sa_p, len, &chr, &pos, &strand);
				chrs[i] = chr;
				positions[i] = pos;
				strands[i] = strand;

				if (!bed_query(s->azf, chr, pos - 1)) {
					all_in_azf = 0;
				}
			}

			if (all_in_azf) {
				for (bwtint_t i = 0; i < n_hits; i++) {
					output_line(s, seq_idx, chrs[i], positions[i], gc_count, strands[i]);
				}
				local_unique++;
			} else {
				local_multi++;
			}

			free(sa_positions);
			free(positions);
			free(strands);
			free(chrs);

		} else if (ret > 1) {
			// 普通多匹配（无 AZF）
			local_multi++;
		} else {
			local_zero++;
		}

		free(q);
		free(name);
		free(seq_str);
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

	/* 解析可选参数 */
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
	if (pos_argc < 4) {
		fprintf(stderr,
		        "Usage: %s [--azf <bed_file>] [--chain <chain_file>] <idxbase> <in.fq> <threads>\n",
		        argv[0]);
		return 1;
	}

	const char *idxbase  = argv[1 + pos_skip];
	const char *fqname   = argv[2 + pos_skip];
	int n_threads = atoi(argv[3 + pos_skip]);

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

	/* 打开 FASTQ */
	gzFile fp = xzopen(fqname, "r");
	if (!fp) {
		fprintf(stderr, "ERROR: Failed to open input '%s'\n", fqname);
		bwa_idx_destroy(idx);
		return 1;
	}

	kseq_t *seq = kseq_init(fp);

	pthread_mutex_t read_lock  = PTHREAD_MUTEX_INITIALIZER;
	pthread_mutex_t count_lock = PTHREAD_MUTEX_INITIALIZER;
	pthread_mutex_t print_lock = PTHREAD_MUTEX_INITIALIZER;

	shared_t shared = {
		.idx = idx,
		.seq = seq,
		.read_lock = &read_lock,
		.seq_counter = 0,
		.unique = 0,
		.multi  = 0,
		.zero   = 0,
		.count_lock = &count_lock,
		.print_lock = &print_lock,
		.azf   = azf,
		.chain = chain,
	};

	/* 创建线程 */
	pthread_t *threads = malloc(n_threads * sizeof(pthread_t));
	if (!threads) {
		perror("malloc");
		return 1;
	}
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
	fprintf(stderr, "TOTAL\t%llu\nUNIQUE\t%llu\nMULTI\t%llu\nZERO\t%llu\n",
	        (unsigned long long)total,
	        (unsigned long long)shared.unique,
	        (unsigned long long)shared.multi,
	        (unsigned long long)shared.zero);

	free(threads);
	kseq_destroy(seq);
	gzclose(fp);
	bwa_idx_destroy(idx);
	bed_destroy(azf);
	chain_destroy(chain);

	return 0;
}
