// gcc -O3 -march=native -pthread  kmer_count_se.c bwt.c bntseq.c bwa.c utils.c ksw.c bwashm.c -lz -lrt -lpthread -o kmer_count_se
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <zlib.h>
#include <pthread.h>

#include "bwa.h"
#include "kseq.h"
#include "utils.h"

KSEQ_INIT(gzFile, gzread)

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
} shared_t;

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

        // 2. 复制必要数据到本地（序列名和序列字符串）
        int len = s->seq->seq.l;
        char *name = strdup(s->seq->name.s);
        char *seq_str = strdup(s->seq->seq.s);
        uint64_t seq_idx = ++s->seq_counter;   // 获取当前序号（从1开始）
        if (!name || !seq_str) {
            perror("strdup");
            exit(1);
        }
        pthread_mutex_unlock(s->read_lock);   // 立即解锁，允许其他线程读取下一条

        // 3. 处理序列（无需再访问共享的 seq）
        uint8_t *q = (uint8_t*)malloc(len);
        if (!q) {
            perror("malloc");
            exit(1);
        }

        int bad = 0;
        int gc_count = 0;
        for (int i = 0; i < len; ++i) {
            q[i] = nst_nt4_table[(int)seq_str[i]];
            if (q[i] > 3) { bad = 1; break; }
            if (q[i] == 1 || q[i] == 2) gc_count++;  // C=1, G=2
        }

        if (bad) {
            local_zero++;
            free(q);
            free(name);
            free(seq_str);
            continue;
        }

        bwtint_t k, l;
        int ret = bwt_match_exact(s->idx->bwt, len, q, &k, &l);

        if (ret > 0) {
            if (ret == 1) {
                // 唯一匹配：计算位置并输出
                bwtint_t sa_pos = bwt_sa(s->idx->bwt, k);
                int is_rev;
                int64_t global_pos = bns_depos(s->idx->bns, sa_pos, &is_rev);
                if (is_rev) global_pos -= len - 1;

                int ref_id = 0;
                while (ref_id < s->idx->bns->n_seqs - 1 &&
                       global_pos >= s->idx->bns->anns[ref_id + 1].offset)
                    ref_id++;

                int64_t chr_pos = global_pos - s->idx->bns->anns[ref_id].offset;
                long long pos_1based = (long long)(chr_pos + 1);

                // 加锁输出
                pthread_mutex_lock(s->print_lock);
                printf("%llu\t%s\t%lld\t%d\t%c\n",
                       (unsigned long long)seq_idx,
                       s->idx->bns->anns[ref_id].name,
                       pos_1based,
                       gc_count,
                       is_rev ? '-' : '+');
                fflush(stdout);  // 确保立即输出
                pthread_mutex_unlock(s->print_lock);

                local_unique++;
            } else {
                local_multi++;
            }
        } else {
            local_zero++;
        }

        free(q);
        free(name);
        free(seq_str);
    }

    // 4. 累加本地计数到全局
    pthread_mutex_lock(s->count_lock);
    s->unique += local_unique;
    s->multi  += local_multi;
    s->zero   += local_zero;
    pthread_mutex_unlock(s->count_lock);

    return NULL;
}

int main(int argc, char *argv[])
{
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <idxbase> <in.fq> <threads>\n", argv[0]);
        return 1;
    }

    const char *idxbase = argv[1];
    const char *fqname  = argv[2];
    int n_threads = atoi(argv[3]);

    // 加载 BWA 索引
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

    // 打开 FASTQ 文件（支持 - 表示 stdin）
    gzFile fp = xzopen(argv[optind + 1], "r");

    if (!fp) {
        fprintf(stderr, "ERROR: Failed to open input '%s'\n", fqname);
        bwa_idx_destroy(idx);
        return 1;
    }

    kseq_t *seq = kseq_init(fp);

    // 初始化互斥锁
    pthread_mutex_t read_lock  = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t count_lock = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t print_lock = PTHREAD_MUTEX_INITIALIZER;

    shared_t shared = {
        .idx = idx,
        .seq = seq,
        .read_lock = &read_lock,
        .unique = 0,
        .multi  = 0,
        .zero   = 0,
        .count_lock = &count_lock,
        .print_lock = &print_lock
    };

    // 创建线程
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

    // 等待所有线程完成
    for (int i = 0; i < n_threads; ++i) {
        pthread_join(threads[i], NULL);
    }

    // 输出统计信息（可选）
    uint64_t total = shared.unique + shared.multi + shared.zero;
    fprintf(stderr, "TOTAL\t%llu\nUNIQUE\t%llu\nMULTI\t%llu\nZERO\t%llu\n",
            (unsigned long long)total,
            (unsigned long long)shared.unique,
            (unsigned long long)shared.multi,
            (unsigned long long)shared.zero);

    // 清理
    free(threads);
    kseq_destroy(seq);
    gzclose(fp);
    bwa_idx_destroy(idx);

    return 0;
}