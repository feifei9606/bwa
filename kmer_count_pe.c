// gcc -O3 -march=native -pthread -o kmer_count_pe kmer_count_pe.c bwt.c bntseq.c bwa.c utils.c ksw.c bwashm.c -lz -lrt -lpthread -o kmer_count_pe
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
} shared_t;

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

        // 1. 锁定并读取一对序列
        pthread_mutex_lock(s->read_lock);

        if (s->interleaved) {
            // 从同一个文件读取两条记录
            int ret1 = kseq_read(s->seq1);
            if (ret1 < 0) {
                pthread_mutex_unlock(s->read_lock);
                break;
            }
            name1 = strdup(s->seq1->name.s);
            seq_str1 = strdup(s->seq1->seq.s);
            len1 = s->seq1->seq.l;

            int ret2 = kseq_read(s->seq1);
            if (ret2 < 0) {
                // 第二条缺失，释放已分配内存并终止
                free(name1); free(seq_str1);
                pthread_mutex_unlock(s->read_lock);
                break;
            }
            name2 = strdup(s->seq1->name.s);
            seq_str2 = strdup(s->seq1->seq.s);
            len2 = s->seq1->seq.l;

            pair_idx = ++s->pair_counter;
        } else {
            // 双文件模式：分别从两个文件读取
            int ret1 = kseq_read(s->seq1);
            int ret2 = kseq_read(s->seq2);
            if (ret1 < 0 || ret2 < 0) {
                pthread_mutex_unlock(s->read_lock);
                break;
            }
            name1 = strdup(s->seq1->name.s);
            seq_str1 = strdup(s->seq1->seq.s);
            len1 = s->seq1->seq.l;

            name2 = strdup(s->seq2->name.s);
            seq_str2 = strdup(s->seq2->seq.s);
            len2 = s->seq2->seq.l;

            pair_idx = ++s->pair_counter;
        }

        if (!name1 || !seq_str1 || !name2 || !seq_str2) {
            perror("strdup");
            exit(1);
        }
        pthread_mutex_unlock(s->read_lock);   // 立即解锁

        // 2. 处理 read 1
        uint8_t *q1 = (uint8_t*)malloc(len1);
        if (!q1) { perror("malloc"); exit(1); }
        int bad1 = 0;
        int gc1 = 0;
        for (int i = 0; i < len1; ++i) {
            q1[i] = nst_nt4_table[(int)seq_str1[i]];
            if (q1[i] > 3) { bad1 = 1; break; }
            if (q1[i] == 1 || q1[i] == 2) gc1++;   // C=1, G=2
        }
        bwtint_t k1, l1;
        int ret_match1 = bad1 ? 0 : bwt_match_exact(s->idx->bwt, len1, q1, &k1, &l1);

        // 记录 read1 的比对信息（仅当唯一匹配时）
        int r1_unique = 0;
        char *r1_chr = NULL;
        long long r1_pos = 0;
        int r1_gc = 0;
        char r1_strand = 0;

        if (ret_match1 == 1) {
            r1_unique = 1;
            // 计算位置
            bwtint_t sa_pos1 = bwt_sa(s->idx->bwt, k1);
            int is_rev1;
            int64_t global_pos1 = bns_depos(s->idx->bns, sa_pos1, &is_rev1);
            if (is_rev1) global_pos1 -= len1 - 1;

            int ref_id1 = 0;
            while (ref_id1 < s->idx->bns->n_seqs - 1 &&
                   global_pos1 >= s->idx->bns->anns[ref_id1 + 1].offset)
                ref_id1++;
            int64_t chr_pos1 = global_pos1 - s->idx->bns->anns[ref_id1].offset;
            r1_pos = chr_pos1 + 1;          // 1‑based
            r1_chr = s->idx->bns->anns[ref_id1].name;
            r1_gc = gc1;
            r1_strand = is_rev1 ? '-' : '+';
        }

        // 3. 处理 read 2
        uint8_t *q2 = (uint8_t*)malloc(len2);
        if (!q2) { perror("malloc"); exit(1); }
        int bad2 = 0;
        int gc2 = 0;
        for (int i = 0; i < len2; ++i) {
            q2[i] = nst_nt4_table[(int)seq_str2[i]];
            if (q2[i] > 3) { bad2 = 1; break; }
            if (q2[i] == 1 || q2[i] == 2) gc2++;
        }
        bwtint_t k2, l2;
        int ret_match2 = bad2 ? 0 : bwt_match_exact(s->idx->bwt, len2, q2, &k2, &l2);

        int r2_unique = 0;
        char *r2_chr = NULL;
        long long r2_pos = 0;
        int r2_gc = 0;
        char r2_strand = 0;

        if (ret_match2 == 1) {
            r2_unique = 1;
            bwtint_t sa_pos2 = bwt_sa(s->idx->bwt, k2);
            int is_rev2;
            int64_t global_pos2 = bns_depos(s->idx->bns, sa_pos2, &is_rev2);
            if (is_rev2) global_pos2 -= len2 - 1;

            int ref_id2 = 0;
            while (ref_id2 < s->idx->bns->n_seqs - 1 &&
                   global_pos2 >= s->idx->bns->anns[ref_id2 + 1].offset)
                ref_id2++;
            int64_t chr_pos2 = global_pos2 - s->idx->bns->anns[ref_id2].offset;
            r2_pos = chr_pos2 + 1;
            r2_chr = s->idx->bns->anns[ref_id2].name;
            r2_gc = gc2;
            r2_strand = is_rev2 ? '-' : '+';
        }

        // 4. 输出一行（包含两个 read 的信息，非唯一字段用 0 填充）
        if(ret_match1 == 1 || ret_match2 == 1){
            local_unique++;
            pthread_mutex_lock(s->print_lock);
            printf("%llu\t%s\t%lld\t%d\t%c\t%s\t%lld\t%d\t%c\n",
                   (unsigned long long)pair_idx,
                   r1_unique ? r1_chr : "0",
                   r1_unique ? r1_pos : 0,
                   r1_unique ? r1_gc : 0,
                   r1_unique ? r1_strand : '0',
                   r2_unique ? r2_chr : "0",
                   r2_unique ? r2_pos : 0,
                   r2_unique ? r2_gc : 0,
                   r2_unique ? r2_strand : '0');
            fflush(stdout);
            pthread_mutex_unlock(s->print_lock);
        } else if(ret_match1 == 0 || ret_match2 == 0){
            local_zero++;
        } else{
            local_multi++;
        }
        // 5. 清理
        free(q1); free(q2);
        free(name1); free(name2);
        free(seq_str1); free(seq_str2);
    }

    // 6. 累加本地计数到全局
    pthread_mutex_lock(s->count_lock);
    s->unique += local_unique;
    s->multi  += local_multi;
    s->zero   += local_zero;
    pthread_mutex_unlock(s->count_lock);

    return NULL;
}

int main(int argc, char *argv[])
{
    const char *idxbase;
    const char *fqname1, *fqname2 = NULL;
    int n_threads;

    // 解析参数（兼容之前的参数格式，但忽略插入片段范围）
    if (argc < 4 || argc > 7) {
        fprintf(stderr,
                "Usage: %s <idxbase> <in.fq> <threads>      (interleaved mode)\n"
                "   or: %s <idxbase> <read1.fq> <read2.fq> <threads> (two-file mode)\n"
                "   (use '-' for stdin)\n", argv[0], argv[0]);
        return 1;
    }

    idxbase = argv[1];
    if (argc == 4 || argc == 6) {
        // interleaved mode
        fqname1 = argv[2];
        n_threads = atoi(argv[3]);
    } else {
        // two-file mode
        fqname1 = argv[2];
        fqname2 = argv[3];
        n_threads = atoi(argv[4]);
    }

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

    // 打开输入文件
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

    // 初始化互斥锁
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

    // 输出统计信息到 stderr
    uint64_t total = shared.unique + shared.multi + shared.zero;
    fprintf(stderr, "TOTAL_READS\t%llu\n", (unsigned long long)total);
    fprintf(stderr, "UNIQUE\t%llu\n", (unsigned long long)shared.unique);
    fprintf(stderr, "MULTI\t%llu\n", (unsigned long long)shared.multi);
    fprintf(stderr, "ZERO\t%llu\n",   (unsigned long long)shared.zero);

    // 清理
    free(threads);
    kseq_destroy(seq1);
    if (seq2) kseq_destroy(seq2);
    gzclose(fp1);
    if (fp2) gzclose(fp2);
    bwa_idx_destroy(idx);

    return 0;
}