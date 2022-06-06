//
// Created by 赵展 on 2021/6/18.
//
#include <iostream>
#include <htslib/sam.h>
#include <chrono>
#include <fstream>
#include <htslib/thread_pool.h>
#include "../concurrentqueue.h"
typedef std::chrono::high_resolution_clock Clock;
#ifdef CYCLING
#define TDEF(x_) static unsigned long long int x_##_t0, x_##_t1;
    #define TSTART(x_) x_##_t0 = __rdtsc();
    #define TEND(x_) x_##_t1 = __rdtsc();
    #define TPRINT(x_, str) printf("%-20s \t%.6f\t M cycles\n", str, (double)(x_##_t1 - x_##_t0)/1e6);
#elif defined TIMING
#define TDEF(x_) chrono::high_resolution_clock::time_point x_##_t0, x_##_t1;
#define TSTART(x_) x_##_t0 = Clock::now();
#define TEND(x_) x_##_t1 = Clock::now();
#define TPRINT(x_, str) printf("%-20s \t%.6f\t sec\n", str, chrono::duration_cast<chrono::microseconds>(x_##_t1 - x_##_t0).count()/1e6);
#else
#define TDEF(x_)
#define TSTART(x_)
#define TEND(x_)
#define TPRINT(x_, str)
#endif
using namespace std;
const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";
int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

const long long Maxn = 1e8;
char Buffer[Maxn+10000];
long long pos=0;
uint8_t Base[16] = {0, 65, 67, 0, 71, 0, 0, 0, 84, 0, 0, 0, 0, 0, 0, 78};
uint8_t BaseRever[16] = {0, 84, 71, 0, 67, 0, 0, 0, 65, 0, 0, 0, 0, 0, 0, 78};

int main(){
    // 时间 ： (fprintf) 7.35s  (stream) 7.03s
    samFile *sin;
    sam_hdr_t *hdr;
    bam1_t *b;
    // /home/old_home/haoz/workspace/data/NC/NC_T_1.sorted.bam
    if ((sin=sam_open("/home/old_home/haoz/workspace/data/NC/NC_T_1.sorted.bam", "r"))==NULL){
        printf("Can`t open this file!\n");
        return 0;
    }
    if ((hdr = sam_hdr_read(sin)) == NULL) {
        return  0;
    }
    if ((b = bam_init1()) == NULL) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }
    htsThreadPool p = {NULL, 0};
    p.pool = hts_tpool_init(8);
    hts_set_opt(sin,  HTS_OPT_THREAD_POOL, &p);
    printf("%s\n","**************************************************************");
    uint8_t* quality;
    uint8_t * seq;
    uint8_t* qname;
    int num=1,ret=0;
    FILE *file=fopen("/Users/zhaozhan/CLionProjects/BAM/my1.txt","w");
    ofstream fout;
    fout.open("./head.fq");
    TDEF(fq)
    TSTART(fq)
    while (sam_read1(sin, hdr, b)>=0) {
        //printf("in this \n");
        if (b->core.flag&2048) continue;
        seq=bam_get_seq(b);
        long long len = strlen(bam_get_qname(b));
        if (pos>Maxn-1)
        {
            //fout.write(Buffer,pos);
            pos=0;
        }
        Buffer[pos++]='@';
        memcpy(Buffer+pos,bam_get_qname(b),len);
        pos+=len;
        Buffer[pos++]='\n';
        if (b->core.flag&16){
            for (int i=0;i<b->core.l_qseq;i++) Buffer[pos+b->core.l_qseq-1-i] = BaseRever[bam_seqi(seq,i)];
            pos+=b->core.l_qseq;
            Buffer[pos++]='\n';Buffer[pos++]='+';Buffer[pos++]='\n';
            quality=bam_get_qual(b);
            for (int i=0;i<b->core.l_qseq;i++){ Buffer[pos+i]=(char)quality[b->core.l_qseq-i-1]+33;}
            pos+=b->core.l_qseq;
        }else{
            for (int i=0;i<b->core.l_qseq;i++) Buffer[pos+i] = Base[bam_seqi(seq,i)];
            pos+=b->core.l_qseq;
            Buffer[pos++]='\n';Buffer[pos++]='+';Buffer[pos++]='\n';
            quality=bam_get_qual(b);
            for (int i=0;i<b->core.l_qseq;i++){ Buffer[pos+i]=(char)quality[i]+33;}
            pos+=b->core.l_qseq;
        }
        Buffer[pos++]='\n';
    }

    //fout.write(Buffer,pos);
    sam_close(sin);
    TEND(fq)
    TPRINT(fq,"change time is : ");
}
