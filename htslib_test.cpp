//
// Created by 赵展 on 2022/3/3.
//
#include <iostream>
//#include "lockfreequeue/concurrentqueue.h"
#include <thread>
#include <htslib/sam.h>
#include <iostream>
#include <htslib/sam.h>
#include <chrono>
#include <fstream>
#include <htslib/thread_pool.h>
#include "src/BamStatus.h"
#include "src/BamReader.h"
#include <sys/stat.h>
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
bool is_empty{false};

int main(int argc,char* argv[]){

    int n_threads_read = 4 ;
    samFile *sin;
    sam_hdr_t *hdr;


    if ((sin=sam_open("./check.bam", "r"))==NULL){
        printf("Can`t open this file!\n");
        return 0;
    }
    if ((hdr = sam_hdr_read(sin)) == NULL) {
        return  0;
    }
    bam1_t *b;
    if ((b = bam_init1()) == NULL) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }
    TDEF(fq)
    TSTART(fq)
    htsThreadPool p = {NULL, 0};
    p.pool = hts_tpool_init(n_threads_read);
    hts_set_opt(sin,  HTS_OPT_THREAD_POOL, &p);
    int num = 0;
    while(sam_read1(sin, hdr, b)>=0){
        num++;
    }
    printf("Bam Number is %d\n",num);
    sam_close(sin);
    TEND(fq)
    TPRINT(fq,"change time is : ");
}

