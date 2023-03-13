//
// Created by 赵展 on 2023/3/10.
//
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <stdint.h>
#include <chrono>
#include "config.h"
#include "BamBlock.h"
#include "Buffer.h"
#include "BamStatus.h"
#include "Duplicate.h"
#include "Overrepresent.h"
#include "CLI/CLI.hpp"
#include <sched.h>
#include <unistd.h>
#include "BamRead.h"
#include "BamCompress.h"
#include "BamCompleteBlock.h"
#include "BamTools.h"
#include "BamWriteCompress.h"
#include "BamWriter.h"
#include "BamReader.h"
#define BLOCK_HEADER_LENGTH 18
#define BLOCK_FOOTER_LENGTH 8


#define DEBUG 0


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
//const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";
//int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };
//uint8_t Base[16] = {0, 65, 67, 0, 71, 0, 0, 0, 84, 0, 0, 0, 0, 0, 0, 78};
//uint8_t BaseRever[16] = {0, 84, 71, 0, 67, 0, 0, 0, 65, 0, 0, 0, 0, 0, 0, 78};

int main(int argc,char* argv[]){

    CLI::App app("RabbitBAM");

    string inputfile;

    string outputfile("./BAMStatus.html");
    int n_thread=1;
    int n_thread_write=1;
    int level = 6;


    CLI::App *api_test = app.add_subcommand("api_test","use api to read and write");
    api_test->add_option("-i",inputfile,"input File name")->required()->check(CLI::ExistingFile);
    api_test->add_option("-o",outputfile,"output File name");
    api_test->add_option("--nr",n_thread,"Read thread number");
    api_test->add_option("--nw",n_thread_write,"Write thread number");
    api_test->add_option("-l,--level",level,"zip level");



    CLI11_PARSE(app, argc, argv);
    if (app.get_subcommands().size()>1){
        printf("you should input one command!!!\n");
        return 0;
    }
    if (strcmp(app.get_subcommands()[0]->get_name().c_str(), "api_test")==0){
        TDEF(fq)
        TSTART(fq)

        printf("Starting Running API Test\n");
        printf("BGZF_MAX_BLOCK_COMPLETE_SIZE is %d\n",BGZF_MAX_BLOCK_COMPLETE_SIZE);
        printf("output File Name is %s\n",outputfile.c_str());
        //if (strcmp(outputfile.substr(outputfile.size()-4).c_str(),"html")==0) outputfile=("./output.fastq");
//        samFile *sin;
//        sam_hdr_t *hdr;
//        samFile *output;
//        if ((sin=sam_open(inputfile.c_str(),"r"))==NULL){
//            printf("Can`t open this file!\n");
//            return 0;
//        }
//
//        if ((output=sam_open(outputfile.c_str(),"wb"))==NULL){
//            printf("Can`t open this file!\n");
//            return 0;
//        }
//        if ((hdr = sam_hdr_read(sin)) == NULL) {
//            return  0;
//        }


        /*
         * 开始创建BamRead 和 BamWriter
         */

        BamReader *reader = new BamReader(inputfile,n_thread);
        BamWriter *writer = new BamWriter(outputfile,reader->getHeader(),n_thread_write,level,200);
        bam1_t *b;
        if ((b = bam_init1()) == NULL) {
            fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
        }
        long long num=0;
        while (reader->getBam1_t(b)){
            num++;
            writer->write(b);
//            if (num%1000 == 0) printf("Bam1_t Num is %d\n",num);
        }
        writer->over();

        cout << "Bam1_t Num : "<< num << endl;


        TEND(fq)
        TPRINT(fq,"time is : ");
        /*
         *  二代数据线程读写比例为 1 ：4
         *  三代数据线程读写比例为 1 ：4
         *
         *  测试下来读写速率为700mb/s 在64线程下
         *
         *  注意需要有额外的四个线程
         */
    }
}



