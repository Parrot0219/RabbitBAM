//
// Created by 赵展 on 2021/3/10.
//
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <zlib.h>
#include <htslib/khash.h>
#include <stdint.h>
#include <chrono>
#include "config.h"
#include "BamBlock.h"
#include "Buffer.h"
#include "BamStatus.h"
#include "Duplicate.h"
#include "Overrepresent.h"
#include "CLI/CLI.hpp"

#include "BamRead.h"
#include "BamCompress.h"
#include "BamCompleteBlock.h"
#define BLOCK_HEADER_LENGTH 18
#define BLOCK_FOOTER_LENGTH 8

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

long long NUM_N[100]={0};
long long NUM_M[100]={0};
long long NUM_TID[100][1000]={0};
void read_pack(BGZF *fp,BamRead *read){
    bam_block * b;
    b=read->getEmpty();
    int count=0;
    while(read_block(fp,b)==0){
        read->inputBlock(b);
//        printf("read block is %d\n",++count);
        b=read->getEmpty();
    }
    read->ReadComplete();
}
void write_pack(Buffer *buffer){
    while(!buffer->is_complete()){
        std::this_thread::sleep_for(chrono::milliseconds(10));
        buffer->output();
    }
}
void compress_pack(BamRead *read,BamCompress *compress){
    pair<bam_block *,int> comp;
    bam_block *un_comp = compress->getEmpty();
    while (1){
        // fg = getRead(comp);
        //printf("%d is not get One compressed data\n",id);
        comp=read->getReadBlock();
        //printf("%d is get One compressed data\n",id);
        if (comp.second<0) {
            //printf("%d is Over\n",id);
            break;
        }
        block_decode_func(comp.first,un_comp);
        read->backBlock(comp.first);

        std::pair<int,int> tmp_pair= find_divide_pos_and_get_read_number(un_comp);
        un_comp->split_pos=tmp_pair.first,un_comp->bam_number=tmp_pair.second;
        compress->inputUnCompressData(un_comp,comp.second);
//        while (!compress->tryinputUnCompressData(un_comp,comp.second)){
//            std::this_thread::sleep_for(std::chrono::milliseconds(1));
//        }
        un_comp = compress->getEmpty();
    }
    compress->CompressThreadComplete();
}
void assign_pack(BamCompress* compress,BamCompleteBlock* completeBlock){
    bam_block *un_comp = nullptr;
    bam_complete_block *assign_block = completeBlock->getEmpty();
    int need_block_len=0,input_length=0;
    int last_use_block_length=0;
    bool isclean = true;
    int ret = -1;
    while (1){
        // fg = getRead(comp);
        //printf("%d is not get One compressed data\n",id);
        if ( isclean && un_comp!=nullptr) {
            compress->backEmpty(un_comp);
        }
//        printf("here?\n");
        un_comp = compress->getUnCompressData();
//        printf("here??\n");
        if (un_comp == nullptr) {
            break;
        }
        /*
         *  放满一整个 bam_complete_block
         */
        printf("here???\n");
//        printf("last use block len is %d\n",last_use_block_length);

        ret = un_comp -> split_pos;
//        ret = find_divide_pos(un_comp);
//        printf("ret number is %d\n",ret);
//        if (ret < 0){
//            printf("unsigned int is wrong\n");
//        }
//            Rabbit_memcpy(&need_block_len,un_comp->data+last_use_block_length,4);
        need_block_len=ret;
//        printf("need block len is %d\n",need_block_len);
//        printf("now_push_length is %d\n",now_push_length);
//        printf("un comp length is %d\n",un_comp->length);
        if (assign_block->length + need_block_len > assign_block->data_size){
            completeBlock->inputCompleteBlock(assign_block);
            assign_block = completeBlock->getEmpty();
        }
        if (ret!=un_comp->length){ // 该分支未经测试
//            printf("Input This\n");

//            printf("un comp length is %d\n",un_comp->length);
            memcpy(assign_block->data+assign_block->length, un_comp->data,ret*sizeof(char));
            assign_block->length+=ret;
            completeBlock->inputCompleteBlock(assign_block);
            assign_block = completeBlock->getEmpty();

            memcpy(assign_block->data+assign_block->length, un_comp->data+ret,(un_comp->length - ret)*sizeof(char));
            assign_block->length += (un_comp->length - ret);
            compress->backEmpty(un_comp);
            un_comp = compress->getUnCompressData();
            memcpy(assign_block->data+assign_block->length, un_comp->data,un_comp->length*sizeof(char));
            assign_block->length += un_comp->length;
            while (find_divide_pos(assign_block) != assign_block->length){
//                printf("find divide pos is %d\n",find_divide_pos(assign_block));
//                printf("assign block length is %d\n",assign_block->length);
                compress->backEmpty(un_comp);
                un_comp = compress->getUnCompressData();
                printf("assign block length is %d\n",assign_block->length);
                printf("assign block data size is %d\n",assign_block->data_size);
                printf("un comp length is %d\n",un_comp->length);
                if (assign_block->length+un_comp->length > assign_block->data_size){
                    change_data_size(assign_block);
                }
                memcpy(assign_block->data+assign_block->length, un_comp->data,un_comp->length*sizeof(char));
                assign_block->length += un_comp->length;
//                printf("assign_block->length is oooooooo %d\n",assign_block->length);
            }
//            printf("OK is Over");
            completeBlock->inputCompleteBlock(assign_block);
            assign_block = completeBlock->getEmpty();
        } else {
//            printf("Input Two?\n");
            if (ret != un_comp->length ) {
//                printf("ai nan ding\n");
                break;
            }
            memcpy(assign_block->data+assign_block->length, un_comp->data,ret*sizeof(char));
            assign_block->length += ret ;
            last_use_block_length = 0;
            isclean = true;
        }
    }
//    printf("here??????????????????????\n");
    if (assign_block->length != 0){
        completeBlock->inputCompleteBlock(assign_block);
    }
    completeBlock->is_over();
    //        if (assign_block->length + un_comp->length >BGZF_MAX_BLOCK_COMPLETE_SIZE){
//            completeBlock->inputCompleteBlock(assign_block);
//            assign_block = completeBlock->getEmpty();
//        }

//        int ret = find_divide_pos(un_comp);
//        if (ret != un_comp->length){
//            printf("ret == %d  block length == %d\n",ret,un_comp->length);
//        }
//        memcpy(assign_block->data+assign_block->length, un_comp->data,un_comp->length*sizeof(char));
//        assign_block->length += un_comp->length;
//        compress->backEmpty(un_comp);

}

void benchmark_pack(BamCompress* compress,BamCompleteBlock* completeBlock){
    bam_block *un_comp = nullptr;
    bam_complete_block *assign_block = completeBlock->getEmpty();
    int need_block_len=0,input_length=0;
    int last_use_block_length=0;
    bool isclean = true;
    int ret = -1;
    int bam_number = 0;
    while (1){
        // fg = getRead(comp);
        //printf("%d is not get One compressed data\n",id);
        if ( isclean && un_comp!=nullptr) {
            compress->backEmpty(un_comp);
        }
//        printf("here?\n");
        un_comp = compress->getUnCompressData();
//        printf("here??\n");
        if (un_comp == nullptr) {
            break;
        }
        /*
         *  放满一整个 bam_complete_block
         */
//        printf("here???\n");
//        printf("last use block len is %d\n",last_use_block_length);

        ret = un_comp -> split_pos;
//        ret = find_divide_pos(un_comp);
//        printf("ret number is %d\n",ret);
//        if (ret < 0){
//            printf("unsigned int is wrong\n");
//        }
//            Rabbit_memcpy(&need_block_len,un_comp->data+last_use_block_length,4);
        need_block_len=ret;
//        printf("need block len is %d\n",need_block_len);
//        printf("now_push_length is %d\n",now_push_length);
//        printf("un comp length is %d\n",un_comp->length);
        if (assign_block->length + need_block_len > assign_block->data_size){
            completeBlock->backEmpty(assign_block);
            assign_block = completeBlock->getEmpty();
        }
        if (ret!=un_comp->length){ // 该分支未经测试
//            printf("Input This\n");

//            printf("un comp length is %d\n",un_comp->length);
            memcpy(assign_block->data+assign_block->length, un_comp->data,ret*sizeof(char));
            assign_block->length+=ret;
            completeBlock->backEmpty(assign_block);
            bam_number+=find_divide_pos_and_get_read_number(assign_block).second;
            assign_block = completeBlock->getEmpty();

            memcpy(assign_block->data+assign_block->length, un_comp->data+ret,(un_comp->length - ret)*sizeof(char));
            assign_block->length += (un_comp->length - ret);
            compress->backEmpty(un_comp);
            un_comp = compress->getUnCompressData();
            memcpy(assign_block->data+assign_block->length, un_comp->data,un_comp->length*sizeof(char));
            assign_block->length += un_comp->length;
            while (find_divide_pos(assign_block) != assign_block->length){
//                printf("find divide pos is %d\n",find_divide_pos(assign_block));
//                printf("assign block length is %d\n",assign_block->length);
                compress->backEmpty(un_comp);
                un_comp = compress->getUnCompressData();
//                printf("assign block length is %d\n",assign_block->length);
//                printf("assign block data size is %d\n",assign_block->data_size);
//                printf("un comp length is %d\n",un_comp->length);
                if (assign_block->length+un_comp->length > assign_block->data_size){
                    change_data_size(assign_block);
                }
                memcpy(assign_block->data+assign_block->length, un_comp->data,un_comp->length*sizeof(char));
                assign_block->length += un_comp->length;
//                printf("assign_block->length is oooooooo %d\n",assign_block->length);
            }

//            printf("OK is Over");
            bam_number+=find_divide_pos_and_get_read_number(assign_block).second;
            completeBlock->backEmpty(assign_block);
            assign_block = completeBlock->getEmpty();
        } else {
//            printf("Input Two?\n");
            if (ret != un_comp->length ) {
                printf("ai nan ding\n");
                break;
            }
//            printf("1\n");
//            printf("assign block length is %d\n",assign_block->length);
//            printf("assign block data size is %d\n",assign_block->data_size);
//            printf("un comp length is %d\n",un_comp->length);
            memcpy(assign_block->data+assign_block->length, un_comp->data,ret*sizeof(unsigned char));

//            printf("2\n");
            assign_block->length += ret ;
            last_use_block_length = 0;
            isclean = true;
            bam_number+=un_comp->bam_number;
        }
    }
//    printf("here??????????????????????\n");
    if (assign_block->length != 0){
        completeBlock->backEmpty(assign_block);
    }
    completeBlock->is_over();
    printf("Bam number is %d\n",bam_number);
    //        if (assign_block->length + un_comp->length >BGZF_MAX_BLOCK_COMPLETE_SIZE){
//            completeBlock->inputCompleteBlock(assign_block);
//            assign_block = completeBlock->getEmpty();
//        }

//        int ret = find_divide_pos(un_comp);
//        if (ret != un_comp->length){
//            printf("ret == %d  block length == %d\n",ret,un_comp->length);
//        }
//        memcpy(assign_block->data+assign_block->length, un_comp->data,un_comp->length*sizeof(char));
//        assign_block->length += un_comp->length;
//        compress->backEmpty(un_comp);

}

void benchmark_bam_pack(BamCompleteBlock* completeBlock){

    bam1_t *b;
    if ((b = bam_init1()) == NULL) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }

    bam_complete_block* un_comp;
    long long ans = 0;
    long long res = 0;
    while (1){
        un_comp = completeBlock->getCompleteBlock();
        if (un_comp == nullptr){
            break;
        }
//        printf("assign over block length is %d\n",un_comp->length);
        int ret;
        while ((ret=(read_bam(un_comp,b,0)))>=0) {
//            printf("One Bam1_t Size is %d\n",ret);
//            printf("This Bam1_t Char Number is %d\n",b->core.l_qseq);
            ans++;
        }
        res++;

        completeBlock->backEmpty(un_comp);
    }
    printf("Bam1_t Number is %lld\n",ans);
    printf("Block  Number is %lld\n",res);
}
int main(int argc,char* argv[]){
    CLI::App app("RabbitBAM");
//    printf("yesyesyesyes\n");
    CLI::App *bam2fq = app.add_subcommand("bam2fq", "BAM format turn to FastQ format");
    CLI::App *bamstatus = app.add_subcommand("bamstatus", "Analyze BAM files");
    CLI::App *benchmark = app.add_subcommand("benchmark", "Performance Testing");
//    printf("yesyesyesyes\n");
    string inputfile;

    string outputfile("./BAMStatus.html");
    int n_thread=1;
    bam2fq->add_option("-i", inputfile, "input File name")->required();
    bam2fq->add_option("-o", outputfile, "output File name");
    bam2fq->add_option("-w,-@,-n,--threads",n_thread,"thread number");

    bamstatus->add_option("-i", inputfile, "input File name")->required();
    bamstatus->add_option("-o", outputfile, "output File name");
    bamstatus->add_option("-w,-@,-n,--threads",n_thread,"thread number");

    benchmark->add_option("-i", inputfile, "input File name")->required();
    benchmark->add_option("-o", outputfile, "output File name");
    benchmark->add_option("-w,-@,-n,--threads",n_thread,"thread number");
//    printf("yesyesyesyes\n");
    CLI11_PARSE(app, argc, argv);
    if (app.get_subcommands().size()>1){
        printf("you should input one command!!!\n");
        return 0;
    }
    if (strcmp(app.get_subcommands()[0]->get_name().c_str(), "benchmark")==0){
        TDEF(fq)
        TSTART(fq)
        printf("Starting Running Benchmark\n");
        printf("BGZF_MAX_BLOCK_COMPLETE_SIZE is %d\n",BGZF_MAX_BLOCK_COMPLETE_SIZE);
        if (strcmp(outputfile.substr(outputfile.size()-4).c_str(),"html")==0) outputfile=("./output.fastq");
        samFile *sin;
        sam_hdr_t *hdr;
        ofstream fout;
        fout.open(outputfile);
        if ((sin=sam_open(inputfile.c_str(),"r"))==NULL){
            printf("Can`t open this file!\n");
            return 0;
        }
        if ((hdr = sam_hdr_read(sin)) == NULL) {
            return  0;
        }
        /*
         *  读取和处理准备
         */
        BamRead read(2000);
        BamCompress compress(2000,n_thread);
        BamCompleteBlock completeBlock(10);

        printf("Malloc Memory is Over\n");
        /*
         * 分析准备
         */
        thread *read_thread = new thread(&read_pack,sin->fp.bgzf,&read);
        thread **compress_thread = new thread *[n_thread];

        for (int i=0;i<n_thread;i++){
            compress_thread[i]=new thread(&compress_pack,&read,&compress);
        }
        thread *assign_thread = new thread(&benchmark_pack,&compress,&completeBlock);
//        thread *consumer_thread = new thread(&benchmark_pack,&completeBlock);
        read_thread->join();
        for (int i=0;i<n_thread;i++) compress_thread[i]->join();
        assign_thread->join();
//        consumer_thread->join();
        sam_close(sin);
        printf("Wait num is %d\n",compress.wait_num);
        TEND(fq)
        TPRINT(fq,"time is : ");
    }

}

//if (strcmp(app.get_subcommands()[0]->get_name().c_str(), "benchmark")==0){
//TDEF(fq)
//TSTART(fq)
//printf("Starting Running Benchmark\n");
//printf("BGZF_MAX_BLOCK_COMPLETE_SIZE is %d\n",BGZF_MAX_BLOCK_COMPLETE_SIZE);
//if (strcmp(outputfile.substr(outputfile.size()-4).c_str(),"html")==0) outputfile=("./output.fastq");
//samFile *sin;
//sam_hdr_t *hdr;
//ofstream fout;
//fout.open(outputfile);
//if ((sin=sam_open(inputfile.c_str(),"r"))==NULL){
//printf("Can`t open this file!\n");
//return 0;
//}
//if ((hdr = sam_hdr_read(sin)) == NULL) {
//return  0;
//}
///*
// *  读取和处理准备
// */
//BamRead read(8000);
//BamCompress compress(4000,n_thread);
//BamCompleteBlock completeBlock(4000);
//
//printf("Malloc Memory is Over\n");
///*
// * 分析准备
// */
//thread *read_thread = new thread(&read_pack,sin->fp.bgzf,&read);
//thread **compress_thread = new thread *[n_thread];
//
//for (int i=0;i<n_thread;i++){
//compress_thread[i]=new thread(&compress_pack,&read,&compress);
//}
//thread *assign_thread = new thread(&assign_pack,&compress,&completeBlock);
////        thread *consumer_thread = new thread(&benchmark_pack,&completeBlock);
//int  consumer_thread_number = 1;
//thread **consumer_thread = new thread*[consumer_thread_number];
//for (int i=0;i<consumer_thread_number;i++) consumer_thread[i] = new thread(&benchmark_bam_pack,&completeBlock);
//read_thread->join();
//for (int i=0;i<n_thread;i++) compress_thread[i]->join();
//assign_thread->join();
////        consumer_thread->join();
//for (int i=0;i<consumer_thread_number;i++) consumer_thread[i]->join();
//sam_close(sin);
//printf("Wait num is %d\n",compress.wait_num);
//TEND(fq)
//TPRINT(fq,"time is : ");
//}

/*
 * check.bam
 * 块数 1338201
 */


