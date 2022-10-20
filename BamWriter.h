//
// Created by 赵展 on 2022/8/17.
//

#ifndef BAMSTATUS_BAMWRITER_H
#define BAMSTATUS_BAMWRITER_Hi
#include "BamTools.h"
#include <thread>
#include "BamCompress.h"
#include "BamWriteCompress.h"
#include "BamCompleteBlock.h"
/*
 *
 * 这个类为对外的输出接口
 *
 * 根据给出的地址，输出BAM文件，并且设置压缩等级
 *
 * 本类使用异步输出，所以需要在析构函数的时候才会完整真正的全部输出
 *
 * 如果需要清空输出，则需要使用over的类内函数进行清空缓存区
 *
 *
 * 需要先输出HDR部分再进行处理
 *
 */

class BamWriter {

public:
    BamWriter(std::string file_name,sam_hdr_t *hdr, int level=6);
    BamWriter(std::string file_name,sam_hdr_t *hdr, int BufferSize, int threadNumber, int level=6);
    ~BamWriter(){
        if (write_block->block_offset>0){
            bam_write_compress->inputUnCompressData(write_block);
        }
        bam_write_compress->WriteComplete();
        for (int i=0;i<n_thread_write;i++) write_compress_thread[i]->join();

    }
    void write(bam1_t* b);

    void over();


private:

    std::thread **write_compress_thread;

    std::thread *write_output_thread;

    BamWriteCompress *bam_write_compress;


    bam_write_block *write_block;

    samFile *output;

    int n_thread_write;



};


#endif //BAMSTATUS_BAMWRITER_H