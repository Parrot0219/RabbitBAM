//
// Created by 赵展 on 2022/7/6.
//

#include "BamWriteCompress.h"
//
// Created by 赵展 on 2022/1/21.
//
#define WriteDeBug 0
BamWriteCompress::BamWriteCompress(int BufferSize,int threadNumber){
    blockNum=0;
    blockInputNum=0;
    isWriteComplete=false;


    compress_bg = 0;
    compress_ed = BufferSize-1;
    compress_size = BufferSize+5;
    compress_data = new bam_write_block*[compress_size];
    for (int i=compress_bg;i<=compress_ed;i++) {
        compress_data[i] = new bam_write_block;
        compress_data[i]->block_length=0;
        compress_data[i]->block_offset=0;
        compress_data[i]->status=0;
        compress_data[i]->block_num=-1;
        compress_data[i]->uncompressed_data=new uint8_t[BGZF_BLOCK_SIZE];
        compress_data[i]->compressed_data=new uint8_t[BGZF_MAX_BLOCK_SIZE];
    }

    blockInputNum=0;
    blockInputPos=0;

    need_compress_bg=1;
    need_compress_ed=0;
    need_compress_size = 2*BufferSize+5;
    need_compress_data = new bam_write_block*[need_compress_size];

    consumer_bg=1;
    consumer_ed=0;
    consumer_size=BufferSize+5;
    consumer_data=new bam_write_block*[consumer_size];
    is_ok = new bool[consumer_size];
    for (int i=0;i<consumer_size;i++) is_ok[i]=false;


    compressThread = threadNumber;

    wait_num=0;
}

/*
 * 获取一个空白的内存块
 * 输入：无
 * 输出：bam_block* : 指向该内存块的指针
 */
bam_write_block* BamWriteCompress::getEmpty(){
    while ((compress_ed+1)%compress_size == compress_bg){
        std::this_thread::sleep_for(std::chrono::nanoseconds(5));
    }
    int num = compress_bg;
    compress_bg=(compress_bg+1)%compress_size;
    return compress_data[num];
}

void BamWriteCompress::inputUnCompressData(bam_write_block* data){

    data->block_num=blockInputNum;
    need_compress_data[blockInputNum%need_compress_size]=data;
    blockInputNum+=1;

}
bam_write_block* BamWriteCompress::getUnCompressData(){

    while(1){
        int num = blockInputNum;

        /*
         * blockInputNum++++++++;
         * blockInputPos++++++++;
         */
       /*
        *
        *  这里必须Strong,不接受伪失败
        *
        */

        while(blockInputPos.compare_exchange_strong(num,num,std::memory_order_relaxed)) {
            std::this_thread::sleep_for(std::chrono::nanoseconds(1));
            if (isWriteComplete &&
                blockInputPos.load(std::memory_order_relaxed) == blockInputNum) {
                return nullptr;
            }
            num=blockInputNum;
        }
        if (num < blockInputNum && blockInputPos.compare_exchange_strong(num,num+1,std::memory_order_relaxed)){
            bam_write_block* res = need_compress_data[num%need_compress_size];
            return res;
        }
    }

}


/*
 * 按照顺序插入输出队列
 * 输入：  bam_block* : 压缩完成的数据
 *        int ：读入的顺序编号
 * 输出：  无
 */

void BamWriteCompress::inputCompressData(bam_write_block* data){
    while (data->block_num != blockNum.load(std::memory_order_acq_rel)) {
        wait_num+=1;
        std::this_thread::sleep_for(std::chrono::nanoseconds(data->block_num-blockNum)/8);
    }
    consumer_data[(consumer_ed + 1) % consumer_size] = data;
    consumer_ed = (consumer_ed + 1) % consumer_size;
    blockNum.store(blockNum.load(std::memory_order_acq_rel)+1,std::memory_order_acq_rel);
}





/*
 * 按照顺序获取解压完成的数据
 * 输入：无
 * 输出：bam_block* 解压完成的数据
 */
bam_write_block* BamWriteCompress::getCompressData(){
    while ((consumer_ed+1)%consumer_size == consumer_bg){
        std::this_thread::sleep_for(std::chrono::nanoseconds(1));
        if (compressThread==0 && (consumer_ed+1)%consumer_size == consumer_bg) return nullptr;
    }
    int num = consumer_bg;
    consumer_bg = (consumer_bg+1)%consumer_size;
    return consumer_data[num];
}
/*
 * 返还使用完毕的内存块
 *
 */
void BamWriteCompress::backEmpty(bam_write_block* data){
//    mtx_compress.lock();
    compress_data[(compress_ed+1)%compress_size] = data;
    compress_ed = (compress_ed+1)%compress_size;
//    mtx_compress.unlock();
}


