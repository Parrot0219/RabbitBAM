//
// Created by 赵展 on 2022/1/21.
//

#include "BamCompress.h"
BamCompress::BamCompress(int BufferSize,int threadNumber){
    blockNum=0;

    compress_bg = 0;
    compress_ed = BufferSize-1;
    compress_size = BufferSize+1;
    compress_data = new bam_block*[compress_size];
    for (int i=compress_bg;i<=compress_ed;i++) compress_data[i] = new bam_block;


    consumer_bg=1;
    consumer_ed=0;
    consumer_size=BufferSize+5;
    consumer_data=new bam_block*[consumer_size];

    compressThread = threadNumber;

    wait_num=0;
}

/*
 * 获取一个空白的内存块
 * 输入：无
 * 输出：bam_block* : 指向该内存块的指针
 */
bam_block* BamCompress::getEmpty(){
    mtx_compress.lock();
    while ((compress_ed+1)%compress_size == compress_bg){
        mtx_compress.unlock();
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
        mtx_compress.lock();
    }
    int num = compress_bg;
    compress_bg=(compress_bg+1)%compress_size;
    mtx_compress.unlock();
    return compress_data[num];
}
/*
 * 按照顺序插入输出队列
 * 输入：  bam_block* : 解压完成的数据
 *        int ：读入的顺序编号
 * 输出：  无
 */

void BamCompress::inputUnCompressData(bam_block* data,int block_num) {
//    std::this_thread::sleep_for(std::chrono::milliseconds(blockNum-block_num));
    mtx_input.lock();
    while (block_num!=blockNum) {
        wait_num+=1;
        mtx_input.unlock();
        std::this_thread::sleep_for(std::chrono::nanoseconds(block_num-blockNum)/2);
        mtx_input.lock();
    }

    consumer_data[(consumer_ed + 1) % consumer_size] = data;
    consumer_ed = (consumer_ed + 1) % consumer_size;
    blockNum++;
    mtx_input.unlock();
//    printf("block Num is %d\n",blockNum);
}


//void BamCompress::inputUnCompressData(bam_block* data,int block_num) {
//    std::this_thread::sleep_for(std::chrono::milliseconds(blockNum-block_num));
//    mtx_input.lock();
//    while (block_num != blockNum) {
//        wait_num+=1;
//        mtx_input.unlock();
//        std::this_thread::sleep_for(std::chrono::milliseconds(blockNum-block_num));
//        mtx_input.lock();
//    }
//
//    consumer_data[(consumer_ed + 1) % consumer_size] = data;
//    consumer_ed = (consumer_ed + 1) % consumer_size;
//    blockNum++;
//    mtx_input.unlock();
////    printf("block Num is %d\n",blockNum);
//}


///*
// * 按照顺序尝试插入输出队列
// * 输入：  bam_block* : 解压完成的数据
// *        int ：读入的顺序编号
// * 输出：  无
// */
//bool BamCompress::tryinputUnCompressData(bam_block* data,int block_num) {
//    if  (block_num != blockNum) {
//        wait_num+=1;
//        return false;
//    }
//    consumer_data[(consumer_ed + 1) % consumer_size] = data;
//    consumer_ed = (consumer_ed + 1) % consumer_size;
//    blockNum++;
//    return true;
////    printf("block Num is %d\n",blockNum);
//}
/*
 * 按照顺序获取解压完成的数据
 * 输入：无
 * 输出：bam_block* 解压完成的数据
 */
bam_block* BamCompress::getUnCompressData(){
    while ((consumer_ed+1)%consumer_size == consumer_bg){
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
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
void BamCompress::backEmpty(bam_block* data){
    compress_data[(compress_ed+1)%compress_size] = data;
    compress_ed = (compress_ed+1)%compress_size;
}
