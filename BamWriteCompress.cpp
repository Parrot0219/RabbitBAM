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
    need_compress_bg=1;
    need_compress_ed=0;
    need_compress_size = BufferSize+5;
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
//    mtx_compress.lock();
    while ((compress_ed+1)%compress_size == compress_bg){
        mtx_compress.unlock();
        std::this_thread::sleep_for(std::chrono::nanoseconds(5));
        mtx_compress.lock();
    }
    int num = compress_bg;
    compress_bg=(compress_bg+1)%compress_size;
//    mtx_compress.unlock();
    return compress_data[num];
}

void BamWriteCompress::inputUnCompressData(bam_write_block* data){
//    mtx_need_compress.lock();
    data->block_num=blockInputNum++;
//    printf("Input Data Block Num : %d\n",data->block_num);
    need_compress_data[(need_compress_ed+1)%need_compress_size]=data;
    need_compress_ed=(need_compress_ed+1)%need_compress_size;
//    if (blockInputNum < 4) {
//        printf("This Time need compress bg : %d need compress ed : %d\n",need_compress_bg,need_compress_ed);
//        printf("Input One Over need compress bg block Num : %d\n",need_compress_data[need_compress_bg]->block_num);
//    }
    if (WriteDeBug){
        printf("In InputUnCompressData: need_compress_bg : %d\n",need_compress_bg);
        printf("In InputUnCompressData: need_compress_ed : %d\n",need_compress_ed);
        printf("Input Un CompressData This Thread ID : %d\n",std::this_thread::get_id());
    }

//    mtx_need_compress.unlock();
}
bam_write_block* BamWriteCompress::getUnCompressData(){
    mtx_need_compress.lock();
    while ((need_compress_ed+1)%need_compress_size == need_compress_bg){
        if (WriteDeBug){
            printf("Get Un Compress Data This Thread ID : %d\n",std::this_thread::get_id());
            printf("need compress bg : %d\n"
                   "need compress ed : %d\n",need_compress_bg,need_compress_ed);
        }
        mtx_need_compress.unlock();
//        printf("Compress Sleep Start\n");
        std::this_thread::sleep_for(std::chrono::nanoseconds(1));
//        printf("Compress Sleep End\n");
        if (isWriteComplete && (need_compress_ed+1)%need_compress_size==need_compress_bg){
            return nullptr;
        }
        mtx_need_compress.lock();

    }
//    printf("need compress bg : %d\n"
//           "need compress ed : %d\n",need_compress_bg,need_compress_ed);
    int num=need_compress_bg;
//    printf("this need compress data block num : %d\n",need_compress_data[num]->block_num);
    need_compress_bg=(need_compress_bg+1)%need_compress_size;
    mtx_need_compress.unlock();

    return need_compress_data[num];
}


/*
 * 按照顺序插入输出队列
 * 输入：  bam_block* : 压缩完成的数据
 *        int ：读入的顺序编号
 * 输出：  无
 */

void BamWriteCompress::inputCompressData(bam_write_block* data){
//    printf("When Input Compress Data the data block num : %d\n",data->block_num);
    while (data->block_num != blockNum.load(std::memory_order_acq_rel)) {
        wait_num+=1;
//        printf("Wait Sleep Start\n");
        std::this_thread::sleep_for(std::chrono::nanoseconds(data->block_num-blockNum)/8);
//        printf("Wait Sleep End\n");
    }
//    printf("comsumer bg : %d\n"
//           "consumer ed : %d\n",consumer_bg,consumer_ed);
    consumer_data[(consumer_ed + 1) % consumer_size] = data;
    consumer_ed = (consumer_ed + 1) % consumer_size;
//    printf("Input Compress Data consumer bg : %d\n"
//           "Input Compress Data consumer ed : %d\n",consumer_bg,consumer_ed);
    blockNum.store(blockNum.load(std::memory_order_acq_rel)+1,std::memory_order_acq_rel);
}





/*
 * 按照顺序获取解压完成的数据
 * 输入：无
 * 输出：bam_block* 解压完成的数据
 */
bam_write_block* BamWriteCompress::getCompressData(){
    while ((consumer_ed+1)%consumer_size == consumer_bg){
//        printf("Write Sleep Start\n");
        std::this_thread::sleep_for(std::chrono::nanoseconds(1));
//        printf("Write Sleep end\n");
        if (compressThread==0 && (consumer_ed+1)%consumer_size == consumer_bg) return nullptr;
    }
    int num = consumer_bg;
//    is_ok[compress_bg]=false;
//    printf("Get Compress Data consumer bg : %d\n"
//           "Get Compress Data consumer ed : %d\n",consumer_bg,consumer_ed);
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


